function [d,adj_model,transform,opts_used]=psg_geo_pwaffine(ref,adj,opts)
% [d,adj_model,transform,opts_used]=psg_geo_projective(ref,adj,opts) finds a 
% piecewise affine model with standardized input and output variables
%
%  See psg_piecwise_notes.doc for details on algorithm for known acut, vcut
%  * Do an exhaustive search of vcut trial directions and acut values
%  * Use the empirical minimum as the starting point for a nonlinear optimization
%  * For the nonlinear optimization, do an uncontstrained search:
%      Consider the tangent plane to the trial vcut, and construct orthonormal coordinates
%      Search in these coordinates: direction to the tangent plane is the
%      new vcut, and distance from the origin to this plane is multiplier for acut
%  * The heuristic is that the optimum is close to the trial vcut, and it is
%    unlikely that one will need to search even close to the equator
%    
% ref: reference coordinates, size=[npts,dim_y]
% adj: coordinates of dataset to be adjusted, size=[npts,dim_x]
% opts: options, can be omitted or empty
%    opts.if_display: 1 to display messages from fmincon
%
% d: goodness of fit, mean squared deviation of adj_model from res, divided by mean squared dev of ref
% adj_model: coordinates of adj, after transformation
% transform is:
%   transform.b: scalar, equal to 1 (scale absorbed in Tpos, Tneg)
%   transform.T: stack of matrices, size [dim_x max(dim_x,dim_y) 2], use (:,:,1) when vcut*adj'>=acut, (:,:,2) when vcut*adj'<=acut
%   transform.c: stack of offsets, size [2 max(dim_x,dim_y) 2], use (:,1) when vcut*adj'>=acut, (:,2) when vcut*adj'<=a
%   transform.vcut: unit vector, as a row of size [1 dim_x], orthog to cut plane
%   transform.acut: cutpoint
%
% opts_used: options used, also auxiliary outputs including inital values found by the exhuastive search step
%   and the best-fit d by standard affine transformation
%
%   See also: PSG_GEOMODELS_TEST, PSG_GEOMODELS_DEFINE, PSG_GEO_GENERAL, FMINCON, REGRESS, PSG_PWAFFINE_APPLY,
%   HSPHERE_SAMPLE, EXTORTHB.
%
if (nargin<=2) opts=struct; end
opts=filldefault(opts,'method','fmin');
opts=filldefault(opts,'if_display',1);
opts=filldefault(opts,'fmin_opts',optimset('fminsearch'));
opts=filldefault(opts,'hsphere_method','axes_and_orthants');
opts=filldefault(opts,'hsphere_nsamps',32);
opts=filldefault(opts,'n_cuts_init',7);
opts=filldefault(opts,'tol_cut',10^-7);
opts=filldefault(opts,'if_nearinit',0); %1 to test near the initialization point
opts=filldefault(opts,'if_keep_inits',0); %1 to keep all values of d during initialization step
%
if opts.if_display==0 %turn off display in fminsearch
    opts.fmin_opts=optimset(opts.fmin_opts,'Display','off');
end
opts_used=opts;
%
dim_x=size(adj,2);
dim_y=size(ref,2);
%
npts=size(adj,1);
if dim_y<dim_x
    ref=[ref,zeros(npts,dim_x-dim_y)];
end
dim_xy=max(dim_x,dim_y);
%
n_cuts_init=opts.n_cuts_init; %this is a parameter for initialization
n_pw=2;
%
% find the piecewise affine transformation
%
% determine initial guesses for vcut
%
opts_hsphere=struct;
opts_hsphere.method=opts.hsphere_method;
opts_hsphere.nsamps=opts.hsphere_nsamps;
opts_hsphere.if_hemisphere=1; %no need to sample opposite end of a radius
[vcut_inits,ou_hsphere]=hsphere_sample(dim_x,opts_hsphere);
n_dirs_init=size(vcut_inits,1);
opts_used.vcut_inits=vcut_inits;
%retrieve any new fields from ou_hsphere
fns=fieldnames(ou_hsphere);
for ifn=1:length(fns)
    fn=fns{ifn};
    opts_used=filldefault(opts_used,cat(2,'hsphere_',fn),ou_hsphere.(fn));
end
%
opts_used.acut_inits=zeros(n_dirs_init,n_cuts_init);
opts_used.d_inits=zeros(n_dirs_init,n_cuts_init);
d_min=Inf;
for i_dir=1:n_dirs_init
    vcut_init=vcut_inits(i_dir,:);
    acut_range=adj*vcut_init';
    cut_vals=linspace(min(acut_range),max(acut_range),n_cuts_init);
    for i_cut=1:n_cuts_init
        acut_init=cut_vals(i_cut);
        [d,transform_init,u_init,ou_va]=psg_geo_pwaffine_va(ref,adj,vcut_init,acut_init,opts);
        opts_used.acut_inits(i_dir,i_cut)=acut_init;
        opts_used.d_inits(i_dir,i_cut)=d;
        if (d<d_min)
            d_min=d;
            %keep best values from initialization           
            u_init_min=u_init;
            acut_init_min=acut_init;
            opts_used.u_init_min=u_init_min;
            opts_used.acut_init_min=acut_init_min;
            opts_used.vcut_init_min=vcut_init;
            opts_used.d_init_min=d;
            opts_used.transform_init_min=transform_init;
            opts_used=filldefault(opts_used,'if_orth',ou_va.if_orth);
        end
    end       
end
if opts.if_keep_inits==0
    opts_used=rmfield(opts_used,'d_inits');
    opts_used=rmfield(opts_used,'vcut_inits');
    opts_used=rmfield(opts_used,'acut_inits');
end
transform=opts_used.transform_init_min;
%
[d_affine,adj_model_affine,transform_affine]=psg_geo_affine(ref,adj,opts);
opts_used.d_affine=d_affine;
opts_used.transform_affine=transform_affine; % for reference
%
%optionally test some nearby points
%
if opts.if_nearinit
    params_nearinit=[zeros(1,dim_x);eye(dim_x)*0.1;-eye(dim_x)*0.1];
    opts_used.params_nearinit=params_nearinit;
    opts_used.d_nearinit=zeros(1,size(params_nearinit,1));
    opts_used.acut_nearinit=zeros(1,size(params_nearinit,1));
    opts_used.vcut_nearinit=zeros(size(params_nearinit,1),dim_x);
    for itest=1:size(params_nearinit,1)
        [opts_used.d_nearinit(itest),transform_nearinit]=psg_geo_pwaffine_obj(params_nearinit(itest,:),ref,adj,u_init_min,acut_init_min);
        opts_used.acut_nearinit(itest)=transform_nearinit.acut;
        opts_used.vcut_nearinit(itest,:)=transform_nearinit.vcut;
    end
end
%
%do a minimization
%
fminsearch_opts=opts.fmin_opts;
if (opts.if_display==0)
    fminsearch_opts=optimset(fminsearch_opts,'Display','off');
end
[p,fval,exitflag,output]=fminsearch(@(p) psg_geo_pwaffine_obj(p,ref,adj,u_init_min,acut_init_min,opts),zeros(1,dim_x),fminsearch_opts);
opts_used.fmin.ssq=fval;
opts_used.fmin.exitflag=exitflag;
opts_used.fmin.fminsearch_opts=fminsearch_opts;
opts_used.fmin.output=output;
[opts_used.d_final,transform]=psg_geo_pwaffine_obj(p,ref,adj,u_init_min,acut_init_min);
%
%final calculation with the best transform
adj_model=psg_pwaffine_apply(transform,adj);
d_num=sum(sum((ref-adj_model).^2,1));
d_den=sum(sum((ref-repmat(mean(ref,1),npts,1)).^2,1));
d=d_num/d_den;
return
end

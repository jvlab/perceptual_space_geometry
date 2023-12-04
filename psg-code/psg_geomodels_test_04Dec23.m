%psg_geomodels_test: test: test fitting geometric models
%
% differs from psg_procrustes_regr_test:
%    * includes a knitted built-in dataset
%    * procrustes with and without scaling
%    * regression model is called 'affine'
%    * for shuffling, reference dataset is permuted rather than adj dataset
%    * shuffling done for all nested models, including the zero model
%    * psg_geomodels_define defines standard model defaults and nestings
%    * results from procrustes, affine, and projective models computed from a master routine
%    * d (goodness of fit) also computed with a denominator equal to variance of the original data
%       (otherwise, surrogates that inflate the variance will artifactually decrease the d)
%    * only one method is used for projective ('fmin')
%    * no fitting of projective to residuals of Procrustes
%    * frozen random numbers for each model and permutations, in case minimization is stochastic
%
% plan is to wrap procrustes in psg_geomodel_procrustes, 
% separate out psg_geo_procrustes (and clean up c), psg_geo_affine, psg_geo_projective
% use this as a platform to add piecewise models
%
%   See also:  PROCRUSTES, REGRESS, PERSP_XFORM_FIND, PSG_GET_COORDSETS, PSG_PROCRUSTES_REGR_TEST,
%     PSG_PROCRUSTES_REGR_DEMO.
%
model_types_def=psg_geomodels_define();
model_types=model_types_def.model_types;
%
nmodels=length(model_types);
%
if ~exist('ref_dim') ref_dim=3; end
if ~exist('adj_dim') adj_dim=3; end
%
if ~exist('if_center') if_center=1; end
if ~exist('nshuff') nshuff=100; end
if ~exist('nshuff_regr') nshuff_regr=100; end
if ~exist('if_cycle') if_cycle=1; end %option for projective fit, but ignored as long as 'fmin' is used in persp_xform_find
if ~exist('opts_read') opts_read=struct; end
if ~exist('opts_rays') opts_rays=struct; end
if ~exist('opts_qpred') opts_qpred=struct; end
if ~exist('opts_persp') opts_persp=struct; end
%
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
%
if_builtin=getinp('use built-in datasets (1: bgca3pt MC and BL; 2: bc55pooled MC','d',[0 2]);
switch if_builtin
    case 1
        if ~exist('ref_file') ref_file='./psg_data/bgca3pt_coords_MC_sess01_10.mat'; end
        if ~exist('adj_file') adj_file='./psg_data/bgca3pt_coords_BL_sess01_10.mat'; end
        ref_dim_max=7;
        adj_dim_max=7;
    case 2
        ref_file='./psg_data/bc55pooled_coords_MC.mat';
        adj_file='./psg_data/bc55pooled9.mat';
        nsets=2;
        opts_read=filldefault(opts_read,'if_log',1);
        opts_read.input_type=[1 2];
        opts_read.data_fullname_def=ref_file;
        opts_read.setup_fullname_def=adj_file;
        [sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred,nsets);
        ref_dim_max=max(sets{1}.dim_list);
        adj_dim_max=max(sets{2}.dim_list);
    otherwise       
        disp('dataset 1 will be reference, dataset 2 will be adjusted to fit.');
        nsets=2;
        opts_read=filldefault(opts_read,'if_log',1);
        [sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred,nsets);
        ref_file=sets{1}.label_long;
        adj_file=sets{2}.label_long;
        ref_dim_max=max(sets{1}.dim_list);
        adj_dim_max=max(sets{2}.dim_list);
end
disp(sprintf('reference dataset: %s',ref_file));
disp(sprintf('dataset to be adjusted: %s',adj_file));
ref_dim=getinp('reference dataset dimension to use','d',[1 ref_dim_max],ref_dim);
adj_dim=getinp('adjusted  dataset dimension to use','d',[1 adj_dim_max],adj_dim);
switch if_builtin
    case 1
        ref=getfield(load(ref_file),sprintf('dim%1.0f',ref_dim));
        adj=getfield(load(adj_file),sprintf('dim%1.0f',adj_dim));
    case 2
        ref=ds{1}{ref_dim};
        adj=ds{2}{adj_dim};
    otherwise
        ref=ds{1}{ref_dim};
        adj=ds{2}{adj_dim};
end
npts=size(ref,1);
disp(sprintf('number of data points: %3.0f',npts));
%model-fitting options
if_center=getinp('1 to center the data','d',[0 1],if_center);
%
tstring=sprintf(' ref %s dim %1.0f adj %s dim %1.0f center %1.0f',...
    ref_file,ref_dim,adj_file,adj_dim,if_center);
%
npts=size(ref,1);
if if_center
    ref=ref-repmat(mean(ref,1),npts,1);
    adj=adj-repmat(mean(adj,1),npts,1);
end
%
d=zeros(nmodels,1);
d_calc_types={'den: surrogate','den:  original'};
%    
d_shuff=zeros(nmodels,nshuff,nmodels-1,length(d_calc_types)); %d1: model, d2: shuffle, d3: nested model, d4: normalization type
transforms=cell(nmodels,1);
adj_model=cell(nmodels,1);
adj_model_check=cell(nmodels,1);
resids=cell(nmodels,1);
d_check=zeros(nmodels,1);
opts_model_used=cell(nmodels,1);
opts_model_shuff_used=cell(nmodels,nshuff,nmodels-1);
%
d_den=sum(sum((ref-repmat(mean(ref,1),npts,1)).^2,1));
if ref_dim<adj_dim
    ref_aug=[ref,zeros(npts,adj_dim-ref_dim)];
else
    ref_aug=ref;
end
if (if_frozen~=0)
    rng('default');
    if (if_frozen<0)
        rand(1,abs(if_frozen));
    end
else
    rng('shuffle');
end
perms=zeros(nshuff,npts);
for ishuff=1:nshuff
    perms(ishuff,:)=randperm(npts);
end
for imodel=1:nmodels
    %
    if (if_frozen~=0)
        rng('default');
        if (if_frozen<0)
            rand(1,abs(if_frozen));
        end
    else
        rng('shuffle');
    end
    model_type=model_types{imodel};
    disp(' ');
    disp(sprintf('model type: %s',model_type))
    opts_model=model_types_def.(model_type).opts;
    model_class=model_types_def.(model_type).class;
    [d(imodel),adj_model{imodel},transforms{imodel},opts_model_used{imodel}]=psg_geo_general(ref,adj,model_class,opts_model);
    %
    %verify the model   
    %
    switch model_class
        case {'mean','procrustes','affine'}
            adj_model_check{imodel}=transforms{imodel}.b*adj*transforms{imodel}.T+repmat(transforms{imodel}.c,npts,1);
        case 'projective'
            %note change of variable names between psg_geo_projective and persp_apply
            adj_model_check{imodel}=persp_apply(transforms{imodel}.T,transforms{imodel}.c,transforms{imodel}.p,adj);
    end
    %calculate residuals and verify d
    resids{imodel}=ref_aug-adj_model{imodel};
    d_check(imodel)=sum(resids{imodel}(:).^2)/d_den;
    disp(sprintf('d: %12.7f   d_check: %12.7f   diff: %12.7f',...
        d(imodel),d_check(imodel),d(imodel)-d_check(imodel)));
    %check the reconstitution
    disp('transform parameters:')
    if isfield(transforms{imodel},'b')
        disp('b');
        disp(transforms{imodel}.b);       
    end
    %
    if isfield(transforms{imodel},'T')
        disp('T');
        disp(transforms{imodel}.T);       
    end
    if isfield(transforms{imodel},'c')
        disp('c');
        disp(transforms{imodel}.c(1,:));
    end
    if isfield(transforms{imodel},'p')
        disp('p transpose');
        disp(transforms{imodel}.p');
    end
    %
    %shuffles
    %
    if nshuff>0 & isempty(opts_model_used{imodel}.warnings)
        nested_types=model_types_def.(model_type).nested;
        for inest=1:length(nested_types)
            nested_type=nested_types{inest};
            nest_ptr=strmatch(nested_type,model_types,'exact');
            disp(sprintf(' doing shuffles with residuals from nested model %2.0f (%s)',nest_ptr,nested_type));
            for ishuff=1:nshuff
                %check this next line -- possibly an issue that the denominator of d has changed?
                shuffled=adj_model{nest_ptr}+(ref_aug(perms(ishuff,:),:)-adj_model{nest_ptr}(perms(ishuff,:),:));
                %d_shuff_orig is calculated with surrogate denominator; we also want to calculate it with original denominator
                %also, turn off fmin display for shuffles
                [d_shuff_orig,adj_model_shuff,transform_shuffle,opts_model_shuff_used{imodel,ishuff,inest}]=...
                    psg_geo_general(shuffled,adj,model_class,setfield(opts_model,'if_display',0));
                resids_shuff=shuffled-adj_model_shuff; %deviation of model from shuffle surrogate
                d_shuff(imodel,ishuff,inest,1)=d_shuff_orig;
                d_shuff(imodel,ishuff,inest,2)=sum(resids_shuff(:).^2)/d_den;
            end
            for id_calc_type=1:length(d_calc_types)
                disp(sprintf('d is >= d_shuff for %5.0f of %5.0f shuffles (%s), d_shuffles: range [%8.5f %8.5f], mean %8.5f s.d. %8.5f',...
                    sum(double(d(imodel)>=d_shuff(imodel,:,inest,id_calc_type))),nshuff,d_calc_types{id_calc_type},...
                    min(d_shuff(imodel,:,inest,id_calc_type)),max(d_shuff(imodel,:,inest,id_calc_type)),...
                    mean(d_shuff(imodel,:,inest,id_calc_type)),std(d_shuff(imodel,:,inest,id_calc_type))));
            end
        end %inest
    end %if shuffled
end

function model_types_def=psg_geomodels_define();
%model_types_def=psg_geomodels_define() sets up the definitions of geometric model types
%
% model_types_def: a structure that defines the models and their hierarchical relationships
%
%   See also: PSG_GEOMODELS_TEST, PSG_GEO_GENERAL.
%
model_types_def=struct;
model_types_def.model_types={'mean','procrustes_noscale','procrustes_scale','affine_nooffset','affine_offset','projective'};
%
model_types_def.mean.opts=struct;
model_types_def.mean.nested={};
%
model_types_def.procrustes_noscale.opts.if_scale=0;
model_types_def.procrustes_noscale.nested={'mean'};
%
model_types_def.procrustes_scale.opts.if_scale=1;
model_types_def.procrustes_scale.nested={'mean','procrustes_noscale'};
%
model_types_def.affine_nooffset.opts.if_offset=0;
model_types_def.affine_nooffset.nested={'mean','procrustes_noscale','procrustes_scale'};
%
model_types_def.affine_offset.opts.if_offset=1;
model_types_def.affine_offset.nested={'mean','procrustes_noscale','procrustes_scale'};
%
model_types_def.projective.opts.method='fmin';
model_types_def.projective.nested={'mean','procrustes_noscale','procrustes_scale','affine_offset'};
%
mnames=model_types_def.model_types;
for im=1:length(mnames)
    mname=mnames{im};
    mclass=mname(1:-1+min(find(cat(2,mname,'_')=='_')));
    model_types_def.(mname).class=mclass;
end
return
end

function [d,adj_model,transform,opts_used]=psg_geo_general(ref,adj,model_class,opts)
% [d,adj_model,transform,opts_used]=psg_geo_general(ref,adj,model_class,opts) finds a general model
% (procrustes, affine, etc) with standardized input and output variables
%
% ref: reference coordinates, size=[npts,ref_dim]
% adj: coordinates of dataset to be adjusted, size=[npts,adj_dim]
% model_class: 'procrustes','affine', or any of the other values in model_types_def.(model).class
% opts: options, can be omitted or empty
%     opts.if_scale: 1 to allow scaling (default: 0)
%
% d: goodness of fit, mean squared deviation of adj_model from res, divided by mean squared dev of ref
% adj_model: coordinates of adj, after transformation
% transform: model-specific, but the following are typical
%   transform.b: scalar, will be 1 if opts.if_scale=0
%   transform.T: matrix of size [adj_dim max(ref_dim,adj_dim)]
%   transform.c: offset, row of size [max(ref_dim, adj_dim) 1], will be zeros if ref and adj are centered
%   adj_model=transform.b*adj*transform.T+repmat(transform.c,npts,1)
% opts_used: options used
%   opts_used.warnings: warnings
%
%   See also: PROCRUSTES, PSG_GEOMODELS_TEST, PSG_GEOMODELS_DEFINE, PSG_GEO_PROCRUSTES, PSG_GEO_AFFINE, PSG_GEO_PROJECTIVE.
%
if (nargin<=3) opts=struct; end
adj_dim=size(adj,2);
ref_dim=size(ref,2);
%
npts=size(adj,1);
%
switch model_class
    case 'procrustes'
        [d,adj_model,transform,opts_used]=psg_geo_procrustes(ref,adj,opts);
    case 'affine'
        [d,adj_model,transform,opts_used]=psg_geo_affine(ref,adj,opts);
    case 'projective'
        [d,adj_model,transform,opts_used]=psg_geo_projective(ref,adj,opts);
    case 'mean' %model is mean of the reference dataset
        d=1;
        adj_model=repmat(mean(ref,1),npts,1);
        if (ref_dim<adj_dim)
            adj_model=[adj_model,zeros(npts,adj_dim-ref_dim)];
        end
        opts_used=opts;
        transform=struct;
        transform.b=1;
        transform.T=zeros(adj_dim,max(ref_dim,adj_dim));
        transform.c=adj_model(1,:);
    otherwise
        wstring=sprintf('model class %s not recognized.',model_class);
        d=1;
        adj_model=zeros(npts,max(ref_dim,adj_dim));
        transform=struct;
        opts_used=opts;
        warning(wstring);
        opts_used.warnings=wstring;
end
if ~isfield(opts_used,'warnings')
    opts_used.warnings=[];
end
return
end

function [d,adj_model,transform,opts_used]=psg_geo_procrustes(ref,adj,opts)
% [d,adj_model,transform,opts_used]=psg_geo_procrustes(ref,adj,opts) finds a procrustes model
% with standardized input and output variables
%
% ref: reference coordinates, size=[npts,ref_dim]
% adj: coordinates of dataset to be adjusted, size=[npts,adj_dim]
% opts: options, can be omitted or empty
%     opts.if_scale: 1 to allow scaling (default: 0)
%
% d: goodness of fit, mean squared deviation of adj_model from res, divided by mean squared dev of ref
% adj_model: coordinates of adj, after transformation
% transform:
%   transform.b: scalar, will be 1 if opts.if_scale=0
%   transform.T: matrix of size [adj_dim max(ref_dim,adj_dim)]
%   transform.c: offset, row of size [max(ref_dim, adj_dim) 1], will be zeros if ref and adj are centered
%   adj_model=transform.b*adj*transform.T+repmat(transform.c,npts,1)
% opts_used: options used
%
%   See also: PROCRUSTES, PSG_GEO_AFFINE, PSG_GEOMODELS_TEST, PSG_GEOMODELS_DEFINE, PSG_GEO_GENERAL.
%
if (nargin<=2) opts=struct; end
opts=filldefault(opts,'if_scale',0);
adj_dim=size(adj,2);
ref_dim=size(ref,2);
opts_used=opts;
%
npts=size(adj,1);
if ref_dim<adj_dim
    ref=[ref,zeros(npts,adj_dim-ref_dim)];
end
%
if opts.if_scale==0
    Scaling=false;
else
    Scaling=true;
end
%
% Procrustes
%
[d,adj_model,transform]=procrustes(ref,adj,'Scaling',Scaling);
transform.c=transform.c(1,:); %clean up redundant rows
return
end

function [d,adj_model,transform,opts_used]=psg_geo_affine(ref,adj,opts)
% [d,adj_model,transform,opts_used]=psg_geo_affine(ref,adj,opts) finds an affine model
% using regression, with standardized input and output variables
%
% ref: reference coordinates, size=[npts,ref_dim]
% adj: coordinates of dataset to be adjusted, size=[npts,adj_dim]
% opts: options, can be omitted or empty
%   opts.if_offset: 1 to allow offset (default)
%
% d: goodness of fit, mean squared deviation of adj_model from res, divided by mean squared dev of ref
% adj_model: coordinates of adj, after transformation
% transform:
%   transform.b: scalar, here always 1 since scale factor is absorbed into T 
%   transform.T: matrix of size [adj_dim max(ref_dim,adj_dim)]
%   transform.c: offset, row of size [max(ref_dim, adj_dim) 1] 
%   adj_model=transform.b*adj*transform.T+repmat(transform.c,npts,1)
% opts_used: options used  
%
% if ref dimension is less than adj dimension, then ref will be padded with
% columns of zeros
%
%   See also:  REGRESS, PSG_GEO_PROCRUSTES, PSG_GEOMODELS_TEST, PSG_GEOMODELS_DEFINE, PSG_GEO_GENERAL.
%
if (nargin<=2) opts=struct; end
opts=filldefault(opts,'if_offset',1);
adj_dim=size(adj,2);
ref_dim=size(ref,2);
opts_used=opts;
%
npts=size(adj,1);
if ref_dim<adj_dim
    ref=[ref,zeros(npts,adj_dim-ref_dim)];
end
%
% regression (affine) fit
%
t_regr=zeros(opts.if_offset+adj_dim,ref_dim);
if (opts.if_offset==1)
    adj_aug=[ones(npts,1),adj];
else
    adj_aug=adj;
end
for iref_dim=1:max(ref_dim,adj_dim) 
    t_regr(:,iref_dim)=regress(ref(:,iref_dim),adj_aug);
end
adj_model=adj_aug*t_regr; %reconstitute based on regression
transform.T=t_regr(1+opts.if_offset:end,:);
if opts.if_offset
    transform.c=t_regr(1,:);
else
    transform.c=zeros(1,max(ref_dim,adj_dim));
end
transform.b=1; %for consistency with Procrustes
d_num=sum(sum((ref-adj_model).^2,1));
d_den=sum(sum((ref-repmat(mean(ref,1),npts,1)).^2,1));
d=d_num/d_den;
return
end

function [d,adj_model,transform,opts_used]=psg_geo_projective(ref,adj,opts)
% [d,adj_model,transform,opts_used]=psg_geo_projective(ref,adj,opts) finds a 
% projective model (perspective transformation)
% with standardized input and output variables
%
% ref: reference coordinates, size=[npts,ref_dim]
% adj: coordinates of dataset to be adjusted, size=[npts,adj_dim]
% opts: options, can be omitted or empty
%     opts.if_display: 1 to display messages from fminsearch
%
% d: goodness of fit, mean squared deviation of adj_model from res, divided by mean squared dev of ref
% adj_model: coordinates of adj, after transformation
% transform:
%   transform.b: scalar, equal to 1 (scale absorbed in T)
%   transform.T: matrix of size [adj_dim max(ref_dim,adj_dim)]
%   transform.c: offset, row of size [max(ref_dim, adj_dim) 1]
%   transform.p: column of size adj_dim
% opts_used: options used
%
%   See also: PSG_GEOMODELS_TEST, PSG_GEOMODELS_DEFINE, PSG_GEO_GENERAL, PERSP_XFORM_FIND, PERSP_APPLY.
%
if (nargin<=2) opts=struct; end
opts=filldefault(opts,'method','fmin');
opts=filldefault(opts,'if_cycle',1); %ignored unless opts.method is set to 'oneshot';
opts=filldefault(opts,'if_display',1);
opts=filldefault(opts,'fmin_opts',optimset('fminsearch'));
%
if opts.if_display==0 %turn off display in fminsearch
    opts.fmin_opts=optimset(opts.fmin_opts,'Display','off');
end
opts_used=opts;
%
adj_dim=size(adj,2);
ref_dim=size(ref,2);
%
npts=size(adj,1);
if ref_dim<adj_dim
    ref=[ref,zeros(npts,adj_dim-ref_dim)];
end
%
% find the perspective transformation
%
[persp,adj_model,ou_persp]=persp_xform_find(adj,ref,opts);
%retrieve any new fields from ou_persp
fns=fieldnames(ou_persp);
for ifn=1:length(fns)
    fn=fns{ifn};
    opts_used=filldefault(opts_used,fn,ou_persp.(fn));
end
%
% intepretration of persp (variable names changed here for compatibility with psg_geo_procrustes and
% psg_geo_affine: 
%
%  persp=[T | p]
%        [-----]
%        [c | 1]
%
% Each row of input (adj) is considered as a homogeneous vector with an augmented coordinate 1 at the end
% This matrix is then post-multiplied by persp, yielding a homogeneous vector.
% Output (adj_model) are the rows of this matrix, divided by the final element.
%   Notes
%     persp only matters up to homogeneity, but its lower right element is fixed at 1.
%     If p=persp(1:end-1,end) c is zero, this is an affine transformation
%     with offset c=persp(end,1:end-1);
%     if p and c are zero, this is a linear transformation by T.
transform.b=1;
transform.T=persp(1:end-1,1:end-1);
transform.c=persp(end,1:end-1);
transform.p=persp(1:end-1,end);
%
d_num=sum(sum((ref-adj_model).^2,1));
d_den=sum(sum((ref-repmat(mean(ref,1),npts,1)).^2,1));
d=d_num/d_den;
%
return
end

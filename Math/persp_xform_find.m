function [persp,y_fit,opts_used]=persp_xform_find(x,y,opts)
% [persp,y_fit,opts_used]=persp_xform_find(x,y,opts) finds a perspective
% transformation from an overdetermined set of corresponding points  in an arbitrary number of dimensions
%
% Method: Estimating Projective Transformation Matrix (Collineation, Homography)
% Zhengyou Zhang
% November 1993; Updated May 29, 2010
% Microsoft Research Techical Report MSR-TR-2010-63
% This uses Method 2 of the above reference, but the data are in rows
%    Note that this is appropriate for an over-determined set, or when viewpoint is inside of data
%
% x: data points to be transformed, size [npts,nd], nd is number of dimensions
% y: target data points, size [npts,nd], nd is number of dimensions
% opts: options
%   opts.method='fmin': use method based on fminsearch, persp_ssq_dif
%   opts.method='oneshot': use Zhang method above (default) with padding if input and output dimensions are unequal
%    if opts.method='oneshot' or 'best':
%      opts.if_cycle: whether to cycle through assignments of which data point
%      in x is treated as the 'first' point and finds the one with the best-fit 
%      in the mean-squared sense, and opts_used.best_point indicates which point yielded
%      the best fit.
%   if opts.method='fmin' or 'best'
%      opts.if_multi: 0 (default) for a single start, 1 to start along each axis, pos and neg
%      opts.fmin_startfactor: multiplier for initialization for distance away from largest value on the axis, defaults to 2
%   opts.method='best': both methods are used, and the results of the best method is returned.
%      opts_used.method_best is the method used
%      opts_used.[fmin,oneshot].[yfit,ssq] are the results of each method
%   opts.fmin_opts: options structure for fminsearch
%
% persp: matrix of size [nd+1 nd+1]
%    persp*[x ones(npts,1)] are the homogeneous coords cooresponding to y
% y_fit: size [npts nd], the non-homogeneous mapping of x via persp, x*persp
% opts_used: options used
%    opts_used.ssq: sum of squares of deviations of y_fit from y
%
%  01Mar26: calls two methods as subroutines, and 'best' option chooses the best
%  11Jun26: fix bug regarding 'best', was only calling 'fmin'
%  18Nov26: begin if_multi option for fmin 
%
%   See also: FIND_XFORM_PROJ_TEST, REGRESS, FILLDEFAULT, PERSP_SSQDIF, PERSP_SSQDIF_FIT.
%
if (nargin<=2) opts=struct; end
opts=filldefault(opts,'method','oneshot');
opts=filldefault(opts,'if_cycle',0);
opts=filldefault(opts,'fmin_opts',struct());
opts=filldefault(opts,'if_multi',0);
opts=filldefault(opts,'fmin_startfactor',2);
opts_used=opts;
%
npts=size(x,1);
nd=size(x,2);
if (npts<=nd+1)
    warning(sprintf('perspective fitting is underdetermined.  need npts>=nd+2, nd=%3.0f npts=%3.0f',nd,npts));
end
switch opts.method   
    case 'oneshot'
        [persp,y_fit,opts_used]=persp_xform_find_oneshot(x,y,opts_used);
    case 'fmin'
        [persp,y_fit,opts_used]=persp_xform_find_fmin(x,y,opts_used);
    case 'best'
        z=struct;
        [z.oneshot.persp,z.oneshot.y_fit,z.oneshot.opts_used]=persp_xform_find_oneshot(x,y,opts_used); %11Jun26, changed from persp_xform_find_fmin
        [z.fmin.persp,z.fmin.y_fit,z.fmin.opts_used]=persp_xform_find_fmin(x,y,opts_used);
        if z.fmin.opts_used.ssq<=z.oneshot.opts_used.ssq
            method_best='fmin';
        else
            method_best='oneshot';
        end
        opts_used.method_best=method_best;
        opts_used.oneshot=z.oneshot;
        opts_used.fmin=z.fmin;
        persp=z.(method_best).persp;
        y_fit=z.(method_best).y_fit;
        opts_used.ssq=z.(method_best).opts_used.ssq;
end %opts.method       
return
end

function [persp,y_fit,opts_used]=persp_xform_find_oneshot(xi,yi,opts)
npts=size(xi,1);
ndx=size(xi,2);
ndy=size(yi,2);
if ndx<ndy
    x=[xi,zeros(npts,ndy-ndx)];
else
    x=xi;
end
if ndy<ndx
    y=[yi,zeros(npts,ndx-ndy)];
else
    y=yi;
end
%Zhang method
opts_used=opts;
nd=size(x,2);
%
x_aug_orig=[x,ones(npts,1)];
y_aug_orig=[y,ones(npts,1)];
%
p=zeros(nd+1,nd+1);
if (opts.if_cycle==0)
    ncycle=1;
else
    ncycle=npts;
end
persp=zeros(nd+1,nd+1);
y_fit=zeros(npts,nd,ncycle);
sumsq=Inf;
opts_used.oneshot.sumsq_trial=zeros(1,ncycle);
for icycle=0:ncycle-1
    cyc=mod(icycle+[0:npts-1],npts)+1;
    x_aug=x_aug_orig(cyc,:);
    y_aug=y_aug_orig(cyc,:);    
    %
    %create the "A" matrix 
    %
    A=zeros((nd+1)*npts,(nd+1)^2-1+npts);
    for ipt=1:npts
        Mi=zeros(nd+1,(nd+1)^2);
        for id=1:nd+1
            Mi(id,(nd+1)*(id-1)+[1:(nd+1)])=x_aug(ipt,:);
        end
        Arows=(nd+1)*(ipt-1)+[1:(nd+1)];
        A(Arows,1:(nd+1)^2)=Mi;
        if (ipt>1)
            A(Arows,(nd+1)^2+ipt-1)=y_aug(ipt,:)';
        end  
    end
    b=zeros((nd+1)*npts,1);
    b(1:nd+1)=y_aug(1,:);
    %
    S=regress(b,A);
    persp_trial=reshape(S(1:(nd+1)^2),nd+1,nd+1);
    y_fit_hom=x_aug_orig*persp_trial;
    y_fit_trial=y_fit_hom(:,1:nd)./repmat(y_fit_hom(:,nd+1),1,nd);
    sumsq_trial=sum((y_fit_trial(:)-y(:)).^2);
    opts_used.oneshot.sumsq_trial(icycle+1)=sumsq_trial;
    if (sumsq_trial<sumsq)
        persp=persp_trial;
        sumsq=sumsq_trial;
        opts_used.oneshot.best_point=icycle+1;
        y_fit=y_fit_trial;
        opts_used.ssq=sumsq;
    end
end
y_fit=y_fit(:,1:ndy);
persp=persp([1:ndx end],[1:ndy end]);
return
end

function [persp,y_fit,opts_used]=persp_xform_find_fmin(x,y,opts)
%search via fmin
%strategy is to search over the column vector that forms the denominator. 
%here, it is c -- the column vector which, if zero, makes the transformation affine.  This is the "p" in psg_geo_projective.
opts_used=opts;
npts=size(x,1);
nd=size(x,2);
%
fminsearch_opts=optimset('fminsearch');
fields=fieldnames(opts.fmin_opts);
for ifn=1:length(fields)
    fn=fields{ifn};
    fminsearch_opts.(fn)=opts.fmin_opts.(fn);
end
if opts.if_multi %set up multiple start points?
    x_ranges=[min(x,[],1);max(x,[],1)];
    nz_list=find(any(x_ranges~=0,1)); %coordinates with nonzero values
    x0_list=zeros(nd,1+2*length(nz_list)); %each column is a possible start
    %start corressponds to 1/viewpoint
    for k=1:length(nz_list)
        x_range=x_ranges(:,nz_list(k)); %min and max
        if x_range(1)==0
            x_range(1)=-1;
        end
        if x_range(2)==0
            x_range(2)=+1;
        end
        if sign(x_range(1))==sign(x_range(2)) %one start point near the origin, one far from it
            x_range_sorted=sign(x_range(1))*sort(abs(x_range));
            x0_list(nz_list(k),2*k  )=opts.fmin_startfactor./x_range_sorted(1);
            x0_list(nz_list(k),2*k+1)=1./(opts.fmin_startfactor*x_range_sorted(2));
        else %opposite signs
            x0_list(nz_list(k),2*k  )=1./(opts.fmin_startfactor*x_range(1));
            x0_list(nz_list(k),2*k+1)=1./(opts.fmin_startfactor*x_range(2));
        end
    end
else
    x0_list=zeros(nd,1);
end
%
fval_best=Inf;
opts_used.fmin.x0_list=x0_list;
opts_used.fmin.fval_list=NaN(1,size(x0_list,2));
for istart=1:size(x0_list,2)
    x0=x0_list(:,istart);
    [c_x0,fval_x0,exitflag_x0,output_x0]=fminsearch(@(c) persp_ssqdif_fit(c,x,y),x0,fminsearch_opts);
    opts_used.fmin.fval_list(:,istart)=fval_x0;
    if fval_x0<fval_best
        fval_best=fval_x0;
        %
        x0_best=x0;
        c=c_x0;
        fval=fval_x0;
        exitflag=exitflag_x0;
        output=output_x0;
    end
end
opts_used.fmin.x0_best=x0_best;
opts_used.fmin.ssq=fval;
opts_used.fmin.exitflag=exitflag;
opts_used.fmin.fminsearch_opts=fminsearch_opts;
opts_used.fmin.output=output;
[ssq,y_fit,a,b]=persp_ssqdif_fit(c,x,y);
opts_used.ssq=ssq;
persp=[a c;b 1];
return
end


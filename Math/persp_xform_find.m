function [persp,y_fit,opts_used]=persp_xform_find(x,y,opts)
% [persp,y_fit,opts_used]=persp_xform_find(x,y,opts) finds a perspective
% transformation from an overdetermined set of corresponding points  in an arbitrary number of dimensions
%
% Method: Estimating Projective Transformation Matrix (Collineation, Homography)
% Zhengyou Zhang
% November 1993; Updated May 29, 2010
% Microsoft Research Techical Report MSR-TR-2010-63
% This uses Method 2 of the above reference, but the data are in rows
%    Note that this is appropriate for an over-determined set
%
% x: data points to be transformed, size [npts,nd], nd is number of dimensions
% y: target data points, size [npts,nd], nd is number of dimensions
% opts: options
%   opts.method='fmin': use method based on fminsearch, persp_ssq_dif
%    opts.fmin_opts: options structure for fminsearch
%   opts.method='oneshot': use Zhang method above (default)
%   if opts.method='oneshot':
%    opts.if_cycle: whether to cycle through assignments of which data point
%    in x is treated as the 'first' point and finds the one with the best-fit 
%    in the mean-squared sense, and opts_used.best_point indicates which point yielded
%    the best fit.
%
% persp: matrix of size [nd+1 nd+1]
%    persp*[x ones(npts,1)] are the homogeneous coords cooresponding to y
% y_fit: size [npts nd], the non-homogeneous mapping of x via persp, x*persp
% opts_used: options used
%    opts_used.ssq: sum of squares of deviations of y_fit from y
%
%   See also: FIND_XFORM_PROJ_TEST, REGRESS, FILLDEFAULT, PERSP_SSQDIF, PERSP_SSQDIF_FIT.
%
if (nargin<=2) opts=struct; end
opts=filldefault(opts,'method','oneshot');
opts=filldefault(opts,'if_cycle',0);
opts=filldefault(opts,'fmin_opts',struct());
opts_used=opts;
%
npts=size(x,1);
nd=size(x,2);
if (npts<=nd+1)
    warning(sprintf('perspective fitting is underdetermined.  need npts>=nd+2, nd=%3.0f npts=%3.0f',nd,npts));
end
switch opts.method
    case 'oneshot'
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
    case 'fmin'
        fminsearch_opts=optimset('fminsearch');
        fields=fieldnames(opts.fmin_opts);
        for ifn=1:length(fields)
            fn=fields{ifn};
            fminsearch_opts.(fn)=opts.fmin_opts.(fn);
        end
        [c,fval,exitflag,output]=fminsearch(@(c) persp_ssqdif_fit(c,x,y),zeros(nd,1),fminsearch_opts);
        opts_used.fmin.ssq=fval;
        opts_used.fmin.exitflag=exitflag;
        opts_used.fmin.fminsearch_opts=fminsearch_opts;
        opts_used.fmin.output=output;
        [ssq,y_fit,a,b]=persp_ssqdif_fit(c,x,y);
        opts_used.ssq=ssq;
        persp=[a c;b 1];
end %opts.method       
return

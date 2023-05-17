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
%   opts.if_cycle: whether to cycle through assignments of which data point
%    in x is treated as the 'first' point; this treats all points equally
%    and may result in a better fit
%
% persp: matrix of size [nd+1 nd+1], total sum of squares normalized to 1
%    persp*[x ones(npts,1)] are the homogeneous coords cooresponding to y
% y_fit: size [npts nd], the non-homogeneous mapping of x via persp, x*persp
% opts_used: options used
%
%   See also: FIND_XFORM_PROJ_TEST, REGRESS, FILLDEFAULT.
%
if (nargin<=2) opts=struct; end
opts=filldefault(opts,'if_cycle',0);
opts_used=opts;
%
npts=size(x,1);
nd=size(x,2);
if (npts<=nd+1)
    warning(sprintf('perspective fitting is underdetermined.  need npts>=nd+2, nd=%3.0f npts=%3.0f',nd,npts));
end
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
    p_raw=reshape(S(1:(nd+1)^2),nd+1,nd+1);
    p_norm=p_raw./sqrt(sum(p_raw(:).^2));
    p=p+p_norm;
end
persp=p/ncycle;
y_fit_hom=x_aug_orig*persp;
y_fit=y_fit_hom(:,1:nd)./repmat(y_fit_hom(:,nd+1),1,nd);
return

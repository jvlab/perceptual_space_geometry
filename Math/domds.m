function [eival,eivec]=domds(distmtx,p)
%
% [eival,eivec]=domds(distmtx,p) does a multidimensional scaling with an arbitrary fixing exponent
%
% Note: The eigenvalues ei need to be divided by 2, prior to square root,
% to obtain coords that recapitulate the distances.  This can be ignored if
% the coords only need to yield distances proportional to those in distmtx.
%
%
% distmtx: the distance matrix
% p: the fixing power, should be <=1 (distances are raised to power p prior to mds, i.e., 
%   the centering is performed on a matrix in which distances are raised to the power 2*p
%
% eival: the eigenvalues
% eivec: the eigenvectors
%
if (nargin <=1) p=1; end
npts=size(distmtx,1);
m=distmtx.^(2.*p);
mc=m-repmat(sum(m,2),1,npts)/npts-repmat(sum(m,1),npts,1)/npts+sum(sum(m))/npts/npts;
[ve,ei]=eig(-mc);
eivec=ve;
eival=diag(ei);
return
end

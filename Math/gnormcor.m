function [x,v,d,crt]=gnormcor(cov,npts)
% [x,v,d,crt]=gnormcor(cov,npts) generates a multichannel correlated Gaussian signal
%
% cov is the desired covariance matrix, size(cov)=[n n]
% npts is the number of points to generate
%
% x: the generated signal
% v: eigenvectors of cov
% d: eigenvalues of cov
% crt: the matrix [n n] by which uncorrelated random numbers are transformed
%
%   See also:  GNORMCOR_LC.
%
n=length(cov);
if (min(min(cov == cov'))==0)
   error('covariance must be a square symmetric matrix.');
end
[v,d]=eig(cov);
if (min(min(d)) < 0)
   error('covariance matrix has negative eigenvalues.')'
end
crt=v*sqrt(d)*v';
x=crt*normrnd(0,1,n,npts);
return


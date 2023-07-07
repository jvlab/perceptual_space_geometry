function [x,v,d,crt]=ellipcor(cov,npts)
% [x,v,d,crt]=ellipcor(cov,npts) generates a multichannel signal in which components
% lie on an ellipse, with covariance specified by cov
%
% cov is the desired covariance matrix, size(cov)=[n n]
% npts is the number of points to generate
%
% x: the generated signal
% v: eigenvectors of cov
% d: eigenvalues of cov
% crt: the matrix [n n] by which uncorrelated random numbers are transformed
%
%   See also:  GNORMCOR.
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
%
gr=normrnd(0,1,n,npts);
gn=gr./repmat(sqrt(sum(gr.^2,1)),n,1);
x=sqrt(n)*crt*gn; %sqrt(n) to take into account that this is an elliptical shell
return


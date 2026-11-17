function y=persp_apply_pwr(a,b,c,x,pwrs)
% y=persp_apply_pwrs(a,b,c,x,pwrs) applies a generalized perspective (projective) transformation to one or more vectors
%  The generalization is that a power law may be applied to the denominator
%
% Args:
%   a (float 2-D array): array of size [dimx, dimy], affine component of the transformation
%
%   b (float 1-D array): array of size [   1, dimy], offset component of the transformation
%
%   c (float 1-D array): array of size [dimx,    1], perspective component of the transformation
%
%   x (float 2-D array): array of size [npts, dimx], the (row) vectors to be transformed; npts is the number of vectors to be transformed
%
%   pwrs: if omitted, denominator is 1+xc, if included, (1+(xc)^pwr(1)).^pwr(2); [1 1] if omitted
%
% Returns:
%   y (float 2-D array): array of size [npts, dimy], its row vectors are the results of transformation
%
% Note: How the transformation is defined
%   Each row of x is considered as a homogeneous vector with an augmented coordinate 1 at the end, creating an array of size [npts, dimx+1].
%   This is a standard projective transformation if pwrs is omitted or equal to [1 1]
%
%   A matrix T, of size [dimx+1, dimy+1] is formed with a in its upper left, b in its lower left, c in its upper right, and 1 in its lower right.
%
%   X is then post-multiplied by T, Y=XT.
%
%   Y is a homogeneous vector.  y is the rows of Y, with each row divided by the final element.
%
%   T only matters up to homogeneity, but its lower right element is fixed at 1.
%
%   If c is zero, this is an affine transformation.
%
%   If c and b are zero, this is a linear transformation.
%
% See also:  PERSP_APPLY, PERSP_XFORM_FIND, PERSP_SSQDIF, PERSP_FIT, PERSP_SSQDIF_FIT.
%
if (nargin<=4)
    pwrs=[1 1];
end
T=[a c;b 1];
X=[x,ones(size(x,1),1)];
Y=X*T;
if any(pwrs~=1)
    yd=Y(:,end)-1;
    yd=(1+yd.^pwrs(1)).^pwrs(2);
    y=Y(:,1:end-1)./repmat(yd,1,size(a,2));
else
    y=Y(:,1:end-1)./repmat(Y(:,end),1,size(a,2));
end
return
end
    
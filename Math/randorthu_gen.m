function [rotm,rsign,basis]=randothu_gen(d,v,dsign,n)
% function [rotm,rsign,basis]=randothu_gen(d,v,dsign,n) makes d-dimensional
% random othogonal matrices that keep one or more vectors fixed
%
% d: dimension
% v: array of size [d k], k<=d, of column vectors to be kept fixed
%   columns of v are assumed to be linearly independent
% dsign: dsign=1 -> force determinant to be 1 (a rotation)
%        dsign=-1 -> force determinant to be -1
%        dsign=0 or omitted -> determinant not forced
% n: number of transformations to generate (1 if omitted)
%
% rotm(d,d,n): each rotm(:,:,i) is orthogonal, and rotm(:,:,k)*v=v
% rsign(i): determinant of rotm(:,:,i)
% basis: auxiliary matrix used in creating rotm
%
% Notes:
%   * Compatibility:
%      d, dsign, rotm(:,:,1), rsign(1), as in randorthu
%      v and basis as in extorthb_gen
%   * Rotations are uniformly distributed in SO(d-k).
%   * Cautions:
%       If k=d, the only transformation that preserves v is the identity.
%         So dsign=-1 will generate an error message, and all of rotm will be the identity.
%       If k=d-1, the only transformatoins that prserve v are the identity
%         and a mirror-flip in the direction orthogonal to the span of v.
%         So dsign=1 will only generate the identiy, dsign=-1 will only
%         generate this mirror, and dsign=0 will generate (randomly) half
%         of each.  Not dealt with as a special case.
%   * n is provided since producing multiple matrices at the same time saves recomputing
%         basis.
% 
%  See also: RANDORTHU, RANDORTH, EXTORTHB_GEN.
%
if (nargin<=2)
    dsign=0;
end
if (nargin<=3)
    n=1;
end
k=size(v,2);
rotm=zeros(d,d,n);
rsign=zeros(n,1);
basis=eye(d);
if k==0
    for m=1:n
        [rotm(:,:,m),rsign(m)]=randorthu(d,dsign);
    end
elseif k==d
    rotm=repmat(eye(d),[1 1 n]);
    rsign=ones(n,1);
    if dsign==-1
        error('Incompatible arguments, cannot create identity matrices with determinant -1');
    end
else
    basis=extorthb_gen(v);
    for m=1:n
        [q,rsign(m)]=randorthu(d-k,dsign);
        z=eye(d);
        z(k+1:d,k+1:d)=q;
        rotm(:,:,m)=basis*z*basis';
    end
end

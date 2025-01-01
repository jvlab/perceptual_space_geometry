function [rotm,rsign]=randorth(d,dsign)

% [rotm,rsign]=randorth(d,dsign) is a random real d-dimensional orthogonal matrix
% made by orthogonalizing d vectors drawn from a Gaussian distribution
%
%  NOTE THAT THE ROTATIONS ARE NOT UNIFORMLY DISTRIBUTED IN SO(d) SINCE
%  matlab's ORTH does odd things.  FOR UNIFORM DISTRIBUTIONS, USE RANDORTHU.
%  THIS IS NOTED ON 1/18/04.
%
% dsign=1 -> force determinant to be 1 (a rotation)
% dsign=-1 -> force determinant to be -1
% dsign=0 or omitted -> determinant not forced
%
% rsign is the sign of the determinant of rotm
%
%    see also: RANDORTHU.
%
msign=0;
if (nargin >=2)
   msign=dsign;
end
rotm=orth(randn(d,d));
mdet=det(rotm);
rsign=sign(mdet);
if ((msign==0) | (rsign*msign>0))
   return
else
   rotm(:,1)=-rotm(:,1);
   rsign=-rsign;
end

function [rotm,rsign]=randorthu(d,dsign)

% [rotm,rsign]=randorthu(d,dsign) is a random real d-dimensional orthogonal matrix
% made by orthogonalizing d vectors drawn from a Gaussian distribution
%
%  ROTATIONS **ARE** UNIFORMLY DISTRIBUTED IN SO(n); GRMSCMDT IS USED
%  RATHER THAN matlab's ORTH. 
%
% dsign=1 -> force determinant to be 1 (a rotation)
% dsign=-1 -> force determinant to be -1
% dsign=0 or omitted -> determinant not forced
%
% rsign is the sign of the determinant of rotm
%
%    see also: RANDORTH, GRMSCMDT.
%
msign=0;
if (nargin >=2)
   msign=dsign;
end
w=grmscmdt(randn(d,d));
%
rotm=w./sqrt(repmat(diag(w'*w)',d,1));
%
mdet=det(rotm);
rsign=sign(mdet);
if ((msign==0) | (rsign*msign>0))
   return
else
   rotm(:,1)=-rotm(:,1);
   rsign=-rsign;
end

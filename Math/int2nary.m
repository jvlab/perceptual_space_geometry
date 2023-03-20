function nary=int2nary(ints,d,pmin)
% nary=int2nary(ints,d,pmin) converts one or more integers 
% to n-ary vectors, with the vectors running along the first unused dimesion of ints
% as integers base d, to integers.  Works also for generalized basis.
%   The least significant bit is in position 1.
%
% ints: the integers
% d:      the base, 2 if omitted or blank
%    if d is a vector, then d(1) is used for the first place, d(2) for the second, etc.
%    and d(end) used if d does not have sufficient length
% pmin:   minimum length of the created dimension.  If omitted, minimum possible value is used
% nary:  the converted array of n-ary numbers
%    ints(k)=nary(k,1)+d(1)*nary(k,2)+d(1)*d(2)*nary(k,3)+...
%
% See also:  NARY2INT.
%
if (nargin<=1) d=2; end
if isempty(d) d=2; end
if (nargin<=2) pmin=0; end
nd=ndims(ints);
if (nd==2) & (size(ints,2)==1)
    nd=1;
end
si=size(ints);
si=si(1:nd);
nary=mod(ints,d(1));
ints=(ints-nary)/d(1);
while (max(ints(:))>0)
    if length(d)>1
        d=d(2:end);
    end
    nextdig=mod(ints,d(1));
    nary=cat(nd+1,nary,nextdig);
    ints=(ints-nextdig)/d(1);
end
% pad if necessary
if size(nary,nd+1)<pmin
    nary=cat(nd+1,nary,zeros([si(1:nd) pmin-size(nary,nd+1)]));
end
return
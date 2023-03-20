function ints=nary2int(narys,d,lastdim)
% ints=nary2int(narys,d,lastdim) converts one or more vectors, considered
% as integers base d, to integers.  Also for "generalized" bases (d a vector)
%   One value of ints is calculated for each element in the first
%   ndims(narys)-1 dimensions of narys.  The last dimension of narys is
%   the base-d expansion, with the least significant bit in position 1.
%
% narys:  an array of n-ary numbers
% d:      the base, 2 if omitted
%    if d is a vector, then d(1) is used for the first place, d(2) for the second, etc.
%    and d(end) used if d does not have sufficient length
% lastdim:  force a value for the last dimension (required if, for example,
%    size(narys)=[10 10] but this is thought of as [10 10 1]
%
% ints: the converted integers.  size(ints)=size(narys)[1:ndims(narys)-1]
%    ints(k)=nary(k,1)+d(1)*nary(k,2)+d(1)*d(2)*nary(k,3)+...
%
% See also:  INT2NARY.
%
if (nargin<=1) d=2; end
sn=size(narys);
if (nargin<=2) lastdim=length(sn); end
if (length(sn)<lastdim)
    sn([(length(sn)+1):lastdim])=1;
end
if (length(sn)>lastdim)
    error(sprintf(' attempt to convert n-ary array of dimension %4.0f to integer over non-last dimension %4.0f.',...
    length(sn),lastdim));
end
dprod=[1 d];
if (length(dprod)>sn(end))
    dprod=dprod(1:sn(end));
end
if (length(dprod)<sn(end))
    dprod=[dprod dprod(end)*ones(1,(sn(end)-length(dprod)))];
end
pwrs=cumprod(dprod);
snres=ones(1,lastdim);
snres(end)=sn(end); %[1 1 1 1... 1 sn(end)]
ints=sum(narys.*repmat(reshape(pwrs,snres),[sn([1:(end-1)]) 1]),lastdim);
%
return


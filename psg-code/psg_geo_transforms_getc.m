function c=psg_geo_transforms_getc(dim_max,T,vcuts,acuts)
% c=psg_geo_transforms_getc(dim_max,T,vcuts,acuts) is a utility to adjust offsets so that 
% piecewise linear transforms match on the cut hyperplane
%
% Note that this assumes that the transformations differ from each other by
% a transformation of the form vcut'*tdif.
% May fail if vcut is along the axes
%
% Method: find a vector vtry whose projection onto each row of vcut is equal to acut
% and ensure same value with each transformation
% vcuts is [ncuts dim_max]
% T is [dim_max dim_max nT], where nT is number of pieces (typically 2^ncuts)
% c is of size [nT dim_max]
%
%  19Jan26: attempt to use EXTORTHB if additional constraints are needed
%
%  See also: GEO_TRANSFORMS_TRANSFORMS_SETUP, EXTORTHB.
%
ncuts=size(vcuts,1);
nextra=dim_max-ncuts;
if nextra==0 %just enough constraints
    vtry=acuts*inv(vcuts');
else
    % NOT correct
    % vcuts_aug=extorthb(vcuts');
    % acuts_aug=[acuts vcuts_aug(1:ncuts,end-nextra+1:end)];
    % vtry=acuts_aug*inv(vcuts_aug);
% correct but may fail if vcuts is parallel to axes
    vaug=eye(dim_max); %need to add constraints
    vtry=[acuts zeros(1,nextra)]*inv([vcuts' vaug(end-dim_max+1:end,end-nextra+1:end)]); 
end
c=zeros(size(T,3),dim_max);
for iT=1:size(T,3)
    c(iT,:)=-vtry*T(:,:,iT);
end
return

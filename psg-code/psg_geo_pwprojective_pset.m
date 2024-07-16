function [p_list,sign_table]=psg_geo_pwprojective_pset(vcut,p0,w)
% [p_list,sign_table]=psg_geo_pwprojective_pset(vcut,p0,w) sets up
% projection parameters that are consistent for one or more coutplanes in a
% piecewise projective model
% 
% vcut: [ncuts dim_x], stack of row vectors of length dim_x, orthogonal to cut planes
%   The lengths of vcut must all be 1.
% p0: array of size size [dim_x 1], starting point for projection params
% w: [1 ncuts], the amount of discontinuity on each cutplane
%   w*vcut is added to p0 if x*vcut'>a and subtracted if x*vcut'<a
%    [see psg_piecwise_notes.doc]
%
% p_list: [dim_x 2^ncuts], the projection parameter for each region
% sign_table: [2^ncuts ncuts]: a table of +/-1's, in which the kth row
%   indicates the contribution of the kth cutplane. 
%   [+ + .... +]
%   [- + .... +]
%   [+ - .... +]
%   [- - .... +]
% ...
%   [- - .... -]
%
%    See also:  PSG_GEO_PWPROJECTIVE_VA, PSG_PWPROJECTIVE_APPLY, PSG_PWPROJECTIVE_TEST.
%
ncuts=size(vcut,1);
n_pw=2^ncuts; %number of regions
sign_table=ones(n_pw,ncuts);
for icut=1:ncuts
    sign_table(:,icut)=repmat([ones(2^(icut-1),1);-ones(2^(icut-1),1)],2^(ncuts-icut),1);
end
sign_weighted=sign_table.*repmat(w,n_pw,1); 
p_aug=sign_weighted*vcut;
p_list=repmat(p0,1,n_pw)+p_aug';
return
end

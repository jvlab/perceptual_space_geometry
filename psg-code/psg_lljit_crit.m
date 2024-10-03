function [jit_crit,lljit,opts_lljit_used]=psg_lljit_crit(pval,coords,typenames,responses,stim_list,opts_lljit)
% [pval,lljit,opts_lljit_used]=psg_lljit_pval(pval,coords,typenames,responses,stim_list,opts_lljit)
% computes critical jitter (per dimension, as RMS) to reduce the log likelihood of a model by log(pval)
%
%  Note that log likelihoods are computed log 2, for consistency with SAW software
%
% For multiple p-values, or o see internals, use psj_lljit.
%
% pval: desired p-value to find critical jitter for
% coords: coordinate data, coordsistim,idim) is the coordinate of stimulus istim on dimension idim
% typenames: cell{nstims,1}, stimulus labels from metadata associated with ds, typically from psg_get_coordsets
% responses: choice data, from choice file, 5 or 6 columns
% stim_list: stimulus labels in choice data, character array, size(stim_list,1)=nstims
% opts_lljit: options
%    See psg_lljit.  Same default values, except for jit_list and pvals
%      opts_lljit.pvals set to pval
%      opts_lljit.jit_list set to 1/4 of minimum nonzero distance between points
%
% jit_crit: critical jitter
% lljit: results structure from psg_lljit
%
% opts_lljit_used: options used
%
%  See also: PSG_LLJIT_DEMO, PSG_LLJIT, COOTODSQ.
%
if (nargin<=6)
    opts_lljit=struct;
end
%
dsq=cootodsq(coords); %compute all pairwise distances
dists=sqrt(sort(unique(dsq)));
dnz=dists(2);
%
opts_lljit.jit_list=dnz/4; %1/4 of smallest nonzero distance
opts_lljit.pvals=pval;
%
ds{1}=coords;
[lljit,opts_lljit_used]=psg_lljit(ds,typenames,responses,stim_list,opts_lljit);
jit_crit=lljit.jit_crits;
return
end

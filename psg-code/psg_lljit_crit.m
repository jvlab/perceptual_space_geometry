function [jit_crit,lljit,opts_lljit_used]=psg_lljit_crit(pval,coords,typenames,responses,stim_list,opts_lljit)
% [pval,lljit,opts_lljit_used]=psg_lljit_pval(pval,coords,typenames,responses,stim_list,opts_lljit)
% computes confidence limits for coordinates of a perceptual space, two ways
% expressed as critical jitter (per dimension, as RMS)
% 
% critical jitter 1: find the amount of jitter that changes the log likelhiood of the model by log(pval)
% critical jitter 2: find the bootstrap distribution of the log likelihood, assume it's a Gaussian,
%   see how much of a change in log likelihood is at the two-tailed p-value of that distribution,
%   and then how much coordinate jitter changes the original model's log likelihood by that amount, based on
%   the dependence of the log likelihood on jitter at the minimum
%
%  Note that log likelihoods are computed log 2, for consistency with SAW software
%
% For multiple p-values, or to see internals, use psg_lljit.
%
% psg_lljit is called twice:  once with no jitters, to compute bootstrap
% mean and variance for log likelihood, and once to regress jitters against
% sqrt(delta-log-likelihood) to determine critical jitters
%
% pval: desired p-value to find critical jitter for
% coords: coordinate data, coordsistim,idim) is the coordinate of stimulus istim on dimension idim
% typenames: cell{nstims,1}, stimulus labels from metadata associated with ds, typically from psg_get_coordsets
% responses: choice data, from choice file, 5 or 6 columns
% stim_list: stimulus labels in choice data, character array, size(stim_list,1)=nstims
% opts_lljit: options
%    See psg_lljit.  Same default values, except for jit_list and pvals
%      opts_lljit.pvals set to pval
%      opts_lljit.jit_list set to minimum nonzero distance between
%      points*(1/16,1/8,1/4)
%
% jit_crit: critical jitters
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
jit_crit=zeros(1,2);
%
%preliminary call to determine bootstrap values
%
ds{1}=coords;
[lljit_prelim,opts_lljit_prelim]=psg_lljit(ds,typenames,responses,stim_list,setfield(opts_lljit,'jit_list',0));
%calculations to determine equivalent p-value for type 2
stdv_mult=norminv(1-pval/2); %two-tailed critical value multiplier for pval
lk2_drop=stdv_mult*lljit_prelim.lk2_bootstdv;% two-tailed
pval_equiv=2^(-lk2_drop); %equivalent p-value for reduction of log likelihood by bootstrap estimate
%
%definitive call to find jitters corresponding to pval and pval_equiv
dsq=cootodsq(coords); %compute all pairwise distances
dists=sqrt(sort(unique(dsq)));
dnz=dists(2);
opts_lljit.pvals=[pval pval_equiv];
opts_lljit.jit_list=dnz*[1/16 1/8 1/4]; %fractions of smallest nonzero distance
[lljit,opts_lljit_used]=psg_lljit(ds,typenames,responses,stim_list,opts_lljit);
%
%critical jitters type 1 and 2
jit_crit=lljit.jit_crits(:)';
%
return
end

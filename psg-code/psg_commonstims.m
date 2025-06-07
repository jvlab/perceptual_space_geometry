function [ds_common,sas_common,opts_used]=psg_commonstims(ds,sas,opts)
% [ds_common,sas_common,opts_used]=psg_commonstims[ds,sas,opts_used] selects and aligns the common stimuli from two datasets
%
% ds: cell array (one cell for each dataset) of coordinate data (typically from psg_get_coordsets)
%  All datasets must have dimension lists beginning at 1 and without gaps
% sas: corresponding cell array of metadata (typically from psg_get_coordsets)
% opts: options
%   if_log: 1 (default) to log
%
% ds_common: datasets with only the stimuli in common between all of the ds
% sas_common: corresponding metadata
% opts_used: options used
%
%   See also:  PSG_ALIGN_COORDSETS, PSG_GEOMODELS_RUN, PSG_GEOMODELS_DEMO, PSG_MAJAXES.
%
if (nargin<=2)
    opts=struct;
end
opts=filldefault(opts,'if_log',1);
opts.min='all';
if opts.if_log
    disp('finding stimuli in common')
end
[sets_align,ds_common,sas_common,ovlp_array,sa_pooled,opts_used]=psg_align_coordsets([],ds,sas,opts);
return

function [ds_common,sas_common,opts_used]=psg_commonstims(ds,sas,opts)
% [ds_common,sas_common,opts_used]=psg_commonstims[ds,sas,opts_used] selects and aligns the common stimuli from two datasets
%
% If no alignment is needed then ds_common and sas_common are the original inputs.
% This overrides the behavior of psg_align_coordsets, which alphabetizes in
%  the order of typenames.
%
% ds: cell array (one cell for each dataset) of coordinate data (typically from psg_get_coordsets)
%  All datasets must have dimension lists beginning at 1 and without gaps
% sas: corresponding cell array of metadata (typically from psg_get_coordsets)
% opts: options
%   if_log: 1 (default) to log
%   force_reorder: 0 (default), 1 will force a reordering even if all stimuli match
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
opts=filldefault(opts,'force_reorder',0);
opts.min='all';
if opts.if_log
    disp('finding stimuli in common')
end
[sets_align,ds_common,sas_common,ovlp_array,sa_pooled,opts_used]=psg_align_coordsets([],ds,sas,opts);
nsets=length(sas);
reorder=opts.force_reorder;
for iset=1:nsets
    if length(sas{iset}.typenames)~=length(sas_common{iset}.typenames)
        reorder=1;
    end
end
if opts.if_log
    if reorder==0
        disp('all typenames are in all sets, no reordering done.');
        sas_common=sas;
        ds_common=ds;
    else
        disp('typenames common to all sets are found and now in alphabetical order.');
    end
end
return

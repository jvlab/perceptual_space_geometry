%psg_align_knit_demo: demonstration of alignment and knitting together of multiple datasets
% that have partially overlapping stimuli
%
% all datasets must have dimension lists beginning at 1 and without gaps
% aligned datasets and metadata (ds_align,sas_align) will have a NaN where there is no match
%
%  See also: PSG_ALIGN_COORDSETS, PSG_COORD_PIPE_PROC, PSG_GET_COORDSETS, PSG_READ_COORDDATA, PROCRUSTES_CONSENSUS_PTL_TEST.
%
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_align') opts_align=struct(); end %for psg_align_coordsets
if ~exist('opts_nonan') opts_nonan=struct; end %for psg_remnan_coordsets
%
disp('This will attempt to knit together two or more coordinate datasets.');
%
opts_read=filldefault(opts_read,'input_type',0); %either experimental data or model
opts_align=filldefault(opts_align,'if_log',1);
opts_nonan=filldefault(opts_nonan,'if_log',1);
%
nsets=getinp('number of datasets','d',[1 100]);
[sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,[],nsets); %get the datasets
[sets_align,ds_align,sas_align,ovlp_array,opts_align_used]=psg_align_coordsets(sets,ds,sas,opts_align); %align the stimulus names
%do a consensus calculation here
[sets_nonan,ds_nonan,sas_nonan,opts_nonan_used]=psg_remnan_coordsets(sets_align,ds_align,sas_align,ovlp_array,opts_nonan); %remove the NaNs


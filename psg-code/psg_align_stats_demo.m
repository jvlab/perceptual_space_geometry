%psg_align_stats_demo: demonstration of alignment and knitting together of multiple datasets
% that have partially overlapping stimuli
%
% Does a consensus alignment of overlapping data that need to be knitted together, i.e., not
% all stimuli are present in each condition, and writes the consensus
% data and metadata file.  Assumes that this is a raw data or model file, no previous entries in pipeline.
%
% Compared to psg_align_knit_demo, this is designed for datasets that have largely the same stimuli, though perhaps missing a few --
%  rather than constructing a space that is of higher dimension than any of the component datasets.
%  The resulting consensus dataset is considered as a "denoised" version of the components, rather than an augmented one.
%  Therefore, compared to psg_align_knit_demo:
%   * Does analysis with and without allowing scale
%   * Computes variance explained, by dataset and stimulus
%   * Uses a shuffle of stimuli within datasets to determine whether variance explained by each dimension is significant
%   * Does NOT allow for creation of a consensus dataset that has higher
%     dimension than any component (but this could be changed in the future)
%   * Does not do visualizations
%   * Does not use ray descriptors
%   * Component datasets not stripped of NaN's
%
% Notes:
%  All datasets must have dimension lists beginning at 1 and without gaps
%  Aligned datasets and metadata (ds_align,sas_align) will have a NaN where there is no match
%  For classes such as mater and domain and aux, there should be the same
%      number of stimuli in each file, since btc_specoords is an identity matrix, with one entry for each stimulus.
%
% 20Nov24: add option for frozen random numbers
% 25Nov24: modularize writing of consensus; clean up variable names
% 29Nov24: added if_normscale (disabled by default)
% 26Oct25: fix a bug affecting component coordinates when scaling allowed
% 26Oct25: begin modularizing with psg_align_stats. modularization to allow for
%   choice of dimension in original data to knit
%   if_normscale (if scaling is allowed, this normalizes final dataset to size of input datasets)
%   choice of dimension to create
% 27Oct25: begin modularizing plot
%
%  See also: PSG_ALIGN_COORDSETS, PSG_COORD_PIPE_PROC, PSG_GET_COORDSETS, PSG_READ_COORDDATA,
%    PROCRUSTES_CONSENSUS, PSG_WRITE_COORDDATA, PSG_COORD_PIPE_UTIL, PSG_ALIGN_KNIT_DEMO, PSG_ALIGN_STATS,
%    PSG_ALIGN_PLOT.
%

%main structures and workflow:
%ds{nsets},                sas{nsets}: original datasets and metadata
%ds_align{nsets},          sas_align{nsets}: datasets with NaN's inserted to align the stimuli
%ds_knitted{ia},            sa_pooled: consensus rotation of ds_align, all stimuli, and metadata; ia=1 for no scaling, ia=2 for scaling
%ds_components{ia}{nsets}, sas_align{nsets}: components of ds_knitted, corresponding to original datasets, but with NaNs -- these are Procrustes transforms of ds_align
%
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_align') opts_align=struct(); end %for psg_align_coordsets
if ~exist('opts_pcon') opts_pcon=struct(); end % for procrustes_consensus
%
if ~exist('shuff_quantiles') shuff_quantiles=[0.01 0.05 0.5 0.95 0.99]; end %quantiles for showing shuffled data
nquantiles=length(shuff_quantiles);
%
if ~exist('label_shorten') label_shorten={'coords','hlid_','odor17_','megamat0','6pt','3pt','sess01_10','sess01_20','__','_'}; end %strings to remove from labels
if ~exist('label_replace') label_replace={''      ,''     ,''       ,''        ,''   ,''   ,''         ,''         ,'_' ,'-'}; end %strings to replace
%
disp('This will attempt to knit together two or more coordinate datasets and do statistics.');
if_normscale=getinp('1 to normalize consensus with scaling to size of data','d',[0 1],0);
opts_pcon.if_normscale=if_normscale;
%
if ~exist('nshuffs') nshuffs=500; end
%
nshuffs=getinp('number of shuffles','d',[0 10000],nshuffs);
%
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
%
opts_read.input_type=1;
opts_align=filldefault(opts_align,'if_log',1);
%
opts_pcon=filldefault(opts_pcon,'allow_reflection',1);
opts_pcon=filldefault(opts_pcon,'allow_offset',1);
opts_pcon=filldefault(opts_pcon,'allow_scale',0); 
opts_pcon=filldefault(opts_pcon,'max_niters',1000); %nonstandard max
%
[sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,[],[],0); %get the datasets
nsets=length(sets); %number of sets actually read
%
% check that all dimensions are present
%
for iset=1:nsets
    if (iset==1)
        dim_list_all=sets{iset}.dim_list;
    else
        dim_list_all=intersect(dim_list_all,sets{iset}.dim_list);
    end
    if length(dim_list_all)~=length(1:max(dim_list_all))
        disp(sprintf('some dimensions are missing in set %1.0f',iset))
    end
end
max_dim_all=max(dim_list_all); %max dimension available across all sets
%
% tally missing stimuli in input datasets and align according to stimuli present
%
nstims_each=zeros(1,nsets);
stims_nan=cell(1,nsets);
disp('before alignment of stimuli')
for iset=1:nsets
    nstims_each(iset)=sas{iset}.nstims;
    stims_nan{iset}=find(isnan(ds{iset}{1}));
    disp(sprintf('set %2.0f: %2.0f stimuli (%2.0f are NaN), label: %s',iset,nstims_each(iset),length(stims_nan{iset}),sets{iset}.label))
end
[sets_align,ds_align,sas_align,ovlp_array,sa_pooled,opts_align_used]=psg_align_coordsets(sets,ds,sas,opts_align); %align the stimulus names
nstims_all=sets_align{1}.nstims;
disp(sprintf('total stimuli: %3.0f',nstims_all));
%
disp('after alignment of stimuli')
stims_nan_align=cell(1,nsets);
stims_each_align=zeros(1,nsets);
for iset=1:nsets
    nstims_each_align(iset)=sas_align{iset}.nstims;
    stims_nan_align{iset}=find(isnan(ds_align{iset}{1}));
    disp(sprintf('set %2.0f: %2.0f stimuli (%2.0f are NaN), label: %s',iset,nstims_each_align(iset),length(stims_nan_align{iset}),sets_align{iset}.label))
end
disp('overlap matrix')
disp(ovlp_array'*ovlp_array);
%
pcon_dim_max=getinp('maximum dimension for the consensus alignment dataset to be created (same dimension used in each component)','d',[1 max_dim_all],max_dim_all);
if_ok=0;
while (if_ok==0)
    dim_list_in=getinp('dimensions to use from each dataset for consensus calculation','d',[1 max_dim_all],[1:max_dim_all]);
    dim_list_out=getinp('dimensions to use for the calculated alignment(must be unique and >= dimension in each dataset)','d',[min(dim_list_in) Inf],dim_list_in);
    if length(dim_list_out)==length(dim_list_in)
        if_ok=all(dim_list_out>=dim_list_in) & (length(dim_list_in)==length(unique(dim_list_in))) & (length(dim_list_out)==length(unique(dim_list_out)));
    end
    if if_ok==0
        disp('alignment dimensions must be unique and exceed corresonding input dataset dimensions')
    end
end
dim_list_in_max=max(dim_list_in); %max dimension to analyze
dim_list_out_max=max(dim_list_out); %max dimension to analyze
if_aug_dim=any(dim_list_out>dim_list_in);
%
pcon_init_method=getinp('method to use for initialization (>0: a specific set, 0 for PCA, -1 for PCA with forced centering, -2 for PCA with forced non-centering','d',[-2 nsets],0);
if if_aug_dim
    disp('At least some model dimensions are to be augmented in the consensus.');
else
    disp('Model dimensions in the consensus will be unchanged.');
end
%
if pcon_init_method<=0
    if_initpca_rot=getinp('1 to rotate initialization to match data (best if dimensions are unchanged), 0 to leave unrotated (best if dimensions are augmented)','d',[0 1],1-if_aug_dim);
else
    if_initpca_rot=0;
end
opts_pcon.if_initpca_rot=if_initpca_rot;
if_log=getinp('1 to log details of consensus calculations','d',[0 1]);
%
if pcon_init_method>0
    opts_pcon.initialize_set=pcon_init_method;
else
    if pcon_init_method==0
        opts_pcon.initialize_set='pca';
    elseif pcon_init_method==-1
        opts_pcon.initialize_set='pca_center';
    else
        opts_pcon.initialize_set='pca_nocenter';
    end
end
opts_pcon.nshuffs=nshuffs;
opts_pcon.if_frozen=if_frozen;
opts_pcon.if_log=if_log;
%
%overlaps indicates same stimulus (from ovlp_array) and also
%that the coordinates are not NaN's
z1=NaN(nstims_all,nsets);
for iset=1:nsets
    z1(:,iset)=ds_align{iset}{1}(:,1); 
end
coords_isnan=reshape(isnan(z1),[nstims_all,nsets]);
opts_pcon.overlaps=ovlp_array.*(1-coords_isnan);
if opts_pcon.if_log
    disp(sprintf('number of overlapping stimuli in component removed because coordinates are NaN'));
    disp(sum(coords_isnan.*ovlp_array,1));
    disp('overlap matrix from stimulus matches')
    disp(ovlp_array'*ovlp_array);
    disp(sprintf('overlapping coords in component datasets with values of NaN that are removed from overlaps'));
    disp(opts_pcon.overlaps'*opts_pcon.overlaps);
end
%
%make short labels
%
dataset_labels=cell(1,nsets);
for iset=1:nsets
    dataset_labels{iset}=sets{iset}.label(1+max([find(sets{iset}.label=='/'),find(sets{iset}.label=='\')]):end);
    if exist('label_shorten')
        for id=1:length(label_shorten)
            dataset_labels{iset}=strrep(dataset_labels{iset},label_shorten{id},label_replace{id});
        end
    end
    if dataset_labels{iset}(end)=='-'
        dataset_labels{iset}=dataset_labels{iset}(1:end-1);
    end
end
%
ra_setup=struct; %use as starting point for results and for plotting
%
ra_setup.dataset_labels=dataset_labels;
ra_setup.stimulus_labels=sa_pooled.typenames;
ra_setup.sets=sets;
ra_setup.data_orig=ds;
ra_setup.metadata_orig=sas;
ra_setup.sa_consensus=sa_pooled; %metadata for ds_consensus
ra_setup.sas_components=sas_align;
%
ra_setup.nstims=nstims_all;
ra_setup.nsets=nsets;
ra_setup.nshuffs=nshuffs;
ra_setup.shuff_quantiles=shuff_quantiles; 
%
ra_setup.dim_list_in_max=dim_list_in_max;
ra_setup.dim_list_out_max=dim_list_out_max;
%
results=ra_setup;
%
results.ds_consensus=cell(1,2);
results.ds_components=cell(1,2);
%
%
results.rmsdev_setwise=zeros(dim_list_in_max,nsets,2); %d1: dimension, d2: set, d3: allow_scale
results.rmsdev_stmwise=zeros(dim_list_in_max,nstims_all,2); %d1: dimension, d2: stim, d3: allow_scale
results.rmsdev_overall=zeros(dim_list_in_max,1,2); %rms distance, across all datasets and stimuli
results.counts_setwise=zeros(1,nsets);
results.counts_stmwise=zeros(1,nstims_all);
%descriptors
results.ds_desc='ds_[consensus|components]: top dim is no scaling vs. scaling (opts_pcon)';
results.rmsdev_desc='d1: dimension, d2: nsets or nstims, d3: no scaling vs. scaling (opts_pcon)';
results.counts_desc='d1: 1, d2: nsets or nstims';
%shuffled values
if (nshuffs>0)
    results.rmsdev_shuff_desc='d4: shuffle, d5: 1: shuffle last coords, 2: shuffle all coords';
    results.rmsdev_setwise_shuff=zeros(dim_list_in_max,nsets,2,nshuffs,2); %d1: dimension, d2: set, d3: allow_scale, d4: shuffle, d5: shuffle all coords or last coord
    results.rmsdev_stmwise_shuff=zeros(dim_list_in_max,nstims_all,2,nshuffs,2); %d1: dimension, d2: stim, d3: allow_scale, d4: shuffle, d5: shuffle all coords or last coord
    results.rmsdev_overall_shuff=zeros(dim_list_in_max,1,2,nshuffs,2); %d1: dimension, d2: n/a, d3: allow_scale, d4: shuffle, d5: shuffle last coord or all coords
end
%
%do calculation for each variant of allow_scale, using same permutations
%for each shuffle
%
results.opts_pcon=cell(1,2);
results.dim_list_in=cell(1,2);
results.dim_list_out=cell(1,2);
ra=cell(1,2);
for allow_scale=0:1
    ia=allow_scale+1;
    opts_pcon.allow_scale=allow_scale;
    results.opts_pcon{ia}=opts_pcon;
    results.dim_list_in{ia}=dim_list_in;
    results.dim_list_out{ia}=dim_list_out;
    %
    [ra{ia},warnings]=psg_align_stats(ds_align,sas_align,dim_list_in,dim_list_out,opts_pcon);
    %
    % quantities independent of allow_scale
    results.counts_overall=ra{ia}.counts_overall;
    results.counts_setwise=ra{ia}.counts_setwise;
    results.counts_stmwise=ra{ia}.counts_stmwise;
    results.rmsavail_setwise=ra{ia}.rmsavail_setwise;
    results.rmsavail_stmwise=ra{ia}.rmsavail_stmwise;
    results.rmsavail_overall=ra{ia}.rmsavail_overall;
    % quantities dependent on allow_scale
    results.rmsdev_setwise(:,:,ia)=ra{ia}.rmsdev_setwise;
    results.rmsdev_stmwise(:,:,ia)=ra{ia}.rmsdev_stmwise;
    results.rmsdev_overall(:,:,ia)=ra{ia}.rmsdev_overall;
    results.ds_consensus{ia}=ra{ia}.ds_knitted;
    results.ds_components{ia}=ra{ia}.ds_components;
    if (nshuffs>0) %shuffled values
        results.rmsdev_setwise_shuff(:,:,ia,:,:)=ra{ia}.rmsdev_setwise_shuff(:,:,1,:,:);
        results.rmsdev_stmwise_shuff(:,:,ia,:,:)=ra{ia}.rmsdev_stmwise_shuff(:,:,1,:,:);
        results.rmsdev_overall_shuff(:,:,ia,:,:)=ra{ia}.rmsdev_overall_shuff(:,:,1,:,:);
    end
end
%
ra_setup.nrows=2;
for allow_scale=0:1
    ia=allow_scale+1;
    ra_setup.row=ia;
    if (ia==1)
        figh=psg_align_stats_plot(ra{ia},ra_setup);
    else
        psg_align_stats_plot(ra{ia},setfield(ra_setup,'figh',figh));
    end
end

%
%save consensus as files?
%
ds_knitted=results.ds_consensus; %for compatibilty with psg_consensus_write_util
psg_consensus_write_util;
%
disp('analysis results in ''results'' structure.');
disp('to save with other options, invoke psg_consensus_write_util without leaving the current workspace');

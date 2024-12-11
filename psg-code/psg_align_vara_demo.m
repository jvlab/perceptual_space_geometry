%psg_align_vara_demo: variance analysis after aligning multiple datasets grouped by condition
%
%   need to bring out tags into heatmap labeling and figure legend
%
% Analyzes variance within and between datasets grouped by condition.
% First step is aligning datasets, and allows for knitting, i.e., i.e., not all stimuli are present in each condition, 
% but if a lot of data are missing, this will confound analysis of variance.
% 
% Finds the global consensus, and the consensus within each condition group.
% Computes variance of individual datasets from global consensus, and from its group consensus.
% The former is a statistic for clustering-by-group, and it is then computed across shuffles and compared.
%
%  Compared to psg_align_knit_demo:
%   * Does analysis with and without allowing scale
%   * Does not do visualizations
%   * Component datasets not stripped of NaN's
%   * By default, removes z field from details in saved data
%  Compared to psg_align_stats_demo:
%   * Shuffling is by dataset, not stimulus
%   * By default, removes z field from details in saved data
%   * For shuffles, note that for rms dev by dataset, datasets are shuffled
%   * Option to use exhaustive set of shuffles (if_exhaust in multi_shuff_groups, saved here as results.if_shuff_all),
%     (shuffles that only differ by the order in which they assign datsets to a group are redundant)
%   * Option to restrict shuffles by tagging the datasets, e.g., so that a
%   dataset from one subject is only shuffled with another dataset from that subject
% 
% Notes:
%  All datasets must have dimension lists beginning at 1 and without gaps
%  Aligned datasets and metadata (ds_align,sas_align) will have a NaN where there is no match
%
% 11Dec24: modularize shuffle generation to multi_shuff_groups, and allow for tagging to restrict shuffles.
%
%  See also: PSG_ALIGN_COORDSETS, PSG_GET_COORDSETS, PSG_READ_COORDDATA,
%    PROCRUSTES_CONSENSUS, PSG_ALIGN_KNIT_DEMO, PSG_ALIGN_STATS_DEMO, PSG_ALIGN_VARA_UTIL, MULTI_SHUFF_ENUM, PSG_ALIGN_VARA_PLOT
%    MULTI_SHUFF_GROUPS.
%

%main structures and workflow:
%ds{nsets},                sas{nsets}: original datasets and metadata
%ds_align{nsets},          sas_align{nsets}: datasets with NaN's inserted to align the stimuli
%ds_knitted{ia},            sa_pooled: consensus rotation of ds_align, all stimuli, and metadata; ia=1 for no scaling, ia=2 for scaling
%ds_components{ia}{nsets}, sas_align{nsets}: components of ds_knitted, corresponding to original datasets, but with NaNs -- these are Procrustes transforms of ds_align
%
%ds_align_gp{igp}{ns(igp)}: as in ds_align, but just for the datasets within a group
%ds_knitted_gp{ia}{igp}:    as in ds_knitted{ia}, but just for the datasets within a group
%ds_components_gp{ia}{igp}{ns(igp)}: as in ds_components{ia}{nsets}, but just for the datasets within a group
%
if_debug=getinp('1 for debugging settings','d',[0 1]);
%
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_align') opts_align=struct(); end %for psg_align_coordsets
if ~exist('opts_pcon') opts_pcon=struct(); end % for procrustes_consensus
if ~exist('opts_multi') opts_multi=struct(); end %for multi_shuff_groups
opts_read.input_type=1;
opts_align=filldefault(opts_align,'if_log',1);
opts_pcon=filldefault(opts_pcon,'allow_reflection',1);
opts_pcon=filldefault(opts_pcon,'allow_offset',1);
opts_pcon=filldefault(opts_pcon,'allow_scale',0); 
if (if_debug)
    opts_pcon=filldefault(opts_pcon,'max_iters',100);
else
    opts_pcon=filldefault(opts_pcon,'max_niters',1000); %nonstandard max
end
%
if ~exist('shuff_quantiles') shuff_quantiles=[0.01 0.05 0.5 0.95 0.99]; end %quantiles for showing shuffled data
nquantiles=length(shuff_quantiles);
if ~exist('nshuffs') 
    if if_debug
        nshuffs=10;
    else
        nshuffs=500;
    end
end
%
if ~exist('label_shorten') label_shorten={'coords','hlid_','odor17_','megamat0','6pt','3pt','sess01_10','sess01_20','__','_'}; end %strings to remove from labels
if ~exist('label_replace') label_replace={''      ,''     ,''       ,''        ,''   ,''   ,''         ,''         ,'_' ,'-'}; end %strings to replace
if ~exist('plot_max_factor') plot_max_factor=1; end
%
if ~exist('if_removez') if_removez=1; end %remove the z field from details to shorten saved files
%
disp('This will attempt to knit together two or more coordinate datasets and do statistics.');
if_normscale=getinp('1 to normalize consensus with scaling to size of data','d',[0 1]);
opts_pcon.if_normscale=if_normscale;
%
nshuffs=getinp('default number of shuffles for random shuffling','d',[0 10000],nshuffs);
%
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
if (if_frozen~=0) 
    rng('default');
    if (if_frozen<0)
        rand(1,abs(if_frozen));
    end
else
    rng('shuffle');
end
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
%get grouping information
%
if_ok=0;
while (if_ok==0)
    ngps=getinp('number of groups','d',[2 nsets]);
    gps=zeros(1,nsets);
    gp_list=cell(1,ngps);
    nsets_gp=zeros(1,ngps); %number of datasets in each group
    sets_avail=[1:nsets];
    for igp=1:ngps
        if ~isempty(sets_avail)
            if igp<ngps
                gp_list{igp}=getinp(sprintf('datasets for group %1.0f',igp),'d',[min(sets_avail) max(sets_avail)],sets_avail);
                gp_list{igp}=intersect(gp_list{igp},sets_avail);
                gps(gp_list{igp})=igp;
                sets_avail=setdiff(sets_avail,gp_list{igp});
            else
                gp_list{igp}=sets_avail;
                gps(sets_avail)=igp;
                sets_avail=[];
            end
        end       
    end
    %are all groups used, and is each stim assigned?
    for igp=1:ngps       
        nsets_gp(igp)=length(gp_list{igp});
        disp(sprintf('group %1.0f',igp))
        for iset=gp_list{igp}
            disp(sprintf(' dataset %2.0f: %s',iset,sets{iset}.label));
        end
    end
    if_ok=1;
    if ~isempty(sets_avail)
        disp('not all datasets assigned');
        if_ok=0;
    end
    if max(gps)<ngps
        disp('not all groups used');
        if_ok=0;
    end
    if (if_ok==1)
        if_ok=getinp('1 if ok','d',[0 1]);
    end
end
nsets_gp_max=max(nsets_gp);
%
% get tag info
% 
if_ok=0;
while (if_ok==0)
    if getinp('1 restrict shuffles within tagged subsets','d',[0 1]);
        opts_multi.tags=getinp(sprintf('%2.0f tags',nsets),'d',[1 nsets],ones(1,nsets));
        for iset=1:nsets
            disp(sprintf('dataset %2.0f (group %2.0f, tag %2.0f): %s',iset,gps(iset),opts_multi.tags(iset),sets{iset}.label));
        end
    else
        opts_multi.tags=[];
    end
    if_ok=getinp('1 if ok','d',[0 1]);
end
%
%create shuffles
%
if_ok=0;
while (if_ok==0)
    opts_multi.if_reduce=1;
    opts_multi.if_ask=1;
    opts_multi.if_justcount=0;
    opts_multi.nshuffs=nshuffs;
    [permute_sets,gp_info,opts_multi_used]=multi_shuff_groups(gps,opts_multi);
    nshuffs=opts_multi_used.nshuffs;
    disp(sprintf(' created %5.0f shuffles for %3.0f datasets',nshuffs,nsets));
    if_ok=getinp('1 if ok','d',[0 1]);
end
%
% tally missing stimuli in input datasets and align according to all stimuli
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
pcon_init_method=getinp('method to use for initialization (>0: a specific set, 0 for PCA, -1 for PCA with forced centering, -2 for PCA with forced non-centering','d',[-2 min(nsets_gp)],0);
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
%reformat the data for consenus calculations
z=cell(pcon_dim_max,1);
for ip=1:pcon_dim_max
    z{ip}=zeros(nstims_all,ip,nsets);
    for iset=1:nsets
        z{ip}(:,:,iset)=ds_align{iset}{ip}(:,[1:ip]); %only include data up to pcon_dim_use
        z{ip}(opts_align_used.which_common(:,iset)==0,:,iset)=NaN; % pad with NaN's if no data
    end
end
%
%overlaps indicates same stimulus (from ovlp_array) and also
%that the coordinates are not NaN's
coords_isnan=reshape(isnan(z{1}),[nstims_all,nsets]);
disp(sprintf('number of overlapping stimuli in component removed because coordinates are NaN'));
disp(sum(coords_isnan.*ovlp_array,1));
opts_pcon.overlaps=ovlp_array.*(1-coords_isnan);
%
disp('overlap matrix from stimulus matches')
disp(ovlp_array'*ovlp_array);
disp(sprintf('overlapping coords in component datasets with values of NaN that are removed from overlaps'));
disp(opts_pcon.overlaps'*opts_pcon.overlaps);
%
results=struct;
consensus=cell(pcon_dim_max,2); %d1: dimnension, d2: allow scale=[0,1]
znew=cell(pcon_dim_max,2);
ts=cell(pcon_dim_max,2);
details=cell(pcon_dim_max,2);
opts_pcon_used=cell(pcon_dim_max,2);
ds_knitted=cell(1,2);
ds_components=cell(1,2); %allow scale or not
%
%these are vector distances, taking all coordinates into account
%
rmsdev_setwise=zeros(pcon_dim_max,nsets,2); %d1: dimension, d2: set, d3: allow_scale
rmsdev_stmwise=zeros(pcon_dim_max,nstims_all,2); %d1: dimension, d2: stim, d3: allow_scale
rmsdev_overall=zeros(pcon_dim_max,1,2); %rms distance, across all datasets and stimuli
%
rmsdev_setwise_gp=zeros(pcon_dim_max,nsets_gp_max,2,ngps); %d1: dimension, d2: set (within group), d3: allow_scale, d4: gp
rmsdev_stmwise_gp=zeros(pcon_dim_max,nstims_all,2,ngps); %d1: dimension, d2: stim, d3: allow_scale, d4: gp
rmsdev_overall_gp=zeros(pcon_dim_max,1,2,ngps); %d1: dimension, d3: allow_scale, d4: gp
%
counts_setwise=zeros(1,nsets);
counts_stmwise=zeros(1,nstims_all);
counts_overall=zeros(1);
%
counts_setwise_gp=zeros(1,nsets_gp_max,ngps);
counts_stmwise_gp=zeros(1,nstims_all,ngps);
counts_overall_gp=zeros(1,1,ngps);
%
%for setwise and stmwise, use NaN so that missing data won't affect averages
rmsdev_setwise_gp_shuff=NaN(pcon_dim_max,nsets_gp_max,2,ngps,nshuffs); %d1: dimension, d2: set, d3: allow_scale, d4: gp, d5: shuffle
rmsdev_stmwise_gp_shuff=NaN(pcon_dim_max,nstims_all,2,ngps,nshuffs); %d1: dimension, d2: stim, d3: allow_scale, d4: gp, d5: shuffle
rmsdev_overall_gp_shuff=zeros(pcon_dim_max,1,2,ngps,nshuffs); %d1: dimension, d3: allow_scale, d4: gp, d5: shuffle
%
%rms variance available in original data
rmsavail_setwise=zeros(pcon_dim_max,nsets);
rmsavail_stmwise=zeros(pcon_dim_max,nstims_all);
rmsavail_overall=zeros(pcon_dim_max,1);
for ip=1:pcon_dim_max
    sqs=sum(z{ip}.^2,2);
    rmsavail_setwise(ip,:)=reshape(sqrt(mean(sqs,1,'omitnan')),[1 nsets]);
    rmsavail_stmwise(ip,:)=reshape(sqrt(mean(sqs,3,'omitnan')),[1 nstims_all]);
    rmsavail_overall(ip,:)=sqrt(mean(sqs(:),'omitnan'));
end
%
for allow_scale=0:1
    ia=allow_scale+1;
    disp(' ')
    disp(sprintf(' calculations with allow_scale=%1.0f',allow_scale));
    if (allow_scale==1)
        disp(sprintf(' if_normscale=%1.0f',if_normscale));
    end
    opts_pcon.allow_scale=allow_scale;
    ds_knitted{ia}=cell(1,pcon_dim_max);
    ds_components{ia}=cell(1,nsets);
    for ip=1:pcon_dim_max
        %find the global consensus, this is independent of shuffle
        [consensus{ip,ia},znew{ip,ia},ts{ip,ia},details{ip,ia},opts_pcon_used{ip,ia}]=procrustes_consensus(z{ip},opts_pcon);
        if if_removez
            details{ip,ia}=rmfield(details{ip,ia},'z');
        end
        disp(sprintf(' creating global Procrustes consensus for dim %2.0f based on component datasets, iterations: %4.0f, final total rms dev per coordinate: %8.5f',...
            ip,length(details{ip,ia}.rms_change),sqrt(sum(details{ip,ia}.rms_dev(:,end).^2))));
        ds_knitted{ia}{ip}=consensus{ip,ia};
        for iset=1:nsets
            ds_components{ia}{iset}{1,ip}=znew{ip}(:,:,iset);
        end
        sqdevs=sum((znew{ip,ia}-repmat(consensus{ip,ia},[1 1 nsets])).^2,2); %squared deviation of consensus from rotated component
        %rms deviation across each dataset, summed over coords, normalized by the number of stimuli in each dataset
        rmsdev_setwise(ip,:,ia)=reshape(sqrt(mean(sqdevs,1,'omitnan')),[1 nsets]);
        counts_setwise=squeeze(sum(~isnan(sqdevs),1))';
        %rms deviation across each stimulus, summed over coords, normalized by the number of sets that include the stimulus
        rmsdev_stmwise(ip,:,ia)=reshape(sqrt(mean(sqdevs,3,'omitnan')),[1 nstims_all]);
        counts_stmwise=(sum(~isnan(sqdevs),3))';
        %rms deviation across all stimuli and coords
        rmsdev_overall(ip,1,ia)=sqrt(mean(sqdevs(:),'omitnan'));
        counts_overall=sum(~isnan(sqdevs(:)));
        %
        %do shuffles, shuffle 0 = unshuffled
        %
        for ishuff=0:nshuffs
            if (ishuff==0)
                perm_use=[1:nsets];
                rmsdev_overall_gplims=repmat([Inf,-Inf],[ngps 1]);
            else
                perm_use=permute_sets(ishuff,:);
            end
            zs=z{ip}(:,:,perm_use); %the datasets in permuted order, with NaN's where stimuli are missing
            %zs(istim,idim,iset)
            for igp=1:ngps
                zg=zeros(nstims_all,ip,nsets_gp(igp));
                for iset_ptr=1:nsets_gp(igp)
                    iset=gp_list{igp}(iset_ptr); %a dataset in this group
                    zg(:,:,iset_ptr)=zs(:,:,iset);
                end
                %eliminate any stimuli that are NaN in all zg and create a
                %pointer (stims_gp) to stimuli in this group to the full stimuls set in sa_pooled
                stims_gp=find(~all(any(isnan(zg),2),3)); %if some coord is NaN in all of the datasets
                zg=zg(stims_gp,:,:); %keep only the stimuli that have data
                %now form a consensus from each group
                overlaps_gp=1-reshape(any(isnan(zg),2),[length(stims_gp),nsets_gp(igp)]); %overlaps within group
                [consensus_gp,znew_gp,ts_gp,details_gp]=procrustes_consensus(zg,setfield(opts_pcon,'overlaps',overlaps_gp));
                r=sqrt(sum(details_gp.rms_dev(:,end).^2));
                %
                sqdevs_gp=sum((znew_gp-repmat(consensus_gp,[1 1 nsets_gp(igp)])).^2,2); %squared deviation of group consensus from rotated component
                rms_setwise_gp=reshape(sqrt(mean(sqdevs_gp,1,'omitnan')),[1 nsets_gp(igp)]);
                rms_stmwise_gp=reshape(sqrt(mean(sqdevs_gp,3,'omitnan')),[1 length(stims_gp)]);
                rms_overall_gp=sqrt(mean(sqdevs_gp(:),'omitnan'));
                if (ishuff==0)
                    rmsdev_setwise_gp(ip,[1:nsets_gp(igp)],ia,igp)=rms_setwise_gp;
                    counts_setwise_gp(1,[1:nsets_gp(igp)],igp)=squeeze(sum(~isnan(sqdevs_gp),1))';
                    %rms deviation across each stimulus, summed over coords, normalized by the number of sets that include the stimulus
                    rmsdev_stmwise_gp(ip,stims_gp,ia,igp)=rms_stmwise_gp;
                    counts_stmwise_gp(1,stims_gp,igp)=(sum(~isnan(sqdevs_gp),3))';
                    %rms deviation across all stimuli and coords
                    rmsdev_overall_gp(ip,1,ia,igp)=rms_overall_gp;
                    counts_overall_gp(1,1,igp)=sum(~isnan(sqdevs_gp(:)));
                    %
                    disp(sprintf('  grp %2.0f: %3.0f datasets, %3.0f of %3.0f stimuli, Procrustes consensus iterations: %4.0f, final total rms dev per coordinate: %8.5f',...
                        igp,nsets_gp(igp),length(stims_gp),nstims_all,length(details_gp.rms_change),r));
                else
                    rmsdev_setwise_gp(ip,[1:nsets_gp(igp)],ia,igp,ishuff)=rms_setwise_gp;
                    rmsdev_stmwise_gp(ip,stims_gp,ia,igp)=rms_stmwise_gp;
                    rmsdev_overall_gp_shuff(ip,1,ia,igp,ishuff)=rms_overall_gp;
                end
                if (ishuff==nshuffs)
                    disp(sprintf('  grp %2.0f: total rms vec distance in data %8.5f; in %5.0f shuffles, range: [%8.5f %8.5f]',...
                        igp,rmsdev_overall_gp(ip,1,ia,igp),nshuffs,...
                        min(rmsdev_overall_gp_shuff(ip,1,ia,igp,:),[],5),max(rmsdev_overall_gp_shuff(ip,1,ia,igp,:),[],5)));
                end
            end %igp
        end %ishuff
        %
    end %ip 
end %ia
results.nstims=nstims_all;
results.nsets=nsets;
results.nshuffs=nshuffs;
results.dim_max=pcon_dim_max;
results.opts_multi=opts_multi_used;
%grouping informtion
results.ngps=ngps;
results.gps=gps;
results.gp_list=gp_list;
results.nsets_gp=nsets_gp;
%available rms variance in original data
results.rmsavail_setwise=rmsavail_setwise;
results.rmsavail_stmwise=rmsavail_stmwise;
results.rmsavail_overall=rmsavail_overall;
%
results.if_shuff_all=opts_multi_used.if_exhaust;
results.if_normscale=if_normscale;
%
results.ds_desc='ds_[knitted|components]: top dim is no scaling vs. scaling';
results.ds_consensus=ds_knitted;
results.ds_components=ds_components;
results.sa_consensus=sa_pooled; %metadata for ds_knitted
results.sas_components=sas_align; %metadata for each of ds_components
%
results.rmsdev_desc='d1: dimension, d2: nsets or nstims, d3: no scaling vs. scaling';
results.rmsdev_setwise=rmsdev_setwise;
results.rmsdev_stmwise=rmsdev_stmwise;
results.rmsdev_overall=rmsdev_overall;
results.counts_desc='d1: 1, d2: nsets or nstims';
results.counts_setwise=counts_setwise;
results.counts_stmwise=counts_stmwise;
results.counts_overall=counts_overall;
%
results.rmsdev_gp_desc='d1: dimension, d2: nsets or nstims, d3: no scaling vs. scaling, d4: group';
results.rmsdev_setwise_gp=rmsdev_setwise_gp;
results.rmsdev_stmwise_gp=rmsdev_stmwise_gp;
results.rmsdev_overall_gp=rmsdev_overall_gp;
% sum of within-gp rms dev explained weighted by number of sets within each group
results.rmsdev_grpwise=sqrt(sum(rmsdev_overall_gp.^2.*repmat(reshape(nsets_gp(:),[1 1 1 ngps]),[pcon_dim_max,1,2,1]),4)/nsets);
%
results.counts_setwise_gp=counts_setwise_gp;
results.counts_stmwise_gp=counts_stmwise_gp;
results.counts_overall_gp=counts_overall_gp;
%
if (nshuffs>0)
    results.rmsdev_gp_shufff_desc='d1: dimension, d2: nsets or nstims, d3: no scaling vs. scaling, d4: group, d5: shuffle';
    results.rmsdev_setwise_gp_shuff=rmsdev_setwise_gp_shuff;
    results.rmsdev_stmwise_gp_shuff=rmsdev_stmwise_gp_shuff;
    results.rmsdev_overall_gp_shuff=rmsdev_overall_gp_shuff;
    results.rmsdev_grpwise_shuff=sqrt(sum(rmsdev_overall_gp_shuff.^2.*repmat(reshape(nsets_gp(:),[1 1 1 ngps 1]),[pcon_dim_max,1,2,1,nshuffs]),4)/nsets);
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
results.dataset_labels=dataset_labels;
results.stimulus_labels=sa_pooled.typenames;
results.sets=sets;
results.data_orig=ds;
results.metadata_orig=sas;
%
%plotting: only uses quantities in results; shuff_quantiles and if_debug re-declared
% so this section is autonomous
%
psg_align_vara_plot;
%
%save consensus as files?
%
psg_consensus_write_util;
%
disp('analysis results in ''results'' structure.');

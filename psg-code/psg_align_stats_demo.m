%psg_align_stats_demo: demonstration of alignment and knitting together of multiple datasets
% that have partially overlapping stimuli
%
% Does a consensus alignment of overlapping data that need to be 'knitted together, i.e., not
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
%
%  See also: PSG_ALIGN_COORDSETS, PSG_COORD_PIPE_PROC, PSG_GET_COORDSETS, PSG_READ_COORDDATA,
%    PROCRUSTES_CONSENSUS, PSG_WRITE_COORDDATA, PSG_COORD_PIPE_UTIL, PSG_ALIGN_KNIT_DEMO.
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
%
if ~exist('nshuffs') nshuffs=500; end
%
nshuffs=getinp('number of shuffles','d',[0 10000],nshuffs);
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
%make permutations for shuffling
%
permutes=cell(1,nsets);
for iset=1:nsets
    stims_nonan=setdiff(1:nstims_each_align,stims_nan_align{iset});
    permutes{iset}=zeros(max_dim_all,nstims_each_align(iset),nshuffs);
    for ishuff=1:nshuffs
        for ip=1:max_dim_all
            permutes{iset}(ip,stims_nonan,ishuff)=stims_nonan(randperm(length(stims_nonan))); %onl shuffle nonans
        end
    end
    permutes{iset}(:,stims_nan_align{iset},:)=repmat(stims_nan_align{iset}(:)',[max_dim_all 1 nshuffs]); %nan's don't get shuffled
    disp(sprintf(' set %2.0f: created shuffles for %3.0f stimuli',iset,length(stims_nonan)));
end
%
pcon_dim_max=getinp('maximum dimension for the consensus alignment dataset to be created (same dimension used in each component)','d',[1 max_dim_all],max_dim_all);
pcon_init_method=getinp('method to use for initialization (>0: a specific set, 0 for PCA, -1 for PCA with forced centering, -2 for PCA with forced non-centering','d',[-2 nsets],0);
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
%
%reformat data for consensus calculation
%
z=cell(pcon_dim_max,1);
for ip=1:pcon_dim_max
    z{ip}=zeros(nstims_all,ip,nsets);
    for iset=1:nsets
        z{ip}(:,:,iset)=ds_align{iset}{ip}(:,[1:ip]); %only include data up to pcon_dim_use
        z{ip}(opts_align_used.which_common(:,iset)==0,:,iset)=NaN; % pad with NaN's if no data
    end
end
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
counts_setwise=zeros(1,nsets);
counts_stmwise=zeros(1,nstims_all);
%
rmsdev_setwise_shuff=zeros(pcon_dim_max,nsets,2,nshuffs,2); %d1: dimension, d2: set, d3: allow_scale, d4: shuffle, d5: shuffle all coords or last coord
rmsdev_stmwise_shuff=zeros(pcon_dim_max,nstims_all,2,nshuffs,2); %d1: dimension, d2: stim, d3: allow_scale, d4: shuffle, d5: shuffle all coords or last coord
rmsdev_overall_shuff=zeros(pcon_dim_max,1,2,nshuffs,2); %d1: dimension, d2: n/a, d3: allow_scale, d4: shuffle, d5: shuffle last coord or all coords
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
    opts_pcon.allow_scale=allow_scale;
    ds_knitted{ia}=cell(1,pcon_dim_max); %21Nov24: was (1,nsets)
    ds_components{ia}=cell(1,nsets);
    for ip=1:pcon_dim_max
        %do unshuffled
        [consensus{ip,ia},znew{ip,ia},ts{ip,ia},details{ip,ia},opts_pcon_used{ip,ia}]=procrustes_consensus(z{ip},opts_pcon);
        disp(sprintf(' creating Procrustes consensus for dim %2.0f based on component datasets, iterations: %4.0f, final total rms dev per coordinate: %8.5f',...
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
        if nshuffs>0
            %shuffles: across last coord or all coords
            zp=z{ip};
            for ist=1:2 %1: shuffle last coord, 2: shuffle all coords
                if (ist==1)
                    dims_to_shuffle=ip;
                else
                    dims_to_shuffle=[1:ip];
                end
                for ishuff=1:nshuffs
                    zshuff=zp; %start from un-shuffled data
                    for iset=1:nsets
                         perms=permutes{iset}(ip,:,ishuff); %permutes{iset}: d1: dimension, d2: stimulus, d3: which shuffle
                         zshuff(:,dims_to_shuffle,iset)=zp(perms,dims_to_shuffle,iset); %zp: d1 is stimulus, d2 is dimension, d3 is set; permute either last or all dimensions
                    end
                    [consensus_shuff,zn_shuff]=procrustes_consensus(zshuff,opts_pcon);
                    sqdevs=sum((zn_shuff-repmat(consensus_shuff,[1 1 nsets])).^2,2); %squared deviation of consensus from rotated component
                    rmsdev_setwise_shuff(ip,:,ia,ishuff,ist)=reshape(sqrt(mean(sqdevs,1,'omitnan')),[1 nsets]);
                    rmsdev_stmwise_shuff(ip,:,ia,ishuff,ist)=reshape(sqrt(mean(sqdevs,3,'omitnan')),[1 nstims_all]);
                    rmsdev_overall_shuff(ip,1,ia,ishuff,ist)=sqrt(mean(sqdevs(:),'omitnan'));
                end %ishuff
            end %ist
        end %nshuff>0
    end %ip
end %ia
results.nstims=nstims_all;
results.nsets=nsets;
results.nshuffs=nshuffs;
results.dim_max=pcon_dim_max;
%available rms variance in original data
results.rmsavail_setwise=rmsavail_setwise;
results.rmsavail_stmwise=rmsavail_stmwise;
results.rmsavail_overall=rmsavail_overall;
%
results.ds_desc='ds_[knitted|components]: top dim is no scaling vs. scaling';
results.ds_consensus=ds_knitted;
results.ds_components=ds_components;
results.sa_consensus=sa_pooled; %metadata for ds_knitted
results.sas_components=sas_align; %metadata for each of ds_components
results.rmsdev_desc='d1: dimension, d2: nsets or nstims, d3: no scaling vs. scaling';
results.rmsdev_setwise=rmsdev_setwise;
results.rmsdev_stmwise=rmsdev_stmwise;
results.rmsdev_overall=rmsdev_overall;
results.counts_desc='d1: 1, d2: nsets or nstims';
results.counts_setwise=counts_setwise;
results.counts_stmwise=counts_stmwise;
results.counts_overall=counts_overall;
if (nshuffs>0)
    results.rmsdev_shuff_desc='d4: shuffle, d5: 1: shuffle last coords, 2: shuffle all coords'
    results.rmsdev_setwise_shuff=rmsdev_setwise_shuff;
    results.rmsdev_stmwise_shuff=rmsdev_stmwise_shuff;
    results.rmsdev_overall_shuff=rmsdev_overall_shuff;
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
%plotting: should only use quantities in results and shuff_quantiles
%
figure;
set(gcf,'NumberTitle','off');
set(gcf,'Name','consensus analysis');
set(gcf,'Position',[100 100 1300 800]);
ncols=4; %unexplained variance by dataset, unexplained variance by stimulus, unexplained variance and shuffles, explained variance and shuffles
rms_plot_max1=max([max(abs(results.rmsdev_setwise(:))),max(abs(results.rmsdev_stmwise(:)))]); %max possible rms variance per category (dataset or stim)
rms_plot_max2=max(abs(results.rmsdev_overall(:)));
if results.nshuffs>0
    rms_plot_max2=max([rms_plot_max2,max(abs(results.rmsdev_overall_shuff(:)))]);
end
rms_plot_max3=max(results.rmsavail_overall);
for allow_scale=0:1
    ia=allow_scale+1;
    if (allow_scale==0)
        scale_string='no scaling';
    else
        scale_string='with scaling';
    end
    %compare rms devs across datasetsdataset
    subplot(2,ncols,allow_scale*ncols+1);
    imagesc(results.rmsdev_setwise(:,:,ia),[0 rms_plot_max1]);
    xlabel('dataset');
    set(gca,'XTick',1:nsets);
    set(gca,'XTickLabel',results.dataset_labels);
    ylabel('dim');
    set(gca,'YTick',1:results.dim_max);
    title(cat(2,'rms var unex, by set, ',scale_string));
    colorbar;
    %compare rms devs across stimuli
    subplot(2,ncols,allow_scale*ncols+2);
    imagesc(results.rmsdev_stmwise(:,:,ia),[0 rms_plot_max1]);
    xlabel('stim');
    set(gca,'XTick',1:nstims_all);
    set(gca,'XTickLabel',results.stimulus_labels);
    ylabel('dim');
    set(gca,'YTick',1:results.dim_max);
    title(cat(2,'rms var unex, by stim, ',scale_string));
    colorbar;
    for iue=1:2 %unexplained or explained
    %overall rms as function of dimension, and compare with shuffles
        if (iue==1)
            var_string='unexplained'; 
            iue_sign=1; %logic to add unexplained variance
            iue_mult=0; %and not include available variance
            ylim=rms_plot_max2;
        else
            var_string='explained';
            iue_sign=-1; %logic to subtract unexplained variance
            iue_mult=1; %and add explained variance
            ylim=rms_plot_max3;
        end
        subplot(2,ncols,allow_scale*ncols+2+iue);
        hl=cell(0);
        hp=plot(1:results.dim_max,iue_mult*results.rmsavail_overall(:,1)+iue_sign*results.rmsdev_overall(:,1,ia),'k');
        hl=[hl,hp];
        ht='consensus, data';
        hold on;
        if (iue==2)
            hp=plot(1:results.dim_max,results.rmsavail_overall(:,1),'b');
            hl=[hl,hp];
            ht=strvcat(ht,'avail');
        end
        if nshuffs>0
            for iq=1:nquantiles
                switch sign(shuff_quantiles(iq)-0.5)
                    case -1
                        linetype=':';
                    case 0
                        linetype='';
                    case 1
                        linetype='--';
                end
                hp_last=plot(1:results.dim_max,iue_mult*results.rmsavail_overall(:,1)+...
                    iue_sign*quantile(results.rmsdev_overall_shuff(:,1,ia,:,1),shuff_quantiles(iq),4),cat(2,'r',linetype));
                hp_all=plot(1:results.dim_max,iue_mult*results.rmsavail_overall(:,1)+...
                    iue_sign*quantile(results.rmsdev_overall_shuff(:,1,ia,:,2),shuff_quantiles(iq),4),cat(2,'m',linetype));
                if iq==round((1+nquantiles)/2)
                    hl=[hl,hp_last,hp_all];
                    ht=strvcat(ht,'cons, last','cons, all');
                 end
            end
        end
        set(gca,'XTick',1:results.dim_max);
        set(gca,'XLim',[0 results.dim_max]);
        xlabel('dim');
        set(gca,'YLim',[0 ylim]);
        ylabel('rms dev');
        title(cat(2,'rms ',var_string,' overall, ',scale_string));
        legend(hl,ht,'Location','Best','FontSize',7);
    end %iue
end
axes('Position',[0.01,0.04,0.01,0.01]); %for text
text(0,0,'consensus analysis','Interpreter','none','FontSize',8);
axis off;
if (nshuffs>0)
    axes('Position',[0.5,0.04,0.01,0.01]); %for text
    text(0,0,cat(2,sprintf('quantiles from %5.0f shuffles: ',nshuffs),sprintf('%6.4f ',shuff_quantiles)),...
        'FontSize',8);
    axis off;
end
%
%save consensus as files?
%
psg_consensus_write_util;
%
disp('analysis results in ''results'' structure.');


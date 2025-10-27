function [ra,ou,warnings]=psg_align_stats(ds_align,sas_align,dim_list_in,dim_list_out,opts_pcon)
% [ra,ou,warnings]=psg_align_stats(ds_align,sas_align,dim_list_in,dim_list_out,opts_pcon)
% knits together sets of coordinates that may have incompletely overlapping
% stimulus lists, and computes statistics of alignments via shuffles
%
% ds_align: cell array of coordinate sets
%   ds_align{iset}{idim} is an array of nstims x idim, containing the coordinates. 
%   coordinates for missing stimuli should be NaN
% sas_align: cell array of metadata
% dim_list_in: list of dimensions to process.  
% dim_list_out: list of dimennsions to create, one entry for each entry in dim_list_in
%    each must be >=corresponding entry in dim_list_in
% opts_pcon: options.  Most in procrustes_consensus, but also
%   nshuffs: number of shuffles (0: statistics not done), defaults to 500
%   if_frozen: 1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed, defaults to 1
%   if_log: 1 to log, defaults to 1
%
% ra: structure of results
% ou: options used
%   ou.eachdim: cell array of structure of options used, for each dimension in dim_list_in
%   ou.pcon_used: options common to all dimensions
% warnings: warnings
% 
%  See also:  PSG_ALIGN_STATS_DEMO, PSG_ALIGN_COORDSETS, PROCRUSTES_CONSENSUS
%
if nargin<=4
    opts_pcon=struct;
end
opts_pcon=filldefault(opts_pcon,'if_frozen',1);
opts_pcon=filldefault(opts_pcon,'nshuffs',500);
opts_pcon=filldefault(opts_pcon,'if_log',1);
%
if opts_pcon.if_log
    disp(sprintf(' calculations with allow_scale=%1.0f, if_normscale=%1.0f',opts_pcon.allow_scale,opts_pcon.if_normscale));
end
%
ra=struct;
ou.pcon_used=struct;
ou.eachdim=cell(1,length(dim_list_in));
warnings=[];
%
nsets=length(ds_align);
dim_list_in_max=max(dim_list_in);
dim_list_out_max=max(dim_list_out);
%
%set up random number generator
%
if_frozen=opts_pcon.if_frozen;
if (if_frozen~=0) 
    rng('default');
    if (if_frozen<0)
        rand(1,abs(if_frozen));
    end
else
    rng('shuffle');
end
%verify that each set has same number of stimuli (possibly including nans)
nstims_all=sas_align{1}.nstims;
for iset=1:nsets
    if sas_align{iset}.nstims~=nstims_all
        wmsg='number of stimuli (including NaN responses for missing stimuli) differ across datasets after alignment, cannot proceed';
        disp(wmsg);
        warnings=strvcat(warnings,wmsg);
        return
    end
end
%
%determine number of stimuli in each set and make permutations for shuffling
%
stims_nan_align=cell(1,nsets);
stims_each_align=zeros(1,nsets);
permutes=cell(1,nsets);
nshuffs=opts_pcon.nshuffs;
for iset=1:nsets
    nstims_each_align(iset)=sas_align{iset}.nstims;
    stims_nan_align{iset}=find(isnan(ds_align{iset}{1}));
    stims_nonan=setdiff(1:nstims_each_align,stims_nan_align{iset});
    permutes{iset}=zeros(dim_list_in_max,nstims_each_align(iset),nshuffs);
    for ishuff=1:nshuffs
        for dptr=1:length(dim_list_in)
            ip=dim_list_in(dptr);
            permutes{iset}(ip,stims_nonan,ishuff)=stims_nonan(randperm(length(stims_nonan))); %only shuffle nonans
        end
    end
    permutes{iset}(:,stims_nan_align{iset},:)=repmat(stims_nan_align{iset}(:)',[dim_list_in_max 1 nshuffs]); %nan's don't get shuffled
    if opts_pcon.if_log
        disp(sprintf(' set %2.0f: created shuffles for %3.0f stimuli',iset,length(stims_nonan)));
    end
end
%
%reformat data for consensus calculation
%
z=cell(dim_list_in_max);
for dptr=1:length(dim_list_in)
    ip=dim_list_in(dptr);
    z{ip}=NaN(nstims_all,ip,nsets);
    for iset=1:nsets
        z{ip}(:,:,iset)=ds_align{iset}{ip}(:,[1:ip]); %only include data up to max dim used
    end
    npad=dim_list_out(dptr)-ip; %pad if we are to compute a consensus in a higher dimension
    if npad>0
        z{ip}(:,ip+[1:npad],:)=0;
    end
end
%
ra.ds_knitted=cell(1,dim_list_out_max);
ra.ds_components=cell(1,nsets); 
for iset=1:nsets
    ra.ds_components{iset}=cell(1,dim_list_out_max);
end
%
%these are vector distances, taking all coordinates into account
%
ra.rmsdev_setwise=zeros(dim_list_in_max,nsets); %d1: dimension, d2: set
ra.rmsdev_stmwise=zeros(dim_list_in_max,nstims_all); %d1: dimension, d2: stim,
ra.rmsdev_overall=zeros(dim_list_in_max,1); %rms distance, across all datasets and stimuli
ra.counts_setwise=zeros(1,nsets);
ra.counts_stmwise=zeros(1,nstims_all);
%
ra.rmsdev_setwise_shuff=zeros(dim_list_in_max,nsets,1,nshuffs,2); %d1: dimension, d2: set, d3: n/a, d4: shuffle, d5: shuffle all coords or last coord
ra.rmsdev_stmwise_shuff=zeros(dim_list_in_max,nstims_all,1,nshuffs,2); %d1: dimension, d2: stim, d3: n/a, d4: shuffle, d5: shuffle all coords or last coord
ra.rmsdev_overall_shuff=zeros(dim_list_in_max,1,1,nshuffs,2); %d1: dimension, d2: n/a, d3: n/a, d4: shuffle, d5: shuffle last coord or all coords
%
%rms variance available in original data
ra.rmsavail_overall=zeros(dim_list_in_max,1);
ra.rmsavail_setwise=zeros(dim_list_in_max,nsets);
ra.rmsavail_stmwise=zeros(dim_list_in_max,nstims_all);
%
%shuffled values
ra.rmsdev_setwise_shuff=zeros(dim_list_in_max,nsets,1,nshuffs,2); %d1: dimension, d2: set, d3: n/a, d4: shuffle, d5: shuffle all coords or last coord
ra.rmsdev_stmwise_shuff=zeros(dim_list_in_max,nstims_all,1,nshuffs,2); %d1: dimension, d2: stim, d3: n/a, d4: shuffle, d5: shuffle all coords or last coord
ra.rmsdev_overall_shuff=zeros(dim_list_in_max,1,1,nshuffs,2); %d1: dimension, d2: n/a, d3: n/a, d4: shuffle, d5: shuffle last coord or all coords
%
for dptr=1:length(dim_list_in)
    ip=dim_list_in(dptr);
    dim_out=dim_list_out(dptr);
    sqs=sum(z{ip}.^2,2);
    ra.rmsavail_setwise(ip,:)=reshape(sqrt(mean(sqs,1,'omitnan')),[1 nsets]);
    ra.rmsavail_stmwise(ip,:)=reshape(sqrt(mean(sqs,3,'omitnan')),[1 nstims_all]);
    ra.rmsavail_overall(ip,:)=sqrt(mean(sqs(:),'omitnan'));
    %do unshuffled
    [consensus,znew,ts,details,ou.opts_pcon_used{ip}]=procrustes_consensus(z{ip},opts_pcon);
    if opts_pcon.if_log
        disp(sprintf(' creating Procrustes consensus for dim %2.0f based on component datasets, iterations: %4.0f, final total rms dev per coordinate: %8.5f',...
            ip,length(details.rms_change),sqrt(sum(details.rms_dev(:,end).^2))));
    end
    ra.ds_knitted{dim_out}=consensus;
    for iset=1:nsets
        ra.ds_components{iset}{dim_out}=znew(:,:,iset);
    end
    sqdevs=sum((znew-repmat(consensus,[1 1 nsets])).^2,2); %squared deviation of consensus from rotated component
    %rms deviation across each dataset, summed over coords, normalized by the number of stimuli in each dataset
    ra.rmsdev_setwise(ip,:)=reshape(sqrt(mean(sqdevs,1,'omitnan')),[1 nsets]);
    ra.counts_setwise=squeeze(sum(~isnan(sqdevs),1))';
    %rms deviation across each stimulus, summed over coords, normalized by the number of sets that include the stimulus
    ra.rmsdev_stmwise(ip,:)=reshape(sqrt(mean(sqdevs,3,'omitnan')),[1 nstims_all]);
    ra.counts_stmwise=(sum(~isnan(sqdevs),3))';
    %rms deviation across all stimuli and coords
    ra.rmsdev_overall(ip,1)=sqrt(mean(sqdevs(:),'omitnan'));
    ra.counts_overall=sum(~isnan(sqdevs(:)));
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
                ra.rmsdev_setwise_shuff(ip,:,1,ishuff,ist)=reshape(sqrt(mean(sqdevs,1,'omitnan')),[1 nsets]);
                ra.rmsdev_stmwise_shuff(ip,:,1,ishuff,ist)=reshape(sqrt(mean(sqdevs,3,'omitnan')),[1 nstims_all]);
                ra.rmsdev_overall_shuff(ip,1,1,ishuff,ist)=sqrt(mean(sqdevs(:),'omitnan'));
            end %ishuff
        end %ist
    end %nshuff>0
%

end
return
end

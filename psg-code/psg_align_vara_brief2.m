%psg_align_vara_brief2: brief script to summarize psg_align_vara_demo
%
%Run after psg_align_vara_demo.
%   tallies whether data is equal to minimum of surrogates
%   also recomputes group consensus and compares them with global consensus
% 
% See also:  PSG_ALIGN_VARA_DEMO, PSG_ALIGN_VARA_BRIEF, PROCRUSTES_CONSENSUS, PROCRUSTES.
%
disp(' ');
% reconstitute opts_pcon
if ~exist('opts_pcon') opts_pcon=struct; end
opts_pcon=filldefault(opts_pcon,'allow_reflection',1);
opts_pcon=filldefault(opts_pcon,'allow_offset',1);
opts_pcon=filldefault(opts_pcon,'allow_scale',0); 
opts_pcon=filldefault(opts_pcon,'max_niters',1000); %nonstandard max
opts_pcon.if_normscale=1;
opts_pcon.initialize_set='pca';
%
if ~exist('quantiles_summ') quantiles_summ=[.001 0.01 .05 .1 .25 .5 .75 .95]; end
%
disp(sprintf('nsets: %4.0f, nshuffs: %5.0f, if_normscale: %1.0f',results.nsets,results.nshuffs,results.if_normscale));
disp('group list');
disp(results.gp_list);
disp('tags');
disp(results.opts_multi.tags)
%
%determine which tasks are present, subjects are present, and whether this is cluster by task or cluster by subject
%
if ~exist('ts_dict')
    ts_dict=struct;
    ts_dict.task_codes={'s','w','u','b'};
    ts_dict.task_strings={'_sess','-wm1000_sess','-gm_sess','-br_sess'};
    ts_dict.subj_codes={'b','c','m','n','s','w','z'};
    ts_dict.subj_strings={'_BL','_CME','_MC','_NF','_SN','_SAW','_ZK'};
end
tasks_max=length(ts_dict.task_strings);
task_list=zeros(results.nsets,1);
subjs_max=length(ts_dict.subj_strings);
subj_list=zeros(results.nsets,1);
for iset=1:results.nsets
    underscores=find(results.sets{iset}.label=='_');
    %file names like \bgca3pt_ZK-wm1000_sess01_10
    if length(underscores)>=3
        label_short=results.sets{iset}.label(underscores(end-2):underscores(end));
        if ~contains(label_short,'-') %special case since similarity has no infix
            task_list(iset)=1;
        end
        for k=2:tasks_max
            if contains(label_short,ts_dict.task_strings{k})
                task_list(iset)=k;
            end
        end
        for k=1:subjs_max
            if contains(label_short,ts_dict.subj_strings{k})
                subj_list(iset)=k;
            end
        end
    end
end
if any(task_list==0) warning('cannot identify a task'); end
if any(subj_list==0) warning('cannot identify a subj'); end
%
group_type='';
task_match=1;
subj_match=1;
match_list=zeros(1,length(results.gp_list));
for k=1:length(results.gp_list)
    if min(task_list(results.gp_list{k}))~=max(task_list(results.gp_list{k}))
        task_match=0;
    else
        match_list(k)=task_list(results.gp_list{k}(1));
    end
    if min(subj_list(results.gp_list{k}))~=max(subj_list(results.gp_list{k}))
        subj_match=0;
    else
        match_list(k)=subj_list(results.gp_list{k}(1));
    end
end
if (task_match+subj_match)~=1 warning('cannot identify group type');end
if task_match==1
    group_type='task';
    gp_nums=task_list;
    codes=ts_dict.task_codes;
else
    group_type='subj';
    gp_nums=subj_list;
    codes=ts_dict.subj_codes;
end
disp(sprintf(' grouping type: %s, %3.0f groups',group_type,length(results.gp_list)))
%
gp_labels=cell(1,length(codes));
for ig_all=1:length(codes) %ig_all is total number of groups, not just those in this datast
    if ismember(ig_all,match_list)
        ig=find(match_list==ig_all);
        gp_labels{ig_all}=sprintf('%s %s: gp %1.0f',group_type,codes{ig_all},ig);
    else
        gp_labels{ig_all}=sprintf('%s %s:  N/A',group_type,codes{ig_all});
    end
end
%
%reassemble the datasets from each group and recompute consensus within each group
rmsdevs_gp_ov=NaN(length(codes),results.dim_max,2); %overall group, dimension, if_scaling
rmsdev_loc=NaN(results.ngps,results.dim_max,2);%; group in this dataset, dimension, if_scaling
for if_scaling=0:1
   %
   disp(' ');
   disp(sprintf('summary, allow_scale=%1.0f',if_scaling));
   disp(' ');
    rmsdev_grpwise=reshape(results.rmsdev_grpwise(:,1,1+if_scaling),[results.dim_max,1])';
    rmsdev_grpwise_shuff=reshape(results.rmsdev_grpwise_shuff(:,1,1+if_scaling,1,:),[results.dim_max,results.nshuffs])';
    pvals=sum(double(rmsdev_grpwise_shuff<=rmsdev_grpwise))/results.nshuffs;
    min_eq_count=double(min(rmsdev_grpwise_shuff,[],1)==rmsdev_grpwise);
    rmsdev_overall=results.rmsdev_overall(:,1,1+if_scaling)';
    disp('data, p-values, min_eq_count, overall rms dev')
    disp('        dimension')
    disp([[NaN 1:results.dim_max];[NaN rmsdev_grpwise];[NaN pvals];[NaN min_eq_count];[NaN rmsdev_overall]]);
    disp(' ');
    for ip=1:results.dim_max
        %for matlabs procrustes
        if if_scaling
            allow_scaling=true;
        else
            allow_scaling=false;
        end
        if opts_pcon.allow_reflection
            allow_reflection=true;
        else
            allow_reflection=false;
        end
        consensus_overall=results.ds_consensus{1+if_scaling}{ip};
        if opts_pcon.allow_offset
            consensus_overall=consensus_overall-repmat(mean(consensus_overall,1),results.nstims,1);
        end
        consensus_within=zeros(results.nstims,ip,results.ngps);
        consensus_within_xform=zeros(results.nstims,ip,results.ngps);
        zw=cell(1,results.ngps);
        for ig=1:results.ngps
            zw{ig}=zeros(results.nstims,ip,length(results.gp_list{ig}));
            for k=1:length(results.gp_list{ig})
                zw{ig}(:,:,k)=results.ds_components{1+if_scaling}{results.gp_list{ig}(k)}{ip};
            end
            consensus_within(:,:,ig)=procrustes_consensus(zw{ig},setfield(opts_pcon,'allow_scale',if_scaling));
            if opts_pcon.allow_offset
                consensus_within(:,:,ig)=consensus_within(:,:,ig)-repmat(mean(consensus_within(:,:,ig),1),results.nstims,1);
            end
            %rotate each consensus_within into the overall consensus
            [d,consensus_within_xform(:,:,ig)]=procrustes(consensus_overall,consensus_within(:,:,ig),...
                'Scaling',allow_scaling,'Reflection',allow_reflection);
        end %ip
        %now tabulate for each dimension, rms distances between
        %within-group consensus after transform, and overall consensus
        for ig=1:results.ngps
            sqdevs=sum((consensus_overall-consensus_within_xform(:,:,ig)).^2,2); %squared distance at each point
            rmsdev_loc(ig,ip,1+if_scaling)=mean(sqrt(sqdevs));
        end
        %
        %then, comcpare between groups in pairs.  Use procrustes_consensus for that,
        %since it treats the two datasets symmetrically
        %
    end %ip
    disp('rms devs of individual-group consensus from overall consensus')
    for ig_all=1:length(codes) %ig_all is total number of groups, not just those in this datast
        if ismember(ig_all,match_list)
            ig=find(match_list==ig_all);
            rmsdevs_gp_ov(ig_all,:,1+if_scaling)=rmsdev_loc(ig,:,1+if_scaling);
        end
        disp(sprintf(' %s, %s',gp_labels{ig_all},sprintf('%9.4f ',rmsdevs_gp_ov(ig_all,:,1+if_scaling))));
    end %ig_all
%
 end

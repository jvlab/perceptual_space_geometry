%psg_procrustes_task: special-purpose script to compare perceptual spaces
%via Procrustes (withoiut scaling) across subjects and paired tasks
%
%  See also: PSG_GET_COORDSETS, PROCRUSTES, PSG_TASK_LOADDATA, PROCRUSTES_CONSENSUS, PSG_DIMSTAT_TASK.
%
%define data selection and read data
%
if ~exist('dlists') dlists=struct; end
dlists=filldefault(dlists,'task_list',{'threshold','similarity','brightness','working_memory','unconstrained_grouping'}); %a subset of expt_grp
dlists=filldefault(dlists,'subj_list',{'bl','cme','mc','nf','sn','saw','zk'}); %could also add cme, saw, but these subjs are more incomplete
dlists=filldefault(dlists,'stimset_list',{'bc6pt','bcpm3pt','bgca3pt','dgea3pt'}); %could also add bc6pt, bcpm3pt, etc
%
if ~exist('paths') paths=struct; end
paths=filldefault(paths,'setup_path','../psg/psg_data/');
paths=filldefault(paths,'psg_path','../psg/psg_data/');
paths=filldefault(paths,'qform_path','../stim/');
%
if ~exist('opts_read') opts_read=struct; end
opts_read.if_auto=1;
opts_read=filldefault(opts_read,'if_log',0);
if ~exist('opts') opts=struct; end
opts=filldefault(opts,'opts_read',opts_read);
%
subjs_have_personal={'df','dt','jd','mc'}; %subjects with personalizd models
subjs_have_custavg={'bl','nf','saw','sn','zk'}; %subjects with customized averge models
qform_pref='btc_allraysfixedb_';
qform_suff='_100surrs_madj.mat';
%
[sets,ds,sas,opts_read_used,paths_used,dlists_used]=psg_task_loaddata(dlists,paths,opts);
ntasks=length(dlists.task_list);
nsubjs=length(dlists.subj_list);
nstimsets=length(dlists.stimset_list);
%
%set up transformations to analyze
%
if ~exist('opts_pcon') opts_pcon=struct(); end %for procrustes_consensus
opts_pcon=filldefault(opts_pcon,'initialize_set','pca');
opts_pcon=filldefault(opts_pcon,'allow_scale',0);
opts_pcon=filldefault(opts_pcon,'allow_offset',1);
%
if ~exist('task_pairs')
    task_pairs={'threshold','similarity';'similarity','brightness';'similarity','working_memory';'similarity','unconstrained_grouping'};
end
ntask_pairs=size(task_pairs,1);
%
task_pair_nums=zeros(ntask_pairs,2);
for itask_pair=1:ntask_pairs
    for k=1:2
        task_pair_nums(itask_pair,k)=strmatch(task_pairs{itask_pair,k},dlists.task_list,'exact');
    end
end
%
dims_avail=NaN(ntasks,nsubjs,nstimsets);
disp(' ');
nstims=NaN(1,nstimsets);
%check consistency of dimensions and stimulus counts
for itask=1:ntasks
    for isubj=1:nsubjs
        for istimset=1:nstimsets
            if ~isempty(ds{itask,isubj,istimset})
                nstims_this=size(ds{itask,isubj,istimset},1);
                dims_avail(itask,isubj,istimset)=length(ds{itask,isubj,istimset});
                if_dimok=1;
                for id=1:dims_avail(itask,isubj,istimset)
                    nstims_this=size(ds{itask,isubj,istimset}{id},1);
                    if isnan(nstims(istimset))
                        nstims(istimset)=nstims_this;
                    end
                    if_dimok=1;
                    if nstims(istimset)~=nstims_this
                        if_dimok=0;
                        disp(sprintf('for task %2.0f (%25s) subject %2.0f (%7s), stim set %2.0f (%10s), coords for dim %2.0f have wrong number of stimuli; condition removed.',...
                            itask,dlists.task_list{itask},isubj,dlists.subj_list{isubj},istimset,dlists.stimset_list{istimset},id));
                    end
                    if size(ds{itask,isubj,istimset}{id},2)~=id
                        disp(sprintf('for task %2.0f (%25s) subject %2.0f (%7s), stim set %2.0f (%10s), coords for dim %2.0f have wrong number of dimensions; condition removed.',...
                            itask,dlists.task_list{itask},isubj,dlists.subj_list{isubj},istimset,dlists.stimset_list{istimset},id));
                        if_dimok=0;
                    end
                end
                if if_dimok==0
                    ds{itask,isubj,istimset}=zeros(0);
                    dims_avail(itask,isubj,istimset)=NaN;
                end
            end
        end %stim
    end %subj
end %task
%
nsets_avail=sum(~isnan(dims_avail(:)));
max_dims_avail=min(dims_avail(:),[],'omitnan');
disp(' ');
disp(sprintf('%4.0f datasets found, out of %4.0f (%4.0f tasks x %4.0f subjects x %4.0f stimulus sets)',nsets_avail,...
    ntasks*nsubjs*nstimsets,ntasks,nsubjs,nstimsets));
disp(sprintf('maximum dimension available across all datasets: %3.0f',max_dims_avail));
disp('options used for procrustes consensus:');
disp(opts_pcon);
%
dvals=cell(ntask_pairs,nsubjs+2,nstimsets); %nsubjs+1: consensus for subjects who did both tasks, nsubjs+2: consensus for subjects who did all tasks
for istimset=1:nstimsets
    disp(sprintf(' '));
    disp(sprintf('analyzing stimulus set %1.0f: %s, %3.0f stimuli',istimset,dlists.stimset_list{istimset},nstims(istimset)));
    %determine which subjects have data for this each task pair and all task pairs
    subjs_havedata_each=cell(1,ntask_pairs);
    subjs_havedata_all=[1:nsubjs];
    for itask_pair=1:ntask_pairs
        subjs_havedata_each{itask_pair}=find(all(dims_avail(task_pair_nums(itask_pair,:),:,istimset)>0,1));  %subjects with data for this pair of conditions and stimulus set
        if ~isempty(subjs_havedata_each{itask_pair})
            subjs_havedata_all=intersect(subjs_havedata_all,subjs_havedata_each{itask_pair});
        end
    end
    ea_string=cell(1,2);
    ea_string{2}='cn all:';
    for isubj_ptr=1:length(subjs_havedata_all) %label for subjects with data for this task pair
        isubj=subjs_havedata_all(isubj_ptr);
        ea_string{2}=cat(2,ea_string{2},' ',dlists.subj_list{isubj});
    end
    %calculate consensus across subjects, for subjects with data for all tasks and for each task
    con_all=cell(1,ntasks);
    con_each=cell(ntask_pairs,ntasks);   
    for itask=1:ntasks
        if any(dims_avail(itask,:,istimset)>0)
            for itask_pair=0:ntask_pairs %0 for all, 1:ntask_pairs for each
                if itask_pair==0
                    subjs_havedata=subjs_havedata_all;
                    con_all{itask}=cell(1,max_dims_avail);
                else
                    if ~ismember(itask,task_pair_nums(itask_pair,:))
                        subjs_havedata=[];
                    else
                        subjs_havedata=subjs_havedata_each{itask_pair};
                        con_each{itask_pair,itask}=cell(1,max_dims_avail);
                    end
                end
                if ~isempty(subjs_havedata) 
                    z=cell(1,max_dims_avail);
                    for id=1:max_dims_avail %select data from the relevant subjects
                        z{id}=zeros(nstims(istimset),id,length(subjs_havedata));
                        for isubj_ptr=1:length(subjs_havedata)
                            isubj=subjs_havedata(isubj_ptr);
                            z{id}(:,:,isubj_ptr)=ds{itask,isubj,istimset}{id};
                        end
                        cons=procrustes_consensus(z{id},opts_pcon);
                        if itask_pair==0
                            con_all{itask}{id}=cons;
                        else
                            con_each{itask_pair,itask}{id}=cons;
                        end
                    end %id
                end
            end %itask_pair
        end %any data for this task?
    end %itask
    %
    for itask_pair=1:ntask_pairs
        subjs_havedata=subjs_havedata_each{itask_pair};
%        task_pairs{itask_pair,:},
%        dlists.subj_list{subjs_havedata}
        if ~isempty(subjs_havedata)
            disp(sprintf('frac of variance explained by procrustes transformation: %s to %s for %s',task_pairs{itask_pair,:},dlists.stimset_list{istimset}));
            disp(cat(2,'                 subj/dim    ',sprintf(' %7.0f',[1:max_dims_avail]))); 
            ea_string{1}='cn each:';
            for isubj_ptr=1:length(subjs_havedata) %label for subjects with data for this task pair
                isubj=subjs_havedata(isubj_ptr);
                ea_string{1}=cat(2,ea_string{1},' ',dlists.subj_list{isubj});
            end
            for isubj_ptr=1:length(subjs_havedata)+2
                if isubj_ptr<=length(subjs_havedata)
                    isubj=subjs_havedata(isubj_ptr);
                    dvals{itask_pair,isubj,istimset}=NaN(1,max_dims_avail);
                  % D = procrustes(X, Y) determines a linear transformation of the points in the matrix Y to best conform them to the points in the matrix X.
                    for id=1:max_dims_avail
                        dvals{itask_pair,isubj,istimset}(id)=procrustes(ds{task_pair_nums(itask_pair,2),isubj,istimset}{id},ds{task_pair_nums(itask_pair,1),isubj,istimset}{id},'Scaling',false);
                    end
                    subj_string=sprintf('               just  %3s        ',dlists.subj_list{isubj});
                    disp(cat(2,subj_string,sprintf(' %7.4f',dvals{itask_pair,isubj,istimset}(:))));
                else %use a consensus subject
                    ea=isubj_ptr-length(subjs_havedata); %1 for consensus across subjects with data for each task, 2 for consensus across subjects with data for all tasks
                    switch ea
                        case 1
                            for id=1:max_dims_avail
                                dvals{itask_pair,nsubjs+ea,istimset}(id)=procrustes(con_each{itask_pair,task_pair_nums(itask_pair,2)}{id},con_each{itask_pair,task_pair_nums(itask_pair,1)}{id},'Scaling',false);
                            end
                        case 2
                            for id=1:max_dims_avail
                                dvals{itask_pair,nsubjs+ea,istimset}(id)=procrustes(con_all{task_pair_nums(itask_pair,2)}{id},con_all{task_pair_nums(itask_pair,1)}{id},'Scaling',false);
                            end
                    end
                    subj_string=sprintf(' %30s ',ea_string{ea});
                    disp(cat(2,subj_string,sprintf(' %7.4f',dvals{itask_pair,nsubjs+ea,istimset}(:))));

                end %
            end %isubj_ptr
        end %have data for this pair of conditions and this stimulus set
    end
end

%psg_procrustes_task: special-purpose script to compare perceptual spaces
%via Procrustes (withoiut scaling) across subjects and paired tasks
%
%  See also: PSG_GET_COORDSETS, PROCRUSTES, PSG_TASK_LOADDATA.
%
%define data selection and read data
%
if ~exist('dlists') dlists=struct; end
dlists=filldefault(dlists,'task_list',{'threshold','similarity','brightness','working_memory','unconstrained_grouping'}); %a subset of expt_grp
dlists=filldefault(dlists,'subj_list',{'bl','cme','mc','nf','sn','saw','zk'}); %could also add cme, saw, but these subjs are more incomplete
dlists=filldefault(dlists,'stim_list',{'bc6pt','bcpm3pt','bgca3pt','dgea3pt'}); %could also add bc6pt, bcpm3pt, etc
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
nstims=length(dlists.stim_list);
%
%set up transformations to analayze
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
dims_avail=NaN(ntasks,nsubjs,nstims);
disp(' ');
for itask=1:ntasks
    for isubj=1:nsubjs
        for istim=1:nstims
            if ~isempty(ds{itask,isubj,istim})
                dims_avail(itask,isubj,istim)=length(ds{itask,isubj,istim});
                if_dimok=1;
                for id=1:dims_avail(itask,isubj,istim)
                    if size(ds{itask,isubj,istim}{id},2)~=id
                        disp(sprintf('for task %2.0f (%25s) subject %2.0f (%7s), stim set %2.0f (%10s), coords for dim %2.0f have wrong number of dimensions; condition removed.',...
                            itask,dlists.task_list{itask},isubj,dlists.subj_list{isubj},istim,dlists.stim_list{istim},id));
                        if_dimok=0;
                    end
                end
                if (if_dimok)
                  
                else
                    ds{itask,isubj,istim}=zeros(0);
                    dims_avail(itask,isubj,istim)=NaN;
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
    ntasks*nsubjs*nstims,ntasks,nsubjs,nstims));
disp(sprintf('maximum dimension available across all datasets: %3.0f',max_dims_avail));
%
dvals=cell(ntask_pairs,nsubjs,nstims);
for istim=1:nstims
    disp(sprintf(' '));
    disp(sprintf('analyzing stimulus set %1.0f -> %s',istim,dlists.stim_list{istim}));
    %determine which subjects have data for this each task pair and all task pairs
    subjs_havedata_each=cell(1,ntask_pairs);
    subjs_havedata_all=[1:nsubjs];
    for itask_pair=1:ntask_pairs
        subjs_havedata_each{itask_pair}=find(all(dims_avail(task_pair_nums(itask_pair,:),:,istim)>0,1));  %subjecs with data for this pair of conditions and stimulus set
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
    for itask_pair=1:ntask_pairs
        subjs_havedata=subjs_havedata_each{itask_pair};
%        task_pairs{itask_pair,:},
%        dlists.subj_list{subjs_havedata}
        if ~isempty(subjs_havedata)
            npts=size(ds{task_pair_nums(itask_pair,1),subjs_havedata(1),istim}{1},1);
            disp(sprintf('transformation: %s to %s for %s',task_pairs{itask_pair,:},dlists.stim_list{istim}));
            disp(cat(2,'         subj/dim         ',sprintf(' %7.0f',[1:max_dims_avail]))); 
            ea_string{1}='cn each:';
            for isubj_ptr=1:length(subjs_havedata) %label for subjects with data for this task pair
                isubj=subjs_havedata(isubj_ptr);
                ea_string{1}=cat(2,ea_string{1},' ',dlists.subj_list{isubj});
            end
            for isubj_ptr=1:length(subjs_havedata)+2
                if isubj_ptr<=length(subjs_havedata)
                isubj=subjs_havedata(isubj_ptr);
                dvals{itask_pair,isubj,istim}=NaN(1,max_dims_avail);
              % D = procrustes(X, Y) determines a linear transformation of the points in the matrix Y to best conform them to the points in the matrix X.
                for id=1:max_dims_avail
                    dvals{itask_pair,isubj,istim}(id)=procrustes(ds{task_pair_nums(itask_pair,2),isubj,istim}{id},ds{task_pair_nums(itask_pair,1),isubj,istim}{id},'Scaling',false);
                end
                subj_string=sprintf('           %3s             ',dlists.subj_list{isubj});
                disp(cat(2,subj_string,sprintf(' %7.4f',dvals{itask_pair,isubj,istim}(:))));
                else
                    ea=isubj_ptr-length(subjs_havedata); %1 for consensus across subjects with data for each task, 2 for consensus across subjects with data for all tasks
                    subj_string=sprintf(' %25s ',ea_string{ea});
                    disp(cat(2,subj_string,sprintf(' %7.4f',zeros(1,max_dims_avail))));
                end 
            end %isubj_ptr
        end %have data for this pair of conditions and this stimulus set
    end
end

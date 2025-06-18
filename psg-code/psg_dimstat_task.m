%psg_dimstat_task: special-purpose script to evaluate models of various dimensions
%compare perceptual spaces based on Künstle, von Luxburg, and Wichmann, JOV 2022
%
%  See also: PSG_GET_COORDSETS, PSG_TASK_LOADDATA, PSG_PROCRUSTES_TASK.
%
%define data selection and read data
%
if ~exist('dlists') dlists=struct; end
dlists=filldefault(dlists,'task_list',{'similarity','brightness','working_memory','unconstrained_grouping'}); %a subset of expt_grp
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
opts.if_choices=1; %also read choice data
%
[sets,ds,sas,opts_read_used,paths_used,dlists_used,choices]=psg_task_loaddata(dlists,paths,opts);
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
                if isempty(choices{itask,isubj,istimset})
                    disp(sprintf('for task %2.0f (%25s) subject %2.0f (%7s), stim set %2.0f (%10s), choices file missing; condition removed.',...
                        itask,dlists.task_list{itask},isubj,dlists.subj_list{isubj},istimset,dlists.stimset_list{istimset}));
                        if_dimok=0;
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
%disp('options used for procrustes consensus:');
%disp(opts_pcon);
for istimset=1:nstimsets
    disp(sprintf(' '));
    disp(sprintf('analyzing stimulus set %1.0f: %s, %3.0f stimuli',istimset,dlists.stimset_list{istimset},nstims(istimset)));
    for itask=1:ntasks
        subjs_havedata=find(all(dims_avail(itask,:,istimset)>0,1)); 
        disp(sprintf('task %25s: %3.0f subjects',dlists.task_list{itask},length(subjs_havedata)));
    end %task
end %subject

%psg_dimstat_task: special-purpose script to evaluate models of various dimensions
%compare perceptual spaces based on Künstle, von Luxburg, and Wichmann, JOV 2022
%
%  See also: PSG_GET_COORDSETS, PSG_TASK_LOADDATA, PSG_PROCRUSTES_TASK, COOTODSQ.
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
acc=nan(max_dims_avail,ntasks,nstimsets,nsubjs,nsubjs); %d1: dimension, d2: task, d3: stimulus set, d4: subject for model, d5: subject for prediction
for istimset=1:nstimsets
    disp(sprintf(' '));
    disp(sprintf('analyzing stimulus set %1.0f: %s, %3.0f stimuli',istimset,dlists.stimset_list{istimset},nstims(istimset)));
    for itask=1:ntasks
        subjs_havedata=find(all(dims_avail(itask,:,istimset)>0,1));
        disp(sprintf('task %25s: %3.0f subjects',dlists.task_list{itask},length(subjs_havedata)));
        if length(subjs_havedata)>=2
            for imod_ptr=1:length(subjs_havedata)
                imod=subjs_havedata(imod_ptr);
                for id=1:max_dims_avail
                    distsq=cootodsq(ds{itask,imod,istimset}{id}); %compute distances squared, suffices for comparisons
                    for iprd_ptr=1:length(subjs_havedata)
                        iprd=subjs_havedata(iprd_ptr);
                        ch=choices{itask,iprd,istimset};
                        nch=size(ch,1); %number of choices
                        nacc=0; % number of accurate predictions
                        ncols=size(ch,2);
                        switch ncols
                            case 5 %triadic comparisons: columns are [ref s1 s2 choices(d(ref,s1)<d(ref,s2)) total trials]
                                %vectorized referencing into distsq
                                d1ptr=ch(:,1)+nstims(itask)*(ch(:,2)-1);
                                d2ptr=ch(:,1)+nstims(itask)*(ch(:,3)-1);
                            case 6 %tetradic comparisons: columns are [s1  s2 s3  s4 choices(d(s1,s2)<d(s3, s4)) total trials]
                                d1ptr=ch(:,1)+nstims(itask)*(ch(:,2)-1);
                                d2ptr=ch(:,3)+nstims(itask)*(ch(:,4)-1);
                            otherwise 
                                warning(sprintf('unexpected format for choice data for subject %s',dlists.subj_list{iprd}));
                        end
                        d1sq=distsq(d1ptr);
                        d2sq=distsq(d2ptr);
                        nacc=sum(sign(d2sq(:)-d1sq(:)).*ch(:,ncols-1))+sum(ch(d1sq(:)>=d2sq(:),ncols));
                        acc(id,itask,istimset,imod,iprd)=nacc/sum(ch(:,ncols));
                    end %subj for predictin
                end %dimension
            end %subj for model data
        end %need at least 2 subjects
    end %task
end %subject
acc_self=zeros(max_dims_avail,ntasks,nstimsets); %mean accuracy of prediction within a subject
acc_xsub=zeros(max_dims_avail,ntasks,nstimsets); %mean accuracy of cross-subject prediction
for istimset=1:nstimsets
    for itask=1:ntasks
        for id=1:max_dims_avail
            acc_resh=reshape(acc(id,istimset,itask,:,:),[nsubjs,nsubjs]);
            acc_self(id,itask,istimset)=mean(diag(acc_resh),'omitnan');
            acc_xsub(id,itask,istimset)=mean(mean(triu(acc_resh)+tril(acc_resh)-diag(diag(acc_resh)),'omitnan'),'omitnan');
        end
    end
end

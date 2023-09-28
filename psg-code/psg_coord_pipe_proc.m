%psg_pipe_coord_proc reads one or more coordinate datasets and processes, either by forming a consensus
% or simple transformations
%
% The output datafile has the fields
%   dim[a], where [a] is a character string 1,...,7,8,9,10,
%   stim_labels, as in the coord data files
%   pipeline, indicating the processing pipeline to create this file
%    pipeline{n} indicates the most recent processing step
%
% all datasets must have dimension lists beginning at 1 and without gaps
%
%  See also: PSG_GET_COORDSETS, PSG_QFORM2COORD_PROC, PSG_READ_COORDDATA, PSG_WTITE_COORDDATA.
%
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_pcon') opts_pcon=struct(); end %for procrustes_consensus
%
pipe_types_long={'done','create a consensus dataset','Procrustes alignment to a given set'};
pipe_types={'','consensus','Procrustes'};
npipe_types=length(pipe_types)-1;
%
disp('This will process one or more coordinate datasets and create new coordinate datasets, via simple transformations or consensus.');
%
nsets=getinp('number of datasets','d',[1 100]);
opts_read.input_type=1;
[sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,[],nsets);
%
%eventually: consensus, and/or align via Procrustes, and/or align via PCA
%
nsets_orig=nsets;
nstims_each=zeros(1,nsets);
labels=cell(1,nsets);
proc_strings=cell(1,nsets);
dim_list_common=[];
for iset=1:nsets
    labels{iset}=sets{iset}.label;
    dim_list=sets{iset}.dim_list;
    nstims_each(iset)=size(ds{1}{1},1);
    disp(sprintf('dataset %2.0f (%30s): %2.0f stimuli, dimension list is',iset,labels{iset},nstims_each(iset)));
    disp(dim_list);
    if iset==1
        dim_list_common=dim_list;
    else
        dim_list_common=intersect(dim_list_common,dim_list);
    end
end
dim_max=max(dim_list_common);
if any(dim_list_common~=[1:length(dim_list_common)]);
    warning('at least one dimension list has a gap');
    if_done=1;
else
    if_done=0;
end
if (nsets>1)
    while(if_done==0)
        for iset=1:nsets
            disp(sprintf('dataset %2.0f is %s',iset,labels{iset}));
        end
        for ipipe=0:npipe_types
            disp(sprintf('%1.0f->%s',ipipe,pipe_types_long{ipipe+1}));
        end
        pipe_type=getinp('choice','d',[0 npipe_types]);
        if pipe_type>0
            switch pipe_types{pipe_type+1}
                %%%%%%%%%%%need to:
                %use psg_write_coorddata at end               
                case 'consensus'
                    con_sets=getinp('datasets to combine','d',[1 nsets],[1:nsets_orig]);
                    con_sets_string=sprintf('%1.0f ',con_sets);
                    if min(nstims_each(con_sets))~=max(nstims_each(con_sets)) %should have been blocked by psg_read_coordsets
                        disp('number of stimuli do not match.');
                    else
                        nstims=min(nstims_each(con_sets));
                        allow_scale=getinp('1 to allow scaling','d',[0 1]);
                        z=cell(1,dim_max);
                        znew=cell(1,dim_max); %datasets after transformatoin
                        consensus=cell(1,dim_max); %consensus
                        ts=cell(1,dim_max); %transformations
                        file_list=cell(0);
                        for id=1:dim_max
                            z{id}=zeros(nstims,id,length(con_sets));
                            for iset_ptr=1:length(con_sets)
                                iset=con_sets(iset_ptr);
                                z{id}(:,:,iset_ptr)=ds{iset}{id};
                                file_list{iset_ptr}=labels{iset};
                            end
                            [consensus{id},znew{id},ts{id},details,opts_pcon_used]=procrustes_consensus(z{id},setfield(opts_pcon,'allow_scale',allow_scale));
                            %    znew(:,:,iset)=ts{iset}.scaling*z(:,:,iset)*ts{iset}.orthog+repmat(ts{iset}.translation,npts,1);
                            disp(sprintf('calculated consensus for %3.0f dimensions; iterations: %4.0f',id,length(details.ts_cum)));
                        end
                        %create metadata for each old dataset rotated into consensus
                        for iset_ptr=1:length(con_sets)
                            nsets=nsets+1; %a new dataset for each file rotated into consensus
                            for id=1:dim_max
                                ds{nsets}{id}=znew{id}(:,:,iset_ptr);
                            end
                            iset=con_sets(iset_ptr);
                            labels{nsets}=sprintf('set %1.0f rotated into consensus of sets %s, allow_scale=%1.0f',iset,con_sets_string,allow_scale);
                            nstims_each(nsets)=nstims;
                            proc_strings{nsets}.type='rotate_to_consensus';
                            proc_strings{nsets}.opts.opts_pcon_used=opts_pcon_used;
                            proc_strings{nsets}.file_list=file_list;
                        end
                        %create medatata for consensus
                        nstims_each(nsets)=nstims;
                        ds{nsets}=consensus;
                        labels{nsets}=sprintf('consensus of sets %s, allow_scale=%1.0f',con_sets_string,allow_scale);
                        proc_strings{nsets}.type='consensus';
                        proc_strings{nsets}.opts.opts_pcon_used=opts_pcon_used;
                        proc_strings{nsets}.file_list=file_list;
                    end
                case 'Procrustes'
            end
        else
            if_done=1;
        end
    end
end
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%
% ndim_max=length(ds{1});
% %
% %assemble the output file
% %
% sout=struct;
% sout.stim_labels=strvcat(sas{1}.typenames);
% sout.pipeline{1}.type='qform2coord';
% sout.pipeline{1}.sets=sets{1};
% sout.pipeline{1}.opts.opts_read_used=opts_read_used{1};
% sout.pipeline{1}.opts.opts_qpred_used=opts_qpred_used{1};
% for idim=1:ndim_max
%     dname=cat(2,'dim',sprintf('%1.0f',idim));
%     sout.(dname)=ds{1}{idim};
% end
% disp(sprintf('output structure with %2.0f stimuli and up to %2.0f dimensions created via %s',size(sout.stim_labels,1),ndim_max,sets{1}.label));
% %
% setup_name=opts_read_used{1}.setup_fullname;
% setup_name=strrep(setup_name,'9.mat','');
% setup_name=strrep(setup_name,'.mat','');
% fname_suggest=strrep(fname_suggest_base,'*',setup_name);
% fname=getinp(sprintf('file name, e.g., %s',fname_suggest_base),'s',[],fname_suggest);
% save(fname,'-struct','sout');
% disp(sprintf('%s written.',fname));
% 

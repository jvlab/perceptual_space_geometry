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
%  See also: PSG_GET_COORDSETS, PSG_QFORM2COORD_PROC, PSG_READ_COORDDATA, PSG_WRITE_COORDDATA.
%
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_pcon') opts_pcon=struct(); end %for procrustes_consensus
if ~exist('opts_write') opts_write=struct(); end %for psg_write_coorddata
%
pipe_types_long={'done (and proceed to write data files)','create a consensus dataset','Procrustes alignment to a given set'};
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
pipelines=cell(1,nsets);
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
                        sets_combined=cell(0);
                        for id=1:dim_max
                            z{id}=zeros(nstims,id,length(con_sets));
                            for iset_ptr=1:length(con_sets)
                                iset=con_sets(iset_ptr);
                                z{id}(:,:,iset_ptr)=ds{iset}{id};
                                file_list{iset_ptr}=labels{iset};
                                sets_combined{iset_ptr}=sets{iset};
                            end
                            [consensus{id},znew{id},ts{id},details,opts_pcon_used]=procrustes_consensus(z{id},setfield(opts_pcon,'allow_scale',allow_scale));
                            disp(sprintf('calculated consensus for %3.0f dimensions; iterations: %4.0f',id,length(details.ts_cum)));
                        end
                        %create data and metadata for each old dataset rotated into consensus
                        for iset_ptr=1:length(con_sets)
                            nsets=nsets+1; %a new dataset for each file rotated into consensus
                            for id=1:dim_max
                                ds{nsets}{id}=znew{id}(:,:,iset_ptr);
                            end
                            iset=con_sets(iset_ptr);
                            sas{nsets}=sas{iset}; %take previous stimulus descriptors
                            labels{nsets}=sprintf('set %1.0f rotated into consensus of sets %s, allow_scale=%1.0f',iset,con_sets_string,allow_scale);
                            nstims_each(nsets)=nstims;
                            pipelines{nsets}.type='rotate_to_consensus';
                            pipelines{nsets}.opts.opts_pcon_used=opts_pcon_used;
                            pipelines{nsets}.files_combined=file_list;
                            pipelines{nsets}.sets_combined=sets_combined;
                            pipelines{nsets}.sets=sets{iset};
                        end
                        %create medatata for consensus
                        nsets=nsets+1; %a new dataset for each file rotated into consensus
                        nstims_each(nsets)=nstims;
                        ds{nsets}=consensus;
                        sas{nsets}=sas{con_sets(1)} %assume all stimulus descriptors are the same
                        labels{nsets}=sprintf('consensus of sets %s, allow_scale=%1.0f',con_sets_string,allow_scale);
                        pipelines{nsets}.type='consensus';
                        pipelines{nsets}.opts.opts_pcon_used=opts_pcon_used;
                        pipelines{nsets}.files_combined=file_list;
                        pipelines{nsets}.sets_combined=sets_combined;
                    end
                case 'Procrustes'
            end
        else
            if_done=1;
        end
    end
end
%write files if asked
for iset=1:nsets_orig
    disp(sprintf('original dataset %2.0f is %s',iset,labels{iset}));
end
for iset=nsets_orig+1:nsets
    disp(sprintf('computed dataset %2.0f is %s',iset,labels{iset}));
end
if_write=-1;
while (if_write~=0)
    if_write=getinp(sprintf('1 to %2.0f to write a file, 0 to end',nsets),'d',[0 nsets]);
    if (if_write>0)
        iset=if_write;
        data_fullname_def=[];
        sout=struct;
        sout.stim_labels=strvcat(sas{iset}.typenames);
        sout.pipeline{1}=pipelines{iset}; %might need to change this to keep track of multiple operations
        opts_used=psg_write_coorddata([],ds{iset},sout,setfield(opts_write,'data_fullname_def',data_fullname_def));
    end
end

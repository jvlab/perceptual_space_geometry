%psg_pipe_coord_proc reads one or more coordinate datasets from experiments 
% processes, either by forming a consensus or simple transformations,
% and simple plots
%
% PCA rotation allows for PC's to be determined based on a subset of stimuli
%  grouped by dimension -- e.g., the first two pc's maximize the variance
%  due to stimuli that contain 'b', and the next two pc's maximize the
%  remaining variance due to stimuli that contain 'g'
%
% The output datafile has the fields
%   dim[a], where [a] is a character string 1,...,7,8,9,10,
%   stim_labels, as in the coord data files
%   pipeline, indicating the processing pipeline to create this file
%    pipeline{n} indicates the most recent processing step
%    pipeline{n}.type: a string, e.g., 'rotate_to_consensus,'consensus','pca_rotation' 
%    pipeline{n}.sets: the 'set' field of a single data file that was used
%          to construct this dataset [may be absent if dataset is a consensus]
%    pipeline{n}.sets_combined:
%          cell array of the 'set' fields of the datasets combined to create this dataset
%    pipeline{n}.file_list:
%          cell array of the 'labels' fields of the datasets combined to create this dataset
%    pipeline{n}.opts: options used
%  for type='rotate_to_consensus', sets is the file that is rotated to the consensus,
%    while sets_combined, file_list indicate the collection of files used to form the consensus
%  for type='pca_rotation', sets_combined' and file_list' are not present, since
%    processing only uses a single dataset
%  for type='consensus', only sets_combined and file_list are present,
%    since all of those files are used together
% 
% To do: Procrustes alignment
%
% all datasets must have dimension lists beginning at 1 and without gaps
%
% 18Feb24: allow for propagtion of pipeline; show available files after each processing step
%
%  See also: PSG_GET_COORDSETS, PSG_QFORM2COORD_PROC, PSG_READ_COORDDATA, PSG_WRITE_COORDDATA, PSG_PLOTCOORDS,
%    SVD.
%
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_pcon') opts_pcon=struct(); end %for procrustes_consensus
if ~exist('opts_write') opts_write=struct(); end %for psg_write_coorddata
if ~exist('opts_plot') opts_plot=struct(); end % for psg_plotcoords
%
pipe_types_long={'done (and proceed to write data files)',...
    'create a consensus dataset from datasets with matching entries',...
    'pca rotation of individual dataset(s)',...
    'procrustes alignment to a given set',...
    'plot coordinates of individual dataset(s) [no new datasets created]'};
pipe_types={'','consensus','pca_rotation','procrustes','plot_coords'};
pipe_minsets=[2 1 1 1];
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
    pipelines{iset}=sets{iset}.pipeline;
    disp(sprintf('dataset %2.0f (%30s): %2.0f stimuli, pipeline length so far: %2.0f, dimension list is ',iset,labels{iset},nstims_each(iset),length(pipelines{iset})));
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
%if (nsets>1)
    while(if_done==0)
        for iset=1:nsets
            disp(sprintf('dataset %2.0f (%30s): %2.0f stimuli, pipeline length so far: %2.0f',iset,labels{iset},nstims_each(iset),length(pipelines{iset})));
            %
            if iset>length(sets)
                sets{iset}.type='processed';
                sets{iset}.dim_list=[1:dim_max];
                sets{iset}.nstims=nstims;
                sets{iset}.label=labels{iset};
                sets{iset}.pipeline=struct(); %will be added later
            end
        end
        for ipipe=0:npipe_types
            disp(sprintf('%1.0f->%s',ipipe,pipe_types_long{ipipe+1}));
        end
        pipe_type=getinp('choice','d',[0 npipe_types]);
        if pipe_type>0
            proc_sets=getinp('datasets to process','d',[1 nsets],[1:nsets_orig]);
            proc_sets_string=sprintf('%1.0f ',proc_sets);
            nproc_sets=length(proc_sets);
            if min(nstims_each(proc_sets))~=max(nstims_each(proc_sets)) %should have been blocked by psg_read_coordsets
                disp('number of stimuli do not match.');
            elseif pipe_minsets(pipe_type)>nproc_sets
                disp(sprintf('need at least %1.0f sets',pipe_minsets(pipe_type)));
            else
                %bookkeeping common to several options
                nstims=min(nstims_each(proc_sets));
                file_list=cell(0);
                for iset_ptr=1:nproc_sets
                    iset=proc_sets(iset_ptr);
                    file_list{iset_ptr}=labels{iset};
                    sets_combined{iset_ptr}=sets{iset};
                end
                switch pipe_types{pipe_type+1}
                    case 'consensus'
                        allow_scale=getinp('1 to allow scaling','d',[0 1]);
                        z=cell(1,dim_max);
                        znew=cell(1,dim_max); %datasets after transformation
                        consensus=cell(1,dim_max); %consensus
                        ts=cell(1,dim_max); %transformations
                        for id=1:dim_max
                            z{id}=zeros(nstims,id,nproc_sets);
                            for iset_ptr=1:nproc_sets
                                iset=proc_sets(iset_ptr);
                                z{id}(:,:,iset_ptr)=ds{iset}{id};
                            end
                            [consensus{id},znew{id},ots{id},details,opts_pcon_used]=procrustes_consensus(z{id},setfield(opts_pcon,'allow_scale',allow_scale));
                            disp(sprintf('calculated consensus for %3.0f dimensions; iterations: %4.0f',id,length(details.ts_cum)));
                        end
                        %create data and metadata for each old dataset rotated into consensus
                        for iset_ptr=1:nproc_sets
                            nsets=nsets+1; %a new dataset for each file rotated into consensus
                            for id=1:dim_max
                                ds{nsets}{id}=znew{id}(:,:,iset_ptr);
                            end
                            iset=proc_sets(iset_ptr);
                            sas{nsets}=sas{iset}; %take previous stimulus descriptors
                            labels{nsets}=sprintf('set %1.0f rotated into consensus of sets %s, allow_scale=%1.0f',iset,proc_sets_string,allow_scale);
                            nstims_each(nsets)=nstims;
                            nextpipe=length(pipelines{iset})+1;
                            pipelines{nsets}{nextpipe}.type='rotate_to_consensus';
                            pipelines{nsets}{nextpipe}.opts.opts_pcon_used=opts_pcon_used;
                            pipelines{nsets}{nextpipe}.files_combined=file_list;
                            pipelines{nsets}{nextpipe}.sets_combined=sets_combined;
                            pipelines{nsets}{nextpipe}.sets=sets{iset};
                        end
                        %create medatata for consensus
                        nsets=nsets+1; %a new dataset for each file rotated into consensus
                        nstims_each(nsets)=nstims;
                        ds{nsets}=consensus;
                        sas{nsets}=sas{proc_sets(1)} %assume all stimulus descriptors are the same
                        labels{nsets}=sprintf('consensus of sets %s, allow_scale=%1.0f',proc_sets_string,allow_scale);
                        nstims_each(nsets)=nstims;
                        nextpipe=1; %de novo dataset
                        pipelines{nsets}{nextpipe}.type='consensus';
                        pipelines{nsets}{nextpipe}.opts.opts_pcon_used=opts_pcon_used;
                        pipelines{nsets}{nextpipe}.files_combined=file_list;
                        pipelines{nsets}{nextpipe}.sets_combined=sets_combined;
                    case 'pca_rotation' 
                        ifok=0;
                        while (ifok==0)
                            dim_groups=[];
                            tstring_dimspecs='';
                            ndgs=0;
                            typenames_sel=cell(1,nproc_sets);
                            typenames_inds=cell(1,nproc_sets);
                            tstring_dimspec=cell(0);
                            sel=cell(0);
                            while sum(dim_groups)<dim_max
                                ndgs=ndgs+1;
                                dim_groups(ndgs)=getinp(sprintf('size of dimension group %1.0f, which starts at dimension %1.0f',ndgs,1+sum(dim_groups)),...
                                    'd',1+[0 dim_max-sum(dim_groups)],dim_max-sum(dim_groups));
                                sel{ndgs}=getinp('selection string or multiple strings, separated by |','s',[],'|');
                                tstring_dimspec{ndgs}=cat(2,'d[',sprintf('%1.0f',[(1+sum(dim_groups(1:end-1))):sum(dim_groups)]),']:',sel{ndgs});
                                tstring_dimspecs=cat(2,tstring_dimspecs,tstring_dimspec{ndgs},' ');
                                for iset_ptr=1:nproc_sets
                                    iset=proc_sets(iset_ptr);                   
                                    [typenames_sel{iset_ptr}{ndgs},typenames_inds{iset_ptr}{ndgs}]=psg_select_util(sel{ndgs},sas{iset});
                                    disp(sprintf('dataset %2.0f (%30s): %2.0f stimuli, %2.0f stimuli selected',...
                                        iset,labels{iset},nstims_each(iset),length(typenames_sel{iset_ptr}{ndgs})));
                                    disp(sprintf('typenames selected for dimension group %1.0f is %s',ndgs,tstring_dimspec{ndgs}));
                                    disp(typenames_sel{iset_ptr}{ndgs});
                                end
                            end
                            tstring_dimspecs=deblank(tstring_dimspecs);
                            disp(sprintf('dimension groups: %s',tstring_dimspecs))
                            ifok=getinp('1 if ok','d',[0 1],0);
                        end
                        %
                        if_center=getinp('1 to center the data','d',[0 1],1);
                        for iset_ptr=1:nproc_sets
                            iset=proc_sets(iset_ptr); 
                            nsets=nsets+1; %a new dataset for each file rotated into consensus
                            disp(sprintf('processing set %2.0f',iset));
                            for id=1:dim_max
                                d_unrot=ds{iset}{id};
                                if (if_center)
                                    d_unrot=d_unrot-repmat(mean(d_unrot,1),nstims,1);
                                end
                                %for each dimension group, determine directions of greatest variance for the
                                %requested stimuli, that are orthogonal to previously-determined directions
                                d_rot=d_unrot;
                                for idg=1:ndgs
                                    d_avail=[sum(dim_groups(1:idg-1))+1:id]; %dimensions available to rotate
                                    if ~isempty(d_avail) ...
                                        d_fill=sum(dim_groups(1:idg-1))+1:min(sum(dim_groups(1:idg)),id); %the dimensions to be modified
                                        d_use=d_rot(typenames_inds{iset_ptr}{idg},d_avail); %data at selected coordinates and stimuli
                                        [u,s,v]=svd(d_use); %d_use=u*s*v', with u and v both orthogonal, so u*s=d_use*v
                                        %now use v to rotate the entire dataset, via the rotation defined by the selected coordinates
                                        d_rot(:,d_avail)=d_rot(:,d_avail)*v;
                                        disp(sprintf('  for dataset dimension %1.0f: rotated coords %1.0f to %1.0f using selection %s',...
                                            id,min(d_avail),max(d_avail),sel{idg}));
                                    end
                                end
                                ds{nsets}{id}=d_rot; %coords for rotated dataset
                            end
                            %create metadata for rotated dataset:  the new set is derived from iset; there is no combination
                            sas{nsets}=sas{iset};
                            labels{nsets}=sprintf('set %1.0f, pca, centering=%1.0f, dimspecs %s',iset,if_center,tstring_dimspecs);
                            nstims_each(nsets)=nstims;
                            nextpipe=length(pipelines{iset})+1;
                            pipelines{nsets}{nextpipe}.type='pca_rotation';
                            pipelines{nsets}{nextpipe}.opts.pca_rotation.if_center=if_center;
                            pipelines{nsets}{nextpipe}.opts.pca_rotation.dim_groups=dim_groups;
                            pipelines{nsets}{nextpipe}.opts.pca_rotation.sel=sel;
                            pipelines{nsets}{nextpipe}.opts.pca_rotation.typenames_sel=typenames_sel{iset_ptr};
                            pipelines{nsets}{nextpipe}.opts.pca_rotation.typenames_inds=typenames_inds{iset_ptr};
                            pipelines{nsets}{nextpipe}.sets=sets{iset};
                        end %next iset_ptr to rotate
                    case 'procrustes'
                        disp('not implemented yet');
                    case 'plot_coords'
                        dim_plot=getinp('dimension to plot','d',[2 dim_max]);
                        dim_select=getinp('dimension selection','d',[1 dim_plot],[1:min(dim_plot,4)]);
                        opts_plot_used=cell(1,nsets);
                        for iset_ptr=1:nproc_sets
                            iset=proc_sets(iset_ptr);
                            tstring_name=sprintf('dataset %2.0f (%30s)',iset,labels{iset});
                            tstring_dim=cat(2,sprintf('dim %1.0f',dim_plot),' [',sprintf(' %1.0f',dim_select),' ]');
                            %
                            figure;
                            set(gcf,'Position',[50 50 1200 800]);
                            set(gcf,'NumberTitle','off');
                            set(gcf,'Name',cat(2,tstring_name,' ',tstring_dim));
                            opts_plot_use=opts_plot;
                            opts_plot_use.axis_handle=gca;
                            if (iset<=nsets_orig)
                                rays_use=rayss{iset};
                            else %create a ray structure if possible
                                if isfield(sas{iset},'btc_specoords') %check if btc_specoords is present and if not use a noray option
                                    rays_use=psg_findrays(sas{iset}.btc_specoords,opts_rays);
                                else
                                    rays_use=[];
                                end
                            end
                            if isempty(rays_use)
                                opts_plot_use.if_use_rays=0;
                            end
                            opts_plot_used{iset}=psg_plotcoords(ds{iset}{dim_plot},dim_select,sas{iset},rays_use,opts_plot_use);
                            axes('Position',[0.01,0.02,0.01,0.01]); %for text
                            text(0,0,cat(2,tstring_name,' ',tstring_dim),'Interpreter','none','FontSize',8);
                            axis off;
                        end %next iset_ptr to plot
                end %switch
            end
        else %pipe_type>0
            if_done=1;
        end
    end %if_done
%end
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
        %propagate pipeline from previous processing steps
        sout.pipeline=pipelines{iset};
        opts_used=psg_write_coorddata([],ds{iset},sout,setfield(opts_write,'data_fullname_def',data_fullname_def));
    end
end

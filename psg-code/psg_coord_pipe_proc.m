%psg_pipe_coord_proc reads one or more coordinate datasets from experiments 
% processes, either by forming a consensus or simple transformations,
% and plots
%
% The output datafile has the fields
%   dim[a], where [a] is a character string 1,...,7,8,9,10,
%   stim_labels, as in the coord data files
%   pipeline, indicating the processing pipeline to create this file
%    pipeline{n} indicates the most recent processing step
%
% To do: Procrustes alignment, centering, PCA into the subspace closest to specific stimuli
%
% all datasets must have dimension lists beginning at 1 and without gaps
%
% 18Feb24: allow for propagtion of pipeline; show available files after each processing step
%
%  See also: PSG_GET_COORDSETS, PSG_QFORM2COORD_PROC, PSG_READ_COORDDATA, PSG_WRITE_COORDDATA, PSG_PLOTCOORDS.
%
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_pcon') opts_pcon=struct(); end %for procrustes_consensus
if ~exist('opts_write') opts_write=struct(); end %for psg_write_coorddata
if ~exist('opts_plot') opts_plot=struct(); end % for psg_plotcoords
%
pipe_types_long={'done (and proceed to write data files)',...
    'create a consensus dataset from datasets with matching entries',...
    'PCA rotation of individual dataset(s)',...
    'Procrustes alignment to a given set',...
    'plot coordinates of individual dataset(s) [no new datasets created]'};
pipe_types={'','consensus','PCA_rotation','Procrustes','plot_coords'};
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
if (nsets>1)
    while(if_done==0)
        for iset=1:nsets
            disp(sprintf('dataset %2.0f (%30s): %2.0f stimuli, pipeline length so far: %2.0f',iset,labels{iset},nstims_each(iset),length(pipelines{iset})));
        end
        for ipipe=0:npipe_types
            disp(sprintf('%1.0f->%s',ipipe,pipe_types_long{ipipe+1}));
        end
        pipe_type=getinp('choice','d',[0 npipe_types]);
        if pipe_type>0
            proc_sets=getinp('datasets to process','d',[1 nsets],[1:nsets_orig]);
            proc_sets_string=sprintf('%1.0f ',proc_sets);
            nproc_sets=length(proc_sets);
            switch pipe_types{pipe_type+1}
                case 'consensus'
                     if min(nstims_each(proc_sets))~=max(nstims_each(proc_sets)) %should have been blocked by psg_read_coordsets
                        disp('number of stimuli do not match.');
                    else
                        nstims=min(nstims_each(proc_sets));
                        allow_scale=getinp('1 to allow scaling','d',[0 1]);
                        z=cell(1,dim_max);
                        znew=cell(1,dim_max); %datasets after transformation
                        consensus=cell(1,dim_max); %consensus
                        ts=cell(1,dim_max); %transformations
                        file_list=cell(0);
                        sets_combined=cell(0);
                        for id=1:dim_max
                            z{id}=zeros(nstims,id,nproc_sets);
                            for iset_ptr=1:nproc_sets
                                iset=proc_sets(iset_ptr);
                                z{id}(:,:,iset_ptr)=ds{iset}{id};
                                file_list{iset_ptr}=labels{iset};
                                sets_combined{iset_ptr}=sets{iset};
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
                    end
                case 'Procrustes'
                    disp('not implemented yet');
                case 'PCA_rotation'
                    ifok=0;
                    while (ifok==0)
                        sel=getinp('selection string or multiple strings, separated by |','s',[]);
                        typenames_sel=cell(nproc_sets,1);
                        typenames_inds=cell(nproc_sets,1);
                        for iset_ptr=1:nproc_sets
                            iset=proc_sets(iset_ptr);                   
                            [typenames_sel{iset_ptr},typenames_inds{iset_ptr}]=psg_select_util(sel,sas{iset});
                            disp(sprintf('dataset %2.0f (%30s): %2.0f stimuli, %2.0f stimuli selected',iset,labels{iset},nstims_each(iset),length(typenames_sel{iset_ptr})));
                            disp('typenames selected:');
                                disp(typenames_sel{iset_ptr});
                        end
                        ifok=getinp('1 if ok','d',[0 1],0);
                        if_center=getinp('1 to center the data','d',[0 1],1);
                        for iset_ptr=1:nproc_sets
                            iset=proc_sets(iset_ptr); 
                            nsets=nsets+1; %a new dataset for each file rotated into consensus
                            for id=1:dim_max
                                d_unrot=ds{iset}{id};
                                %ds{nsets}{id}=znew{id}(:,:,iset_ptr);
                            end
                        end
                    end %next iset_ptr to rotate
                    disp('not implemented yet')
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
        %propagate pipeline from previous processing steps
        sout.pipeline=pipelines{iset};
        opts_used=psg_write_coorddata([],ds{iset},sout,setfield(opts_write,'data_fullname_def',data_fullname_def));
    end
end

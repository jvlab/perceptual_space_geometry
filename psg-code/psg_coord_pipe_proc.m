%psg_pipe_coord_proc reads one or more coordinate datasets from experiments 
% processes, either by forming a consensus or simple transformations,
% and simple plots
%
% consensus forms a consensus of multiple datasets with same stimuli, and
%  also the Procrustes rotation of each set into the consensus
%
%The datasets that are combined must have the same stimuli in the same order, and the same number of dimensions.
%All datasets produced have the same stimuli in the same order, so no metadata file is created.
%
%For knitting together datasets that have only partially-overlapping
%stimuli, or stimuli  in a different order, see psg_align_knit_demo.
%
% Note that the pipeline saved for consensus is short, and does not have
%  the transformations to rotate into consensus.  For a full pipeline, use psg_align_knit_demo.
% procrustes finds the procrustes rotation of one dataset into a reference dataset,
%  and saves the transformation.
% pca_rotation allows for PC's to be determined based on a subset of stimuli
%  grouped by dimension -- e.g., the first two pc's maximize the variance
%  due to stimuli that contain 'b', and the next two pc's maximize the
%  remaining variance due to stimuli that contain 'g'
% nest causes the coords for a lower dimension to be a subset of the coords for a higher dimension
%
% The rotations applied by 'consensus' and 'pca_rotation' are not saved,
% but can be recoverered by doing a Procrustes rotatoin of the original
% dataset into the pca_rotation or consensus dataset, and then looking at
% pipeline.opts.transforms{dim} (which is saved)
%
% The output datafile has the fields
%   dim[a], where [a] is a character string 1,...,7,8,9,10,
%   stim_labels, as in the coord data files
%   pipeline, indicating the processing pipeline to create this file
%    pipeline.type: a string, e.g., 'rotate_to_consensus,'consensus','pca_rotation' 
%    pipeline.sets: the 'set' field of data file used to construct this dataset [may be absent if dataset is a consensus]
%    pipeline.sets_combined: cell array of 'set' fields of the datasets combined to create this dataset
%    pipeline.file_list: cell array of 'labels' fields of the datasets combined to create this dataset
%    pipeline.opts: options used
%  for type='rotate_to_consensus', sets is the file that is rotated to the consensus,
%    while sets_combined, file_list indicate the collection of files used to form the consensus
%  for type='pca_rotation' and 'nest', sets_combined' and file_list' are not present, since processing only uses a single dataset
%  for type='consensus', only sets_combined and file_list are present, since all of those files are used together
%
% 18Feb24: allow for propagtion of pipeline; show available files after each processing step
% 20Feb24: add pca rotation with option of restricting to specific stimuli based on stimulus names
% 21Feb24: fix pipeline so pipeline.sets.pipeline, pipeline.sets_combined.pipeline shows previous processing
% 21Feb24: tweak some dialog defaults, give warning if not enough stimuli selected for pca rotation
% 23Feb24: allow for model datasets to be read
% 26Feb24: modularize psg_coord_pipe_util
% 20Mar24: use psg_coords_fillin if datasets have missing dimensions
% 22Mar24: add Procrustes alignment
% 25Mar24: add simple transformations via psg_get_transform
% 26May24: allow for choice of initialization of consensus
% 26Apr25: add option to apply a geometric transformation from a file; change pipeline type to be string in pipe_types
% 12Jun25: add option to do pca around centroid but then restore centroid
% 12Jun25: add nest
% 
%  See also: PSG_GET_COORDSETS, PSG_QFORM2COORD_PROC, PSG_READ_COORDDATA, PSG_WRITE_COORDDATA, PSG_PLOTCOORDS,
%    PSG_COORD_PIPE_UTIL, SVD, PSG_COORDS_FILLIN, PSG_GEO_PROCRUSTES, PSG_GET_TRANSFORM, PSG_ALIGN_KNIT_DEMO,
%    PSG_GET_GEOTRANFORMS, PSG_GEOMODELS_APPLY.
%
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_pcon') opts_pcon=struct(); end %for procrustes_consensus
if ~exist('opts_procrustes') opts_procrustes=struct(); end %for simple procrustes
if ~exist('opts_write') opts_write=struct(); end %for psg_write_coorddata
if ~exist('opts_plot') opts_plot=struct(); end % for psg_plotcoords
if ~exist('opts_transform') opts_transform=struct(); end % for psg_get_transform
if ~exist('opts_geot') opts_geot=struct(); end; % for psg_get_geotransforms
%
pipe_types_long={'done (and proceed to write data files)',...
    'create a consensus dataset from datasets with matching entries',...
    'pca rotation of individual dataset(s), optionally centered',...
    'procrustes alignment of one or more sets to a given set',...
    'apply simple linear transformations and offset',...
    'plot coordinates of individual dataset(s) [no new datasets created]',...
    'apply a geometric transformation read from an auxiliary file',...
    'make coordinates strictly nested'};
pipe_types={'','consensus','pca_rotation','procrustes','linear_transformation','plot_coords','geometric_transformation','nest'}; %this string will appear in the metadata pipeline
pipe_minsets=[2 1 1 1 1 1 1];
npipe_types=length(pipe_types)-1;
%
disp('This will process one or more coordinate datasets and create new coordinate datasets, via simple transformations or consensus.');
%
nsets_signed=getinp('number of datasets (negative to use dialog box, data only)','d',[-100 100]);
[sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,[],nsets_signed);
nsets=length(sets); %number of sets actually read
for iset=1:nsets
    [sets{iset},ds{iset}]=psg_coords_fillin(sets{iset},ds{iset});
end
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
    pipelines{iset}=struct;
    disp(sprintf('dataset %2.0f (%30s): %2.0f stimuli, largest contained dimension is ',iset,labels{iset},nstims_each(iset)));
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
            disp(sprintf('dataset %2.0f (%30s): %2.0f stimuli',iset,labels{iset},nstims_each(iset)));
            %
            if iset>length(sets)
                sets{iset}.type='processed';
                sets{iset}.dim_list=[1:dim_max];
                sets{iset}.nstims=nstims;
                sets{iset}.label=labels{iset};
                sets{iset}.pipeline=pipelines{iset}; %which should be derived from earlier pipelines
           end
        end
        for ipipe=0:npipe_types
            disp(sprintf('%1.0f->%s',ipipe,pipe_types_long{ipipe+1}));
        end
        pipe_type=getinp('choice','d',[0 npipe_types]);
        if pipe_type>0
            proc_sets=getinp('datasets to process','d',[1 nsets]);
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
                pipe_type_selected=pipe_types{pipe_type+1};
                switch pipe_type_selected
                    case 'consensus'
                        pcon_init_method=getinp('method to use for initialization (>0: a specific set, 0 for PCA, -1 for PCA with forced centering, -2 for PCA with forced non-centering, 1->legacy','d',[-2 nsets],0);
                        if pcon_init_method>0
                            opts_pcon.initiailze_set=pcon_init_method;
                        else
                            if pcon_init_method==0
                                opts_pcon.initialize_set='pca';
                            elseif pcon_init_method==-1
                                opts_pcon.initialize_set='pca_center';
                            else
                                opts_pcon.initialize_set='pca_nocenter';
                            end
                        end
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
                        %create data and metadata for each original dataset rotated into consensus
                        for iset_ptr=1:nproc_sets
                            nsets=nsets+1; %a new dataset for each file rotated into consensus
                            for id=1:dim_max
                                ds{nsets}{id}=znew{id}(:,:,iset_ptr);
                            end
                            iset=proc_sets(iset_ptr);
                            sas{nsets}=sas{iset}; %take previous stimulus descriptors
                            labels{nsets}=sprintf('set %1.0f rotated into consensus of sets %s, allow_scale=%1.0f',iset,proc_sets_string,allow_scale);
                            nstims_each(nsets)=nstims;
                            pipelines{nsets}=psg_coord_pipe_util('rotate_to_consensus',setfield([],'opts_pcon_used',opts_pcon_used),sets{iset},file_list,sets_combined);
                        end
                        %create metadata for consensus
                        nsets=nsets+1; %a new dataset for each file rotated into consensus
                        nstims_each(nsets)=nstims;
                        ds{nsets}=consensus;
                        sas{nsets}=sas{proc_sets(1)}; %assume all stimulus descriptors are the same
                        labels{nsets}=sprintf('consensus of sets %s, allow_scale=%1.0f',proc_sets_string,allow_scale);
                        pipelines{nsets}=psg_coord_pipe_util(pipe_type_selected,setfield([],'opts_pcon_used',opts_pcon_used),[],file_list,sets_combined);
                    case 'linear_transformation'
                        transform=psg_get_transform(dim_max,opts_transform); %specify a transform
                        transforms=cell(1,dim_max);
                        for id=1:dim_max
                            transforms{id}.T=transform.T(1:id,1:id);
                            transforms{id}.b=transform.b;
                            transforms{id}.c=transform.c(1:id);
                        end
                        for iset_ptr=1:nproc_sets
                            nsets=nsets+1; %a new dataset for each transformed file
                            ds_new=cell(1,dim_max);
                            disp(sprintf(' applying specified transformation to set %1.0f (%s)',iset,labels{iset}));
                            for id=1:dim_max
                                ds_new{id}=transforms{id}.b*ds{iset}{id}*transforms{id}.T+repmat(transforms{id}.c,nstims,1); %truncate the transformation if needed
                            end
                            ds{nsets}=ds_new;
                            sas{nsets}=sas{iset}; %take previous stimulus descriptors
                            labels{nsets}=sprintf('set %1.0f linearly transformed',iset);
                            nstims_each(nsets)=nstims;
                            pipelines{nsets}=psg_coord_pipe_util(pipe_type_selected,...
                                setfields([],{'transforms'},{transforms}),... %save the transformations
                                sets{iset});
                        end
                    case 'geometric_transformation'
                        if exist('opts_geot_used')
                            if getinp('1 to use same transform file but different kind of transform','d',[0 1]);
                                opts_geot=setfield(opts_geot_used,'model_type',NaN);
                            end
                        end
                        [transforms_avail,dims_avail,geot_desc,opts_geot_used]=psg_get_geotransforms(opts_geot); %get a geometric transform from a file
                        if isempty(transforms_avail)
                            disp('no transform retrieved.')
                        else
                            disp(sprintf('transform retrieved, type %s, class %s, for dims %s',opts_geot_used.model_type,opts_geot_used.model_class,sprintf(' %2.0f',dims_avail)));
                            disp(sprintf('file: %s',opts_geot_used.filename))
                            disp(sprintf('description: %s',geot_desc))
                            transforms=cell(1,dim_max);
                            for id=1:dim_max
                                if ismember(id,dims_avail)
                                    transforms{id}=transforms_avail{id};
                                else
                                    disp(sprintf(' dimension %2.0f skipped, no transformation available',id));
                                end
                            end
                            for iset_ptr=1:nproc_sets
                                nsets=nsets+1; %a new dataset for each transformed file
                                ds_new=cell(1,dim_max);
                                disp(sprintf(' applying specified geometric transformation to set %1.0f (%s)',iset,labels{iset}));
                                for id=1:dim_max
                                    if ismember(id,dims_avail)
                                        ds_new{id}=psg_geomodels_apply(opts_geot_used.model_class,ds{iset}{id},transforms{id});
                                    else
                                        ds_new{id}=ds{iset}{id};
                                    end
                                end
                                ds{nsets}=ds_new;
                                sas{nsets}=sas{iset}; %take previous stimulus descriptors
                                labels{nsets}=sprintf('set %1.0f geometrically transformed (%s)',iset,geot_desc);
                                nstims_each(nsets)=nstims;
                                pipelines{nsets}=psg_coord_pipe_util(pipe_type_selected,...
                                    setfields([],{'transforms'},{transforms}),... %save the transformations
                                    sets{iset});
                            end
                        end
                    case 'procrustes'
                        ref_set=getinp('reference dataset','d',[1 nsets]);
                        ref_set_string=sprintf('%1.0f',ref_set);
                        allow_scale=getinp('1 to allow scaling','d',[0 1]);
                        ref_dim_choice=getinp('reference dimension (0 to match dimenson of adjusted set)','d',[0 dim_max],0);
                        if (ref_dim_choice==0)
                            ref_dim_string='same dim';
                        else
                            ref_dim_string=sprintf('dim %2.0f',ref_dim_choice);
                        end
                        file_list_ref=cell(1);
                        file_list_ref{1}=labels{ref_set};
                        sets_combined_ref=cell(1);
                        sets_combined_ref{1}=sets{ref_set};
                        %create data and metadata for each original dataset transformed by Procrustes
                        for iset_ptr=1:nproc_sets
                            iset=proc_sets(iset_ptr);
                            nsets=nsets+1; %a new dataset for each transformed file
                            ds_new=cell(1,dim_max);
                            disp(sprintf(' Procrustes-transforming set %1.0f (%s) to align with reference set %1.0f (%s)',...
                                iset,labels{iset},ref_set,labels{ref_set}));
                            d_fits=zeros(1,dim_max);
                            scales=zeros(1,dim_max);
                            transforms=cell(1,dim_max);
                            for id=1:dim_max
                                if (ref_dim_choice==0)
                                    id_ref=id;
                                else
                                    id_ref=ref_dim_choice;
                                end
                                [d_fit(id),ds_new{id},transforms{id},opts_procrustes_used]=...
                                    psg_geo_procrustes(ds{ref_set}{id_ref},ds{iset}{id},setfield(opts_procrustes,'if_scale',allow_scale)); %finds a procrustes model
                                scales(id)=transforms{id}.b;
                                disp(sprintf('for dim %2.0f of adj set aligned to dim %2.0f of ref set: scale factor %9.5f, d(procrustes fit) is %8.5f',...
                                    id,id_ref,scales(id),d_fit(id)));
                            end
                            ds{nsets}=ds_new;
                            sas{nsets}=sas{iset}; %take previous stimulus descriptors
                            labels{nsets}=sprintf('set %1.0f rotated into %s of %1.0f, allow_scale=%1.0f',iset,ref_dim_string,ref_set,allow_scale);
                            nstims_each(nsets)=nstims;
                            pipelines{nsets}=psg_coord_pipe_util(pipe_type_selected,...
                                setfields([],{'opts_procrustes_used','transforms'},{opts_procrustes_used,transforms}),... %save the transformations
                                sets{iset},file_list_ref,sets_combined_ref);
                        end
                    case 'pca_rotation' 
                        %
                        %get selection strings to select which stimuli are considered for successive rotations
                        % 
                        ifok=0;
                        while (ifok==0)
                            dim_groups=[];
                            tstring_dimspecs='';
                            ndgs=0;
                            typenames_sel=cell(1,nproc_sets);
                            typenames_inds=cell(1,nproc_sets);
                            tstring_dimspec=cell(0);
                            sel=cell(0);
                            ifok=1;
                            selcount_range=zeros(0,2);
                            while sum(dim_groups)<dim_max
                                ndgs=ndgs+1;
                                selcount_min=Inf;
                                selcount_max=-Inf;
                                selcount_needed=dim_max-sum(dim_groups([1:(ndgs-1)]));
                                dim_groups(ndgs)=getinp(sprintf('size of dimension group %1.0f, which starts at dimension %1.0f',ndgs,1+sum(dim_groups)),...
                                    'd',[1 dim_max-sum(dim_groups)],dim_max-sum(dim_groups));
                                sel{ndgs}=getinp('selection string or multiple strings, separated by |','s',[],'|');
                                tstring_dimspec{ndgs}=cat(2,'d[',sprintf('%1.0f ',[(1+sum(dim_groups(1:end-1))):sum(dim_groups)]),']:',sel{ndgs});
                                tstring_dimspecs=cat(2,tstring_dimspecs,tstring_dimspec{ndgs},' ');
                                for iset_ptr=1:nproc_sets
                                    iset=proc_sets(iset_ptr);                   
                                    [typenames_sel{iset_ptr}{ndgs},typenames_inds{iset_ptr}{ndgs}]=psg_select_util(sel{ndgs},sas{iset});
                                    disp(sprintf('typenames selected for dimension group %1.0f is %s',ndgs,tstring_dimspec{ndgs}));
                                    disp(typenames_sel{iset_ptr}{ndgs});
                                    disp(sprintf('dataset %2.0f (%30s): %2.0f stimuli, %2.0f stimuli selected',...
                                        iset,labels{iset},nstims_each(iset),length(typenames_sel{iset_ptr}{ndgs})));
                                    selcount_min=min(selcount_min,length(typenames_sel{iset_ptr}{ndgs}));
                                    selcount_max=max(selcount_max,length(typenames_sel{iset_ptr}{ndgs}));
                                end
                                selcount_range(ndgs,:)=[selcount_min,selcount_max];
                                if selcount_range(ndgs,1)<selcount_needed %do we have enough stimuli?
                                    disp(sprintf('not enough stimuli selected for this group: min available %2.0f, min needed to avoid singularity %2.0f',selcount_range(ndgs,1),selcount_needed));
                                    ifok=0;
                                end
                            end
                            tstring_dimspecs=deblank(tstring_dimspecs);
                            disp(sprintf('dimension groups: %s',tstring_dimspecs))
                            for idg=1:ndgs
                                disp(sprintf('dimension group %2.0f: string %10s, range of number of stimuli selected: [%3.0f %3.0f]',...
                                    idg,sel{idg},selcount_range(idg,:)));
                            end
                            ifok=getinp('1 if ok','d',[0 1],ifok);
                        end
                        %
                        if_center=getinp('1 to center the data, -1 to do pca around centroid and then restore','d',[-1 1],1);
                        %
                        %process each dataset
                        %
                        for iset_ptr=1:nproc_sets
                            iset=proc_sets(iset_ptr); 
                            nsets=nsets+1; %a new dataset for each file rotated
                            disp(sprintf('processing set %2.0f',iset));
                            for id=1:dim_max
                                d_unrot=ds{iset}{id};
                                if (if_center~=0)
                                    centroid=mean(d_unrot,1);
                                    d_unrot=d_unrot-repmat(centroid,nstims,1);
                                end
                                if if_center==-1
                                    centroid_restore=centroid;
                                else
                                    centroid_restore=zeros(1,id);
                                end
                                %for each dimension group (idg), determine directions of greatest variance for the
                                %stimuli selected by sel{idg}, that are orthogonal to previously-determined directions
                                %
                                d_rot=d_unrot;
                                for idg=1:ndgs
                                    d_avail=[sum(dim_groups(1:idg-1))+1:id]; %dimensions available to rotate
                                    if ~isempty(d_avail)
                                        d_fill=sum(dim_groups(1:idg-1))+1:min(sum(dim_groups(1:idg)),id); %the dimensions to be modified
                                        d_use=d_rot(typenames_inds{iset_ptr}{idg},d_avail); %data at selected coordinates and stimuli
                                        [u,s,v]=svd(d_use); %d_use=u*s*v', with u and v both orthogonal, so u*s=d_use*v
                                        %now use v to rotate the entire dataset, via the rotation defined by the selected coordinates
                                        d_rot(:,d_avail)=d_rot(:,d_avail)*v;
                                        %
                                        s_diag=diag(s(1:min(size(s)),1:min(size(s))));
                                        var_avail=sum(s_diag);
                                        d_fixed=min(dim_groups(idg),id-sum(dim_groups(1:idg-1))); %dimensions set by this group that won't be later changed
                                        var_remain=sum(s_diag(d_fixed+1:end));
                                        var_each=s_diag(1:d_fixed);
                                        var_string=sprintf('variance available %7.4f  remaining %7.4f (frac: %7.4f), on each dim:',var_avail,var_remain,var_remain/var_avail);
                                        var_string=cat(2,var_string,sprintf('%7.4f ',var_each));
                                        disp(sprintf('  for dim %1.0f: rotated coords %1.0f to %1.0f using selection %10s: %s',...
                                            id,min(d_avail),max(d_avail),sel{idg},var_string));
                                    end
                                end
                                ds{nsets}{id}=d_rot+centroid_restore; %coords for rotated dataset
                            end
                            %create metadata for rotated dataset:  the new set is derived from iset; there is no combination
                            sas{nsets}=sas{iset};
                            if if_center>=0
                                labels{nsets}=sprintf('set %1.0f, pca, centering=%1.0f, dimspecs %s',iset,if_center,tstring_dimspecs);
                            else
                                labels{nsets}=sprintf('set %1.0f, pca, centering=%1.0f (pca around centroid, centroid then restored), dimspecs %s',iset,if_center,tstring_dimspecs);
                            end
                            nstims_each(nsets)=nstims;
                            pca_rotation=struct;
                            pca_rotation.if_center=if_center;
                            pca_rotation.dim_groups=dim_groups;
                            pca_rotation.sel=sel;
                            pca_rotation.typenames_sel=typenames_sel{iset_ptr};
                            pca_rotation.typenames_inds=typenames_inds{iset_ptr};
                            pipelines{nsets}=psg_coord_pipe_util(pipe_type_selected,setfield([],'pca_rotation',pca_rotation),sets{iset});
                        end %next iset_ptr to rotate
                  case 'nest' 
                        %
                        %get dimension groups (as in pca_rotation)
                        % 
                        ifok=0;
                        while (ifok==0)
                            dim_groups=[];
                            tstring_dimspecs='';
                            ndgs=0;
                            tstring_dimspec=cell(0);
                            ifok=1;
                            while sum(dim_groups)<dim_max
                                ndgs=ndgs+1;
                                dim_groups(ndgs)=getinp(sprintf('size of dimension group %1.0f, which starts at dimension %1.0f',ndgs,1+sum(dim_groups)),...
                                    'd',[1 dim_max-sum(dim_groups)],dim_max-sum(dim_groups));
                                tstring_dimspec{ndgs}=cat(2,'d[',sprintf('%1.0f ',[(1+sum(dim_groups(1:end-1))):sum(dim_groups)]),']');
                                tstring_dimspecs=cat(2,tstring_dimspecs,tstring_dimspec{ndgs},' ');
                            end
                            tstring_dimspecs=deblank(tstring_dimspecs);
                            disp(sprintf('dimension groups: %s',tstring_dimspecs))
                            ifok=getinp('1 if ok','d',[0 1],ifok);
                        end
                        %
                        %process each dataset
                        %
                        for iset_ptr=1:nproc_sets
                            iset=proc_sets(iset_ptr); 
                            nsets=nsets+1; %a new dataset for each file nested
                            disp(sprintf('processing set %2.0f',iset));
                            ds{nsets}=cell(1,dim_max);
                            for id=1:dim_max
                                dgp=min(find(id<=cumsum(dim_groups)));
                                dim_source=sum(dim_groups(1:dgp));
                                ds{nsets}{id}=ds{iset}{dim_source}(:,1:id);
                                disp(sprintf(' dim %2.0f of new set %2.0f created from the first %2.0f coordinates of dim %2.0f of set %2.0f',id,nsets,id,dim_source,iset));
                            end
                            %create metadata for nested dataset:  the new set is derived from iset; there is no combination
                            sas{nsets}=sas{iset};
                            labels{nsets}=sprintf('set %1.0f, nested, dimspecs %s',iset,tstring_dimspecs);
                            nstims_each(nsets)=nstims;
                            nest=struct;
                            nest.dim_groups=dim_groups;
                            pipelines{nsets}=psg_coord_pipe_util(pipe_type_selected,setfield([],'nest',nest),sets{iset});
                        end %next iset_ptr                       
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

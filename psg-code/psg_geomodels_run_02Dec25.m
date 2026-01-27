%psg_geomodels_run: fitting geometric models, across a range of dimensions
%
% derived from psg_geomodels_test, but steps throgh a range of dimensions for reference and adjusted dataset,
% and puts results in a results structure
%
% differs from psg_procrustes_regr_test:
%    * includes a knitted built-in dataset
%    * procrustes with and without scaling
%    * regression model is called 'affine'
%    * for shuffling, reference dataset is permuted rather than adj dataset
%    * shuffling done for all nested models, including the zero model ("mean")
%    * psg_geomodels_define defines standard model defaults and nestings
%    * results from procrustes, affine, and projective models computed from a master routine
%    * d (goodness of fit) also computed with a denominator equal to variance of the original data
%       (otherwise, surrogates that inflate the variance will artifactually decrease the d)
%    * only one method is used for projective ('fmin')
%    * no fitting of projective to (unshuffled) residuals of Procrustes 
%    * frozen random numbers for each model and permutations, in case minimization is stochastic
%    * does nesting by dimension (for lower dimensions of adj dataset)
%
% NB: nesting by the dimension of the adjusted dataset only makes sense if
% the adjusted dataset is built up, one dimension at a time. This will
% always be the case if the adjusted dataset is created by MDS of a
% distance matrix (primary btc datasets), or by PCA of a response matrix
% (hlid datasets created directly from data or via affine merging with hlid_orn*merge)
% *but it will not be the case if the coordinates from each dimension are separately created by a consensus procedure*
%
% 28May24: allow for removal of specific models, test adequacy of input dimension
% 11Jun24: also use psg_geomodels_apply 
% 12Jun24: add if_pwaffine_details
% 12Jun24: begin to allow testing of lower dimensions as a nested model (if_nestbydim)
% 06Jun25: uses psg_align_coordsets via psg_commonstims to determine the overlapping stimuli for adj and ref datasets
% 02Dec25: frozen random numbers also applies to randperm
%
%   See also:  PROCRUSTES, PSG_GET_COORDSETS, PSG_PROCRUSTES_REGR_TEST,
%     PSG_PROCRUSTES_REGR_DEMO, PSG_GEO_GENERAL, PSG_GEOMODELS_DEFINE,
%     PSG_GEO_PROCRUSTES, PSG_GEO_AFFINE, PSG_GEO_PROJECTIVE, PERSP_APPLY, PSG_GEOMEODELS_TEST,
%     PSG_GEOMODELS_SUMM, PSG_GEOMODELS_APPLY, PSG_ALIGN_COORDSETS, PSG_COMMONSTIMS.
%
if ~exist('if_model_select') if_model_select=1;end
model_types_def=psg_geomodels_define(if_model_select);
model_types=model_types_def.model_types;
%
nmodels=length(model_types);
%
if ~exist('ref_dim_list') ref_dim_list=[2 3]; end
if ~exist('adj_dim_list') adj_dim_list=[2 3]; end
%
if ~exist('if_center') if_center=1; end
if ~exist('nshuff') nshuff=100; end
if ~exist('if_cycle') if_cycle=1; end %option for projective fit, but ignored as long as 'fmin' is used in persp_xform_find
if ~exist('opts_read') opts_read=struct; end
if ~exist('opts_rays') opts_rays=struct; end
if ~exist('opts_qpred') opts_qpred=struct; end
%
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
%
if_builtin=getinp('>0 to use built-in datasets (1-> bgca3pt MC and BL; 2-> bc55pooled MC) or 0 to specify','d',[0 2]);
switch if_builtin
    case 1
        if ~exist('ref_file') ref_file='./psg_data/bgca3pt_coords_MC_sess01_10.mat'; end
        if ~exist('adj_file') adj_file='./psg_data/bgca3pt_coords_BL_sess01_10.mat'; end
        ref_dim_max=7;
        adj_dim_max=7;
    case 2
        ref_file='./psg_data/bc55pooled_coords_MC.mat';
        adj_file='./psg_data/bc55pooled9.mat';
        nsets=2;
        opts_read=filldefault(opts_read,'if_log',1);
        opts_read.input_type=[1 2];
        opts_read.data_fullname_def=ref_file;
        opts_read.setup_fullname_def=adj_file;
        [sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred,nsets);
        ref_dim_max=max(sets{1}.dim_list);
        adj_dim_max=max(sets{2}.dim_list);
    otherwise       
        disp('dataset 1 will be reference, dataset 2 will be adjusted to fit.');
        nsets=2;
        opts_read=filldefault(opts_read,'if_log',1);
        [sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred,nsets);
        ref_file=sets{1}.label_long;
        adj_file=sets{2}.label_long;
        ref_dim_max=max(sets{1}.dim_list);
        adj_dim_max=max(sets{2}.dim_list);
end
%
if if_builtin~=1
    [ds,sas]=psg_commonstims(ds,sas); %only look at stimuli shared by adj and ref datasets
end
%model-fitting options
if_center=getinp('1 to center the data','d',[0 1],if_center);
%
disp(sprintf('reference dataset: %s',ref_file));
disp(sprintf('dataset to be adjusted: %s',adj_file));
%
ref_dim_list=getinp('list of reference dataset dimensions to use','d',[1 ref_dim_max],ref_dim_list);
adj_dim_list=getinp('list of adjusted  dataset dimensions to use','d',[1 adj_dim_max],adj_dim_list);
if_nestbydim=getinp('1 to also do nesting by dimension, within model type','d',[0 1]); 
if_geomodel_check=getinp('1 to also recheck computation of transform via psg_geomodels_apply','d',[0 1],0);
if_pwaffine_details=getinp('1 to show details in summary for piecewise-affine minimizations','d',[0 1],0);
%
ref_dim_list=sort(ref_dim_list);
adj_dim_list=sort(adj_dim_list);
results=cell(max(ref_dim_list),max(adj_dim_list));
%
for iref_ptr=1:length(ref_dim_list)
    for iadj_ptr=1:length(adj_dim_list)
        ref_dim=ref_dim_list(iref_ptr);
        adj_dim=adj_dim_list(iadj_ptr);
        disp(' ');
        disp(sprintf('analyzing with ref dim %3.0f and adj dim %3.0f',ref_dim,adj_dim));
        disp(' ');
        %
        r=struct;
        r.model_types_def=model_types_def;
        r.ref_dim=ref_dim;
        r.adj_dim=adj_dim;
        r.ref_file=ref_file;
        r.adj_file=adj_file;
        r.d_shuff_dims='d1: model, d2: shuffle, d3: nested model, d4: normalization type';
        r.surrogate_count_dims='d1: model, d2: nested model, d3: normalization type';
        r.nshuff=nshuff;
        r.if_center=if_center;
        r.if_cycle=if_cycle;
        switch if_builtin
            case 1
                ref=getfield(load(ref_file),sprintf('dim%1.0f',ref_dim));
                adj=getfield(load(adj_file),sprintf('dim%1.0f',adj_dim));
            case 2
                ref=ds{1}{ref_dim};
                adj=ds{2}{adj_dim};
            otherwise
                ref=ds{1}{ref_dim};
                adj=ds{2}{adj_dim};
        end
        npts=size(ref,1);
        disp(sprintf('number of data points: %3.0f',npts));
        %
        tstring=sprintf(' ref %s dim %1.0f adj %s dim %1.0f center %1.0f',...
            ref_file,ref_dim,adj_file,adj_dim,if_center);
        %
        npts=size(ref,1);
        if if_center
            ref=ref-repmat(mean(ref,1),npts,1);
            adj=adj-repmat(mean(adj,1),npts,1);
        end
        %
        d=zeros(nmodels,1);
        d_nonorth=zeros(nmodels,1);
        d_calc_types={'den: surrogate','den:  original'};
        %    
        transforms=cell(nmodels,1);
        transforms_nonorth=cell(nmodels,1);
        adj_model=cell(nmodels,1);
        adj_model_check=cell(nmodels,1);
        adj_model_nonorth=cell(nmodels,1);
        resids=cell(nmodels,1);
        d_check=zeros(nmodels,1);
        %
        opts_model_used=cell(nmodels,1);
        opts_model_nonorth_used=cell(nmodels,1);
        %
        %quantities for analysis of one model nested in another
        d_shuff=zeros(nmodels,nshuff,nmodels-1,length(d_calc_types)); %d1: model, d2: shuffle, d3: nested model, d4: normalization type
        opts_model_shuff_used=cell(nmodels,nshuff,nmodels-1);
        model_lastnested=zeros(nmodels,1);
        surrogate_count=zeros(nmodels,nmodels-1,length(d_calc_types));
        %
        %quantities for analysis of a model with a loewr adj dim nested in a higher adj dim
        d_shuff_nestdim=zeros(nmodels,nshuff,iadj_ptr-1,length(d_calc_types)); %d1: model, d2: shuffle, d3: nested dim, 4: normalization type
        opts_model_shuff_used_nestdim=cell(nmodels,nshuff,iadj_ptr-1); 
        surrogate_count_nestdim=zeros(nmodels,iadj_ptr-1,length(d_calc_types));
        nestdim_list=adj_dim_list(1:iadj_ptr-1);
        %
        d_den=sum(sum((ref-repmat(mean(ref,1),npts,1)).^2,1));
        if ref_dim<adj_dim
            ref_aug=[ref,zeros(npts,adj_dim-ref_dim)];
        else
            ref_aug=ref;
        end
        %
        perms=zeros(nshuff,npts);
        if (if_frozen~=0)
            rng('default');
            if (if_frozen<0)
                rand(1,abs(if_frozen));
            end
        else
            rng('shuffle');
        end
        for ishuff=1:nshuff
            perms(ishuff,:)=randperm(npts);
        end
        for imodel=1:nmodels
            %
            if (if_frozen~=0)
                rng('default');
                if (if_frozen<0)
                    rand(1,abs(if_frozen));
                end
            else
                rng('shuffle');
            end
            model_type=model_types{imodel};
            disp(' ');
            disp(sprintf('model type: %s, adj dim %2.0f ref dim %2.0f',model_type,adj_dim,ref_dim))
            if adj_dim<model_types_def.(model_type).min_inputdims
                disp(sprintf(' model type skipped, requires input dimension of at least %2.0f',model_types_def.(model_type).min_inputdims));
            else
                opts_model=model_types_def.(model_type).opts;
                model_class=model_types_def.(model_type).class;
                [d(imodel),adj_model{imodel},transforms{imodel},opts_model_used{imodel}]=psg_geo_general(ref,adj,model_class,opts_model);
                %
                %verify the model   
                %
                switch model_class
                    case {'mean','procrustes','affine'}
                        adj_model_check{imodel}=transforms{imodel}.b*adj*transforms{imodel}.T+repmat(transforms{imodel}.c,npts,1);
                    case 'projective'
                        %note change of variable names between psg_geo_projective and persp_apply
                        adj_model_check{imodel}=persp_apply(transforms{imodel}.T,transforms{imodel}.c,transforms{imodel}.p,adj);
                    case 'pwaffine'
                        adj_model_check{imodel}=psg_pwaffine_apply(transforms{imodel},adj);
                        %also compare orthogonal and non-orthogonal versions
                        [d_nonorth(imodel),adj_model_nonorth{imodel},transforms_nonorth{imodel},opts_model_nonorth_used{imodel}]=psg_geo_general(ref,adj,model_class,setfield(opts_model,'if_orth',0));
                end
                if if_geomodel_check
                    adj_model_geo=psg_geomodels_apply(model_class,adj,transforms{imodel});
                    maxdev=max(abs(adj_model_check{imodel}(:)-adj_model_geo(:)));
                    if maxdev>0
                        warning(sprintf('geomodel check fails: maxdev=%16.14f',maxdev));
                    end
                end
                %verify model
                model_check=max(abs(adj_model{imodel}(:)-adj_model_check{imodel}(:)));
                disp(sprintf('model check: %12.5f',model_check));
                %calculate residuals and verify d
                resids{imodel}=ref_aug-adj_model{imodel};
                d_check(imodel)=sum(resids{imodel}(:).^2)/d_den;
                disp(sprintf('d: %12.7f   d_check: %12.7f   diff: %12.7f',...
                    d(imodel),d_check(imodel),d(imodel)-d_check(imodel)));
                %check the reconstitution
                disp('transform parameters:')
                if isfield(transforms{imodel},'b') disp('b'); disp(transforms{imodel}.b); end
                %
                if isfield(transforms{imodel},'T') disp('T'); disp(transforms{imodel}.T); end
                if isfield(transforms{imodel},'c') disp('c'); disp(transforms{imodel}.c); end
                %for projective
                if isfield(transforms{imodel},'p') disp('p transpose'); disp(transforms{imodel}.p'); end
                %for piecewise
                if isfield(transforms{imodel},'vcut') disp('vcut'); disp(transforms{imodel}.vcut); end
                if isfield(transforms{imodel},'acut') disp('acut'); disp(transforms{imodel}.acut); end
                %
                %shuffles
                %
                if nshuff>0 & isempty(opts_model_used{imodel}.warnings)
                    nested_types=model_types_def.(model_type).nested;
                    for inest=1:length(nested_types)
                        nested_type=nested_types{inest};
                        nest_ptr=strmatch(nested_type,model_types,'exact');
                        disp(sprintf(' doing shuffles with residuals from nested model %2.0f (%s)',nest_ptr,nested_type));
                        for ishuff=1:nshuff
                            %if perms is trival, shuffled=ref_aug; otherwise the residuals from the nested model are shuffled
                            %note that the denominator used to normalize d has changed.
                            shuffled=adj_model{nest_ptr}+(ref_aug(perms(ishuff,:),:)-adj_model{nest_ptr}(perms(ishuff,:),:));
                            %d_shuff_orig is calculated with surrogate denominator; we also want to calculate it with original denominator
                            %also, turn off fmin display for shuffles
                            [d_shuff_orig,adj_model_shuff,transform_shuffle,opts_model_shuff_used{imodel,ishuff,inest}]=...
                                psg_geo_general(shuffled,adj,model_class,setfield(opts_model,'if_display',0));
                            resids_shuff=shuffled-adj_model_shuff; %deviation of model from shuffle surrogate
                            d_shuff(imodel,ishuff,inest,1)=d_shuff_orig;
                            d_shuff(imodel,ishuff,inest,2)=sum(resids_shuff(:).^2)/d_den;
                        end
                        for id_calc_type=1:length(d_calc_types)
                            surrogate_count(imodel,nest_ptr,id_calc_type)=sum(double(d(imodel)>=d_shuff(imodel,:,inest,id_calc_type)));
                            disp(sprintf('d is >= d_shuff for %5.0f of %5.0f shuffles (%s), d_shuffles: range [%8.5f %8.5f], mean %8.5f s.d. %8.5f',...
                                surrogate_count(imodel,nest_ptr,id_calc_type),nshuff,d_calc_types{id_calc_type},...
                                min(d_shuff(imodel,:,inest,id_calc_type)),max(d_shuff(imodel,:,inest,id_calc_type)),...
                                mean(d_shuff(imodel,:,inest,id_calc_type)),std(d_shuff(imodel,:,inest,id_calc_type))));
                        end
                        model_lastnested(imodel)=nest_ptr;
                    end %inest
                    %are lower values of adj_dim available to test?
                    if if_nestbydim==1
                        for iadj_ptr_nest=1:iadj_ptr-1
                            adj_dim_nest=adj_dim_list(iadj_ptr_nest);
                            if adj_dim_nest<model_types_def.(model_type).min_inputdims
                                disp(sprintf('   skipping nesting model type %s with adj dim %2.0f in adj dim %2.0f, with ref dim %2.0f',...
                                    model_type,adj_dim_nest,adj_dim,ref_dim))
                            else
                                disp(sprintf(' evaluating nesting model type %s with adj dim %2.0f in adj dim %2.0f, with ref dim %2.0f',...
                                    model_type,adj_dim_nest,adj_dim,ref_dim))
                                    %recover lower-dim data
                                switch if_builtin
                                    case 1
                                        adj_nest=getfield(load(adj_file),sprintf('dim%1.0f',adj_dim_nest));
                                    case 2
                                        adj_nest=ds{2}{adj_dim_nest};
                                    otherwise
                                        adj_nest=ds{2}{adj_dim_nest};
                                end
                                if if_center
                                    adj_nest=adj_nest-repmat(mean(adj_nest,1),npts,1);
                                end
                                transform_nest=results{ref_dim,adj_dim_nest}.transforms{imodel};
                                adj_model_nest=psg_geomodels_apply(model_class,adj_nest,transform_nest);
                                if size(adj_model_nest,2)<size(ref_aug,2)
                                    %this can happen if size(adj_model_nest,2)<=size(ref,2) but size(ref,2)<size(adj_model,2)
                                    %so ref was augmented to match size(adj_model) but adj_model_nest
                                    %was only fit to a target of dimension size(ref,2)
                                    adj_model_nest=[adj_model_nest,zeros(npts,size(ref_aug,2)-size(adj_model_nest,2))];
                                end
                                %now do shuffles, find model params, and tabulate d
                                for ishuff=1:nshuff
                                   %if perms is trival, shuffled=ref_aug; otherwise the residuals from the nested model are shuffled
                                   %note that the denominator used to normalize d has changed.
                                   shuffled=adj_model_nest+(ref_aug(perms(ishuff,:),:)-adj_model_nest(perms(ishuff,:),:));
                                   %d_shuff_orig is calculated with surrogate denominator; we also want to calculate it with original denominator
                                   %also, turn off fmin display for shuffles
                                   [d_shuff_orig,adj_model_shuff,transform_shuffle,opts_model_shuff_used_nestdim{imodel,ishuff,iadj_ptr_nest}]=...
                                      psg_geo_general(shuffled,adj,model_class,setfield(opts_model,'if_display',0));
                                   resids_shuff=shuffled-adj_model_shuff; %deviation of model from shuffle surrogate
                                   d_shuff_nestdim(imodel,ishuff,iadj_ptr_nest,1)=d_shuff_orig;
                                   d_shuff_nestdim(imodel,ishuff,iadj_ptr_nest,2)=sum(resids_shuff(:).^2)/d_den;
                                end %ishuff
                            for id_calc_type=1:length(d_calc_types)
                                surrogate_count_nestdim(imodel,iadj_ptr_nest,id_calc_type)=sum(double(d(imodel)>=d_shuff_nestdim(imodel,:,iadj_ptr_nest,id_calc_type)));
                                disp(sprintf('d is >= d_shuff for %5.0f of %5.0f shuffles (%s), d_shuffles: range [%8.5f %8.5f], mean %8.5f s.d. %8.5f',...
                                    surrogate_count_nestdim(imodel,iadj_ptr_nest,id_calc_type),nshuff,d_calc_types{id_calc_type},...
                                    min(d_shuff_nestdim(imodel,:,iadj_ptr_nest,id_calc_type)),max(d_shuff_nestdim(imodel,:,iadj_ptr_nest,id_calc_type)),...
                                    mean(d_shuff_nestdim(imodel,:,iadj_ptr_nest,id_calc_type)),std(d_shuff_nestdim(imodel,:,iadj_ptr_nest,id_calc_type))));
                              end
                            end %enough dimensions for model?
                        end %iadj_ptr_nest
                    end %if_nestbydim
                end %if shuffled
            end %model input dimension test
        end
        %final summary
        disp(' ');
        disp('summary')
        disp(sprintf('reference dataset: %s',ref_file));
        disp(sprintf('dataset to be adjusted: %s',adj_file));
        disp(sprintf('reference dataset dimension: %2.0f',ref_dim));
        disp(sprintf('adjusted  dataset dimension: %2.0f',adj_dim));
        for imodel=1:nmodels
            model_type=model_types{imodel};
            model_class=model_types_def.(model_type).class;
            disp(sprintf('model %20s: class: %s d: %8.5f',model_type,model_class,d(imodel)));
            if ~isempty(opts_model_used{imodel}) %in case model did not pass dimension test
                if (nshuff>0)
                    nested_types=model_types_def.(model_type).nested;
                    for inest=1:length(nested_types)
                        nested_type=nested_types{inest};
                        nest_ptr=strmatch(nested_type,model_types,'exact');
                        for id_calc_type=1:length(d_calc_types)
                            disp(sprintf('nested %20s: d is >= d_shuff for %5.0f of %5.0f shuffles (%s), d_shuffles: range [%8.5f %8.5f], mean %8.5f s.d. %8.5f',...
                                nested_type,surrogate_count(imodel,nest_ptr,id_calc_type),nshuff,d_calc_types{id_calc_type},...
                                min(d_shuff(imodel,:,inest,id_calc_type)),max(d_shuff(imodel,:,inest,id_calc_type)),...
                                mean(d_shuff(imodel,:,inest,id_calc_type)),std(d_shuff(imodel,:,inest,id_calc_type))));
                        end
                    end %inest
                    if if_nestbydim==1
                        for iadj_ptr_nest=1:iadj_ptr-1
                            adj_dim_nest=adj_dim_list(iadj_ptr_nest);
                            if adj_dim_nest>=model_types_def.(model_type).min_inputdims
                                for id_calc_type=1:length(d_calc_types)
                                    disp(sprintf('%27s: d is >= d_shuff for %5.0f of %5.0f shuffles (%s), d_shuffles: range [%8.5f %8.5f], mean %8.5f s.d. %8.5f',...
                                        sprintf('dn %s dim %2.0f',model_type,adj_dim_nest),surrogate_count_nestdim(imodel,iadj_ptr_nest,id_calc_type),nshuff,d_calc_types{id_calc_type},...
                                        min(d_shuff_nestdim(imodel,:,iadj_ptr_nest,id_calc_type)),max(d_shuff_nestdim(imodel,:,iadj_ptr_nest,id_calc_type)),...
                                        mean(d_shuff_nestdim(imodel,:,iadj_ptr_nest,id_calc_type)),std(d_shuff_nestdim(imodel,:,iadj_ptr_nest,id_calc_type))));
                                end %id_calc_type
                             end %dims OK
                        end
                    end %if_nestbydim
                end
                if strcmp(model_class,'pwaffine') & (if_pwaffine_details==1)
                    disp(sprintf('standard minimization for %s:',model_type))
                    disp(opts_model_used{imodel}.fmin);
                    disp(opts_model_used{imodel}.fmin.output);
                    disp(transforms{imodel});
                    %
                    disp(sprintf('minimization for %s with if_orth=0:',model_type))
                    disp(opts_model_nonorth_used{imodel}.fmin);
                    disp(opts_model_nonorth_used{imodel}.fmin.output);
                    disp(transforms_nonorth{imodel});
                end
            end
        end
        %
        r.d=d;
        %
        r.transforms=transforms;
        r.opts_model_used=opts_model_used;
        r.d_shuff=d_shuff;
        r.surrogate_count=surrogate_count;
        %
        if (if_nestbydim)
            r.nestdim_list=nestdim_list;
            r.opts_model_shuff_used_nestdim=opts_model_shuff_used_nestdim;
            r.d_shuff_nestdim=d_shuff_nestdim;
            r.surrogate_count_nestdim=surrogate_count_nestdim;
        end
        %
        results{ref_dim,adj_dim}=r;
    end %adj_dim_ptr
end %ref_dim_ptr
disp('can save ''results'' structure with result summary');


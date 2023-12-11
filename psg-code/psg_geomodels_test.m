%psg_geomodels_test: test: test fitting geometric models
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
%    * no fitting of projective to residuals of Procrustes
%    * frozen random numbers for each model and permutations, in case minimization is stochastic
%
% plan is to wrap procrustes in psg_geomodel_procrustes, 
% separate out psg_geo_procrustes (and clean up c), psg_geo_affine, psg_geo_projective
% use this as a platform to add piecewise models
%
% 05Dec23: begin to add piecewise affine
% 07Dec23: piecewise affine added, summary added
% 11Dec23: add test of non-orthogonal option in psg_geo_pwaffine
%
%   See also:  PROCRUSTES, PSG_GET_COORDSETS, PSG_PROCRUSTES_REGR_TEST,
%     PSG_PROCRUSTES_REGR_DEMO, PSG_GEO_GENERAL, PSG_GEOMODELS_DEFINE,
%     PSG_GEO_PROCRUSTES, PSG_GEO_AFFINE, PSG_GEO_PROJECTIVE, PERSP_APPLY.
%
model_types_def=psg_geomodels_define();
model_types=model_types_def.model_types;
%
nmodels=length(model_types);
%
if ~exist('ref_dim') ref_dim=3; end
if ~exist('adj_dim') adj_dim=3; end
%
if ~exist('if_center') if_center=1; end
if ~exist('nshuff') nshuff=100; end
if ~exist('nshuff_regr') nshuff_regr=100; end
if ~exist('if_cycle') if_cycle=1; end %option for projective fit, but ignored as long as 'fmin' is used in persp_xform_find
if ~exist('opts_read') opts_read=struct; end
if ~exist('opts_rays') opts_rays=struct; end
if ~exist('opts_qpred') opts_qpred=struct; end
if ~exist('opts_persp') opts_persp=struct; end
%
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
%
if_builtin=getinp('use built-in datasets (1: bgca3pt MC and BL; 2: bc55pooled MC','d',[0 2]);
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
disp(sprintf('reference dataset: %s',ref_file));
disp(sprintf('dataset to be adjusted: %s',adj_file));
ref_dim=getinp('reference dataset dimension to use','d',[1 ref_dim_max],ref_dim);
adj_dim=getinp('adjusted  dataset dimension to use','d',[1 adj_dim_max],adj_dim);
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
%model-fitting options
if_center=getinp('1 to center the data','d',[0 1],if_center);
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
d_calc_types={'den: surrogate','den:  original'};
%    
d_shuff=zeros(nmodels,nshuff,nmodels-1,length(d_calc_types)); %d1: model, d2: shuffle, d3: nested model, d4: normalization type
transforms=cell(nmodels,1);
adj_model=cell(nmodels,1);
adj_model_check=cell(nmodels,1);
resids=cell(nmodels,1);
d_check=zeros(nmodels,1);
opts_model_used=cell(nmodels,1);
opts_model_shuff_used=cell(nmodels,nshuff,nmodels-1);
model_lastnested=zeros(nmodels,1);
surrogate_count=zeros(nmodels,nmodels-1,length(d_calc_types));
%
d_den=sum(sum((ref-repmat(mean(ref,1),npts,1)).^2,1));
if ref_dim<adj_dim
    ref_aug=[ref,zeros(npts,adj_dim-ref_dim)];
else
    ref_aug=ref;
end
if (if_frozen~=0)
    rng('default');
    if (if_frozen<0)
        rand(1,abs(if_frozen));
    end
else
    rng('shuffle');
end
perms=zeros(nshuff,npts);
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
    disp(sprintf('model type: %s',model_type))
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
            [d_nonorth,adj_model_nonorth,transforms_nonorth,opts_model_nonorth_used]=psg_geo_general(ref,adj,model_class,setfield(opts_model,'if_orth',0));
    end
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
                %check this next line -- possibly an issue that the denominator of d has changed?
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
    end %if shuffled
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
    disp(sprintf('model %20s:  d: %8.5f',model_type,d(imodel)));
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
    end
    if strcmp(model_type,'pwaffine')
        disp(sprintf('standard minimization for piecewise affine:'))
        disp(opts_model_used{imodel}.fmin);
        disp(opts_model_used{imodel}.fmin.output);
        disp(transforms{imodel});
        %
        disp(sprintf('minimization for piecewise affine with if_orth=0:'))
        disp(opts_model_nonorth_used.fmin);
        disp(opts_model_nonorth_used.fmin.output);
        disp(transforms_nonorth);
    end
end

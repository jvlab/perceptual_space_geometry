function [results,opts_geofit_used]=psg_geomodels_fit(d_ref,d_adj,opts_geofit)
% [results,opts_geofit_used]=psg_geomodels_fit(d_ref,d_adj,opts_geofit) fits geometric models
% (scale, affine, projective, piecewise affine) across a range of dimensions, and do statistis
%
% This is the modularized version of psg_geomodels_run, in turn derived from psg_geomodels_test, 
% steps throgh a range of dimensions for reference and adjusted dataset
%
% d_adj: dataset to be adjusted, k-dimensional model is in d_adj{k}, one row per stimulus
% d_ref: reference dataset; k-dimensional model is in d_ref{k}
% opts_geofit: a structure with these fields:
%   model_types_def:  model type definitions, typically returned by psg_geomodels_define
%      specifies which models are to be fit; if empty, requested interactively
%   ref_dim_list: reference dataset dimensions to analyze, defaults to [2 3]
%   adj_dim_list: adjusted dataset dimensions to analyze, defaults to ref_dim_list
%   if_center: 1 (default) to center the data
%   if_frozen: 1 (default) to use fronzen random numbers, 0 for random each time, <0 to specify a seed)
%   if_log: 1 (default) to log
%   if_summary: 1 (default) to show a summary
%   nshuffs: number of shuffles, defaults to 100
%   if_nestbydim: 1 (default: 0) to also do statistics for neting by dimension within adjusted dataset
%      nesting by the dimension of the adjusted dataset only makes sense if
%      the adjusted dataset is built up, one dimension at a time. This will
%      always be the case if the adjusted dataset is created by MDS of a
%      distance matrix (primary btc datasets), or by PCA of a response matrix
%      *but it will not be the case if the coordinates from each dimension are 
%      created by a consensus procedure, or rotated
%  [not typically needed]
%   if_geomodel_check: 1 (default: 0) to recheck computation of transform via psg_geomodels_apply
%   if_pwaffine_details: 1 (default: 0) to show details in summary for piecewise-affine minimizations)
%   persp_method: 'fmin' (default) or 'oneshot' (Zhang method, see persp_xform_find)
%   persp_if_cycle: 1 (default), or 0, variants of Zhang method
%
%  results: cell(max(ref_dim_list),max(adj_dim_list)), a structure with results
%  opts_geofit_used: options used, and warnings field
%
%   See also:  PSG_GEOMODELS_RUN, PSG_GEO_GENERAL, PSG_GEOMODELS_DEFINE.
%
if (nargin<=2)
    opts_geofit=struct;
end
opts_geofit=filldefault(opts_geofit,'model_types_def',[]);
opts_geofit=filldefault(opts_geofit,'ref_dim_list',[2 3]);
opts_geofit=filldefault(opts_geofit,'adj_dim_list',opts_geofit.ref_dim_list);
opts_geofit=filldefault(opts_geofit,'if_center',1);
opts_geofit=filldefault(opts_geofit,'if_frozen',1);
opts_geofit=filldefault(opts_geofit,'if_log',1);
opts_geofit=filldefault(opts_geofit,'if_summary',1);
opts_geofit=filldefault(opts_geofit,'nshuffs',100);
opts_geofit=filldefault(opts_geofit,'if_nestbydim',0);
opts_geofit=filldefault(opts_geofit,'persp_method','fmin');
opts_geofit=filldefault(opts_geofit,'persp_if_cycle',1);
opts_geofit=filldefault(opts_geofit,'if_geomodel_check',0);
opts_geofit=filldefault(opts_geofit,'if_geomodel_check_tol',10^-6); %tolerance for if_geomodel_check
opts_geofit=filldefault(opts_geofit,'if_pwaffine_details',0);
%
if ~isstruct(opts_geofit.model_types_def)
    opts_geofit.model_types_def=psg_geomodels_define(1);
end
%for compatibility with psg_geomodels_run
ref_dim_list=unique(opts_geofit.ref_dim_list);
adj_dim_list=unique(opts_geofit.adj_dim_list);
if_center=opts_geofit.if_center;
if_frozen=opts_geofit.if_frozen;
nshuff=opts_geofit.nshuffs;
model_types_def=opts_geofit.model_types_def;
if_nestbydim=opts_geofit.if_nestbydim;
if_geomodel_check=opts_geofit.if_geomodel_check;
if_pwaffine_details=opts_geofit.if_pwaffine_details;
%
opts_geofit_used=opts_geofit;
opts_geofit_used.warnings=[];
%
model_types=model_types_def.model_types;
nmodels=length(model_types);
%
results=cell(max(ref_dim_list),max(adj_dim_list));
%
for iref_ptr=1:length(ref_dim_list)
    for iadj_ptr=1:length(adj_dim_list)
        ref_dim=ref_dim_list(iref_ptr);
        adj_dim=adj_dim_list(iadj_ptr);
        %
        ref=d_ref{ref_dim};
        adj=d_adj{adj_dim};
        if (size(ref,2)~=ref_dim | size(adj,2)~=adj_dim) %make sure dimensions are correct
            wmsg=sprintf('ref dim (%2.0f) or adj dim (%2.0f) not found',ref_dim,adj_dim);
            opts_geofit_used.warnings=strvcat(opts_geofit_used.warnings,wmsg);
            if opts_geofit.if_log
                disp(wmsg);
            end
        else
            r=struct;
            r.model_types_def=model_types_def;
            r.ref_dim=ref_dim;
            r.adj_dim=adj_dim;
            r.d_shuff_dims='d1: model, d2: shuffle, d3: nested model, d4: normalization type';
            r.surrogate_count_dims='d1: model, d2: nested model, d3: normalization type';
            r.opts_geofit=opts_geofit;  
            npts=size(ref,1);
            if opts_geofit.if_log
               disp(' ');
               disp(sprintf('analyzing ref dim %3.0f and adj dim %3.0f, number of data points: %3.0f',ref_dim,adj_dim,npts));
            end
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
                tstring=sprintf('model type: %s, adj dim %2.0f ref dim %2.0f',model_type,adj_dim,ref_dim);
                if opts_geofit.if_log
                    disp(tstring);
                end
                if adj_dim<model_types_def.(model_type).min_inputdims
                    disp(sprintf(' model type skipped, requires input dimension of at least %2.0f',model_types_def.(model_type).min_inputdims));
                else
                    opts_model=model_types_def.(model_type).opts;
                    model_class=model_types_def.(model_type).class;
                    opts_model.if_cycle=opts_geofit.persp_if_cycle;
                    opts_model.method=opts_geofit.persp_method;
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
                        if abs(maxdev)>opts_geofit.if_geomodel_check_tol
                            wmsg=sprintf('%s: geomodel check fails: maxdev=%16.14f, tol: %16.14f',tstring,maxdev,opts_geofit.if_geomodel_check_tol);
                            opts_geofit_used.warnings=strvcat(opts_geofit_used.warnings,wmsg);
                            if opts_geofit.if_log
                                disp(wmsg);
                            end
                        end
                    end
                    %verify model
                    model_check=max(abs(adj_model{imodel}(:)-adj_model_check{imodel}(:)));
                    %calculate residuals and verify d
                    resids{imodel}=ref_aug-adj_model{imodel};
                    d_check(imodel)=sum(resids{imodel}(:).^2)/d_den;
                    if opts_geofit.if_log
                        disp(sprintf('model check: %12.5f',model_check));
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
                    end
                    %
                    %shuffles
                    %
                    if nshuff>0 & isempty(opts_model_used{imodel}.warnings)
                        nested_types=model_types_def.(model_type).nested;
                        for inest=1:length(nested_types)
                            nested_type=nested_types{inest};
                            nest_ptr=strmatch(nested_type,model_types,'exact');
                            if opts_geofit.if_log
                                disp(sprintf(' doing shuffles with residuals from nested model %2.0f (%s)',nest_ptr,nested_type));
                            end
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
                                if opts_geofit.if_log
                                    disp(sprintf('d is >= d_shuff for %5.0f of %5.0f shuffles (%s), d_shuffles: range [%8.5f %8.5f], mean %8.5f s.d. %8.5f',...
                                        surrogate_count(imodel,nest_ptr,id_calc_type),nshuff,d_calc_types{id_calc_type},...
                                        min(d_shuff(imodel,:,inest,id_calc_type)),max(d_shuff(imodel,:,inest,id_calc_type)),...
                                        mean(d_shuff(imodel,:,inest,id_calc_type)),std(d_shuff(imodel,:,inest,id_calc_type))));
                                end
                            end
                            model_lastnested(imodel)=nest_ptr;
                        end %inest
                        %are lower values of adj_dim available to test?
                        if if_nestbydim==1
                            for iadj_ptr_nest=1:iadj_ptr-1
                                adj_dim_nest=adj_dim_list(iadj_ptr_nest);
                                if adj_dim_nest<model_types_def.(model_type).min_inputdims
                                    if opts_geofit.if_log
                                        disp(sprintf('   skipping nesting model type %s with adj dim %2.0f in adj dim %2.0f, with ref dim %2.0f',...
                                            model_type,adj_dim_nest,adj_dim,ref_dim))
                                    end
                                else
                                    if opts_geofit.if_log
                                        disp(sprintf(' evaluating nesting model type %s with adj dim %2.0f in adj dim %2.0f, with ref dim %2.0f',...
                                            model_type,adj_dim_nest,adj_dim,ref_dim))
                                    end
                                    %recover lower-dim data
                                    adj_nest=d_adj{adj_dim_nest};
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
                                    if opts_geofit.if_log
                                        disp(sprintf('d is >= d_shuff for %5.0f of %5.0f shuffles (%s), d_shuffles: range [%8.5f %8.5f], mean %8.5f s.d. %8.5f',...
                                            surrogate_count_nestdim(imodel,iadj_ptr_nest,id_calc_type),nshuff,d_calc_types{id_calc_type},...
                                            min(d_shuff_nestdim(imodel,:,iadj_ptr_nest,id_calc_type)),max(d_shuff_nestdim(imodel,:,iadj_ptr_nest,id_calc_type)),...
                                            mean(d_shuff_nestdim(imodel,:,iadj_ptr_nest,id_calc_type)),std(d_shuff_nestdim(imodel,:,iadj_ptr_nest,id_calc_type))));
                                    end
                                  end
                               end %enough dimensions for model?
                            end %iadj_ptr_nest
                        end %if_nestbydim
                     end %if shuffled
                end %model input dimension test
            end
            %final summary
            if opts_geofit.if_summary
                disp(' ');
                disp('summary')
                disp(sprintf('reference dataset dimension: %2.0f',ref_dim));
                disp(sprintf('adjusted  dataset dimension: %2.0f',adj_dim));
                for imodel=1:nmodels
                    model_type=model_types{imodel};
                    model_class=model_types_def.(model_type).class;
                    disp(' ');
                    disp(sprintf('model %20s: class: %s d: %8.5f',model_type,model_class,d(imodel)));
                    if ~isempty(opts_model_used{imodel}) %in case model did not pass dimension test
                        if (nshuff>0)
                            nested_types=model_types_def.(model_type).nested;
                            for inest=1:length(nested_types)
                                nested_type=nested_types{inest};
                                nest_ptr=strmatch(nested_type,model_types,'exact');
                                for id_calc_type=1:length(d_calc_types)
                                    disp(sprintf('nested %23s: d is >= d_shuff for %5.0f of %5.0f shuffles (%s), d_shuffles: range [%8.5f %8.5f], mean %8.5f s.d. %8.5f',...
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
                                                sprintf('dn %20s dim %2.0f',model_type,adj_dim_nest),surrogate_count_nestdim(imodel,iadj_ptr_nest,id_calc_type),nshuff,d_calc_types{id_calc_type},...
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
            end
            r.d=d;
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
        end %dimensions ok
    end %adj_dim_ptr
end %ref_dim_ptr
return
end

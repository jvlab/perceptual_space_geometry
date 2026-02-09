function opts_used=psg_geomodels_plot(results,opts)
% opts_used=psg_geomodels_plot(results,opts) creates summaryy plots of goodness
% of fit of geometric models, including comparison of models with nested models, and comparison of models nested by dimension
%
% results: output of psg_geomodels_fit or psg_geomodels_run containing goodness of fit of one or more geometric models
%     results{ref,adj} contains the analysis for transforming an input set with adj dimensions to an output set with ref dimensions
% 
% opts: options
%   if_nestbymodel_showall: 1 to show all nested-by-model comparisons, 0 for only maximally nested models
%   sig_level: significance level
%   if_showsig: which significance flags to show for d (goodness of fit): 0: none, 1: based on original denom, 2 based on shuffle denom, 3: both (default)
%   if_showquant: 1 to show quantile at significance level sig_level (defaults to 0)
%   colors_models: colors to use for model, used in cyclic order, default- {'k','b','c','m','r',[1 0.5 0],[0.7 0.7 0],'g',[.5 .5 .5],[.5 0 0]}
%   sig_symbols: symbols to mark significant values, sig_symbols{1} for original denom, sig_symbols{2} for shuffle denom, default-{'+','x'}
%   sig_symsize: symbol size for significant values, default=14
%   lw_model: line width for model, default=2
%   lw_nest: line width for nested model, default=2
%   lw_quant: line width for quantile, default=1
%
% opts_used: options used
%    also:
%    opts_used.fig_handles is a cell array of handles to the figures created
%    opts_used.fig_names is a cell array of names of the figures created; no special chars, suitable for use in a file name
%
% This is a modularized version of psg_geomodels_summ, with the following main differences:
%   results must be passed, not read from a file
%   results cannot be a cell array of subsidiary results structures ('mode 2' of psg_geomodels_summ)
%   results may be present on any subset of the (ref,adj) pairs
%
%   See also: RS_SAVE_FIGS, PSG_GEOMODELS_SUMM, PSG_GEOMODELS_FIT, PSG_GEOMODELS_RUN, PSG_GEOMODELS_DEFINE, SURF_AUGVEC, HLID_MDS_COORDS_GEOMODELS.
%
if nargin<=1
    opts=struct();
end
%plot format options
opts=filldefault(opts,'sig_symbols',{'+','x'});
opts=filldefault(opts,'sig_symsize',14);
opts=filldefault(opts,'quant_lines',{'--','-.'});
opts=filldefault(opts,'colors_mn',{'k','b'});
opts=filldefault(opts,'colors_models',{'k','b','c','m','r',[1 0.5 0],[0.7 0.7 0],'g',[.5 .5 .5],[.5 0 0]});
opts=filldefault(opts,'lw_model',2); %line width for a model
opts=filldefault(opts,'lw_nest',2); %line width for a nested model
opts=filldefault(opts,'lw_quant',1); %line width for quantiles
opts=filldefault(opts,'if_omnicolors',1); %1 to use colors from omnibus plots in comparison plots
%
opts=filldefault(opts,'if_nestbymodel_showall',0); %set to 1 for all comparisons
opts=filldefault(opts,'sig_level',0.05);
opts=filldefault(opts,'if_showsig',3); % which significance flags to show(0: none, 1: orig, 2: shuff, 3: both','d',[0 3],3);
opts=filldefault(opts,'if_showquant',0); %1 to show quantile at requested significance level (sig_level)
%
opts=filldefault(opts,'if_log',0);
%
norm_labels={'orig','shuff'}; %denominator used for normalization of d in significance calculations
%
showsigs(1)=mod(opts.if_showsig,2);
showsigs(2)=double(opts.if_showsig>=2);
%
fig_counter=0;
opts.fig_handles=cell(0);
opts.fig_names=cell(0);
%
have_data=zeros(size(results));
for ref_dim=1:size(results,1)
    for adj_dim=1:size(results,2);
        have_data(ref_dim,adj_dim)=isfield(results{ref_dim,adj_dim},'model_types_def');
    end
end
ref_dim_list=find(any(have_data,2)'>0);
adj_dim_list=find(any(have_data,1)>0);
if opts.if_log
    disp('reference set dimension list:')
    disp(ref_dim_list);
    disp('adjusted  set dimension list')
    disp(adj_dim_list);
end
%
r=results{ref_dim_list(1),adj_dim_list(1)};
if isfield(r,'d_shuff')
    nshuff=size(r.d_shuff,2);
else
    nshuff=0;
end
%
model_types_def=r.model_types_def;
n_calc_types=size(r.surrogate_count,3); %typically 2       
if_nestbydim=double(isfield(r,'nestdim_list'));
clear r
%
model_types=model_types_def.model_types;
if opts.if_log
    disp('model types');
    disp(model_types);
    disp(sprintf('nest by dimension analysis present: %1.0f',if_nestbydim));
end
compares=zeros(0,2); %list of nested comparisons (each row is a pointer into model_types for model, and nested model)
%build up complete list of comparisons
for imodel=1:length(model_types)
    model_type=model_types{imodel};
    nested_list=model_types_def.(model_type).nested;
    for inest_ptr=1:length(nested_list)
        nested_type=nested_list{inest_ptr};
        inest=strmatch(nested_type,model_types,'exact');
        compares=[compares;[imodel inest]];
        if opts.if_log
            disp(sprintf('model %30s contains %30s as a nested model',model_type,model_types{inest})); %do this to verify proper match
        end
    end %inset_ptr
    if length(nested_list)==0
        if opts.if_log
            disp(sprintf('model %30s has no nested models',model_type));
        end
    end
end %imodel
compares=sortrows(compares,[1 2],{'ascend' 'descend'}); %sort nested models downward
%find critical comparisons, indicate with a flag in third column
compares=[compares,zeros(size(compares,1),1)];
if opts.if_log
    disp('determining which model nestings are maximal')
end
for icompare=1:size(compares,1)
    imodel=compares(icompare,1);
    inest=compares(icompare,2);
    model_type=model_types{imodel};
    nested_type=model_types{inest};
    %a comparison with a nested model is not critical if that nested model
    %is already contained in an intermediate nested model
    nesteds=compares(find(compares(:,1)==imodel),2); %what this model can be compared to
    nesteds_already=nesteds(nesteds>=inest); %working down the list
    intermed=compares(find(compares(:,2)==inest),1);
    crit=double(isempty(intersect(intermed,nesteds_already)));
    compares(icompare,3)=crit;
    if (opts.if_log)
        disp(sprintf('model %30s contains %30s, maximal nesting: %1.0f',model_type,nested_type,crit))
    end
end %icompare
%
%collect values of d across all models
%
d_models=nan(length(ref_dim_list),length(adj_dim_list),length(model_types)); %d1: ref_dim_ptr, d2: adj_dim_ptr, dim3: model
d_models_nestdim=       nan(length(ref_dim_list),length(adj_dim_list),length(model_types),nshuff,max(adj_dim_list)-1,n_calc_types); %d1: ref_dim_ptr, d2: adj_dim_ptr, dim3: model, d4: shuf, d5: nestdim, d6: n_calc_types
surrogate_count_nestdim=zeros(length(ref_dim_list),length(adj_dim_list),length(model_types)       ,max(adj_dim_list)-1,n_calc_types); %d1: ref_dim_ptr, d2: adj_dim_ptr, dim3: model,           d4: nestdim, d5: n_calc_types
%
for imodel=1:length(model_types)
    model_type=model_types{imodel};
    %
    %collect values of d across all dimensions and models
    %
    for ref_dim_ptr=1:length(ref_dim_list)
        ref_dim=ref_dim_list(ref_dim_ptr);
        for adj_dim_ptr=1:length(adj_dim_list)
            adj_dim=adj_dim_list(adj_dim_ptr);
            rz=results{ref_dim,adj_dim};
            if have_data(ref_dim,adj_dim)
                if adj_dim>=model_types_def.(model_type).min_inputdims
                    d_models(ref_dim_ptr,adj_dim_ptr,imodel)=rz.d(imodel);
                end
                if if_nestbydim
                    if size(rz.nestdim_list,2)>0
                        d_models_nestdim(ref_dim_ptr,adj_dim_ptr,:,:,rz.nestdim_list,:)=...
                            reshape(rz.d_shuff_nestdim,[1 1 length(model_types) nshuff length(rz.nestdim_list) n_calc_types]);
                        surrogate_count_nestdim(ref_dim_ptr,adj_dim_ptr,:,rz.nestdim_list,:)=...
                            reshape(rz.surrogate_count_nestdim,[1 1 length(model_types) length(rz.nestdim_list) n_calc_types]);
                    end
                end %if_nestbydim
            end
        end %adj_dim_ptr
    end %ref_dim_ptr
end
%
%make omnibus plot
%
if length(adj_dim_list)==1
    xlim_adj=adj_dim_list+[-0.5 0.5];
else
    xlim_adj=adj_dim_list([1 end]);
end
if length(ref_dim_list)==1
    ylim_ref=ref_dim_list+[0-0.5 0.5];
else
    ylim_ref=ref_dim_list([1 end]);
end
figname_tstring='all-models';
set(gcf,'Position',[100 100 1200 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name',figname_tstring);
hold on;
for imodel=1:length(model_types)
    h_model=surf_augvec(adj_dim_list,ref_dim_list,d_models(:,:,imodel));
    set(h_model,'FaceColor','none');
    set(h_model,'EdgeColor',opts.colors_models{1+mod(imodel-1,length(opts.colors_models))});
    set(h_model,'LineWidth',opts.lw_model);
end
xlabel('adj dim');
ylabel('ref dim');
set(gca,'XLim',xlim_adj);
set(gca,'YLim',ylim_ref);
set(gca,'XTick',adj_dim_list);
set(gca,'YTick',ref_dim_list);
zlabel('d');
set(gca,'ZLim',[0 1]);
grid on;
box on;
view(3);
hl=legend(strvcat(model_types));
set(hl,'Interpreter','none');
%
fig_counter=fig_counter+1;
opts.fig_handles{1,fig_counter}=gcf;
opts.fig_names{1,fig_counter}=figname_tstring;
%
%plot nested model comparisons
%
disp(' ');
for icompare=1:size(compares,1)
    imodel=compares(icompare,1);
    inest=compares(icompare,2);
    %d_shuff dim 3 is indexed by serial order of nested model
    allnested=unique(compares(find(compares(:,1)==imodel),2));
    inest_dshuff_index=find(allnested==inest); %26Jan26: where to find nested-model data in d_shuff dim 3
    % inest
    % inest_dshuff_index
    %
    model_type=model_types{imodel};
    nested_type=model_types{inest};
    if opts.if_showsig>0
        text_string=sprintf('nshuff %5.0f p=%5.3f  normalization:',nshuff,opts.sig_level);
    else
        text_string=sprintf('nshuff %5.0f',nshuff);
    end
    if compares(icompare,3)==1 | opts.if_nestbymodel_showall==1
        disp(sprintf('model %30s contains %30s: analyzing',model_type,nested_type))
        %
        d_model=d_models(:,:,imodel);
        d_nest=d_models(:,:,inest);
        if_sig=zeros(length(ref_dim_list),length(adj_dim_list),n_calc_types); %d3 is normalization type
        d_shuffs=zeros(length(ref_dim_list),length(adj_dim_list),nshuff,n_calc_types);
        %
        %collect values of d across all dimensions
        nest_ptr=strmatch(nested_type,model_types_def.model_types,'exact'); %added 26Jan26, pointer to surrogates are based on search order
        for ref_dim_ptr=1:length(ref_dim_list)
            ref_dim=ref_dim_list(ref_dim_ptr);
            for adj_dim_ptr=1:length(adj_dim_list)
                adj_dim=adj_dim_list(adj_dim_ptr);
                rz=results{ref_dim,adj_dim};
                if have_data(ref_dim,adj_dim)
                    for calc_type=1:n_calc_types
                        if adj_dim>=model_types_def.(model_type).min_inputdims
                            if_sig(ref_dim_ptr,adj_dim_ptr,calc_type)=double(rz.surrogate_count(imodel,nest_ptr,calc_type)<opts.sig_level*nshuff); %inest->nest_ptr, 26Jan26
                        end
                        d_shuffs(ref_dim_ptr,adj_dim_ptr,:,calc_type)=reshape(rz.d_shuff(imodel,:,inest_dshuff_index,calc_type),[1 1 nshuff]); %inest->inest_dshuff_index, 26Jan26
                    end
                end
            end
        end
        figure;
        figname_tstring=cat(2,'model-',model_type,'-nested-',nested_type);
        figname_tstring=strrep(figname_tstring,' ','-');
        figname_tstring=strrep(figname_tstring,'_','-');
        figname_tstring=strrep(figname_tstring,'procrustes','proc');
        figname_tstring=strrep(figname_tstring,'affine','aff');
        figname_tstring=strrep(figname_tstring,'projective','proj');
        %
        set(gcf,'Position',[100 100 1200 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',figname_tstring);
        hold on;
        h_model=surf_augvec(adj_dim_list,ref_dim_list,d_model);
        set(h_model,'FaceColor','none');
        set(h_model,'LineWidth',opts.lw_model);
        h_nest=surf_augvec(adj_dim_list,ref_dim_list,d_nest);
        set(h_nest,'FaceColor','none');
        set(h_nest,'LineWidth',opts.lw_nest);
        %
        if (opts.if_omnicolors) %use colors from ombnibus ;plots?
            set(h_model,'EdgeColor',opts.colors_models{1+mod(imodel-1,length(opts.colors_models))});
            set(h_nest,'EdgeColor',opts.colors_models{1+mod(inest-1,length(opts.colors_models))});
        else
            set(h_model,'EdgeColor',opts.colors_mn{1});
            set(h_nest,'EdgeColor',opts.colors_mn{2});
        end
        %
        legend_labels=strvcat(model_type,nested_type);
        %plot quantile at significance level
        if opts.if_showquant
            qstring=sprintf('p=%6.3f, den: ',opts.sig_level);
            for calc_type=1:n_calc_types
                if showsigs(calc_type)
                    h_quant=surf_augvec(adj_dim_list,ref_dim_list,quantile(d_shuffs(:,:,:,calc_type),opts.sig_level,3));
                    set(h_quant,'FaceColor','none');
                    set(h_quant,'EdgeColor',get(h_model,'EdgeColor'));
                    set(h_quant,'LineStyle',opts.quant_lines{calc_type});
                    set(h_quant,'LineWidth',opts.lw_quant);
                    legend_labels=strvcat(legend_labels,cat(2,qstring,norm_labels{calc_type}));
                end
            end
        end
        %plot significance flags
        for calc_type=1:n_calc_types
            if showsigs(calc_type)
                for ref_dim_ptr=1:length(ref_dim_list)
                    for adj_dim_ptr=1:length(adj_dim_list)
                        if if_sig(ref_dim_ptr,adj_dim_ptr,calc_type)
                            hsig=plot3(adj_dim_list(adj_dim_ptr),ref_dim_list(ref_dim_ptr),d_model(ref_dim_ptr,adj_dim_ptr),opts.sig_symbols{calc_type});
                            set(hsig,'MarkerSize',opts.sig_symsize);
                            set(hsig,'Color',get(h_model,'EdgeColor')); %color of model
                        end
                    end %adj_dim_ptr
                end %ref_dim_ptr
                text_string=cat(2,text_string,sprintf('%s denom (%s) ',norm_labels{calc_type},opts.sig_symbols{calc_type}));
            end
        end
        xlabel('adj dim');
        ylabel('ref dim');
        set(gca,'XLim',xlim_adj);
        set(gca,'YLim',ylim_ref);
        set(gca,'XTick',adj_dim_list);
        set(gca,'YTick',ref_dim_list);
        zlabel('d');
        set(gca,'ZLim',[0 1]);
        grid on;
        box on;
        view(3);
        hl=legend(legend_labels);
        set(hl,'Interpreter','none');
        %
        axes('Position',[0.01,0.05,0.01,0.01]); %for text
        text(0,0,text_string,'Interpreter','none');
        axis off;
        %
        figname=figname_tstring;
        fig_counter=fig_counter+1;
        opts.fig_handles{1,fig_counter}=gcf;
        opts.fig_names{1,fig_counter}=figname;
    end
end %icompare
%
% plot nested comparisons by dimension
%
if (if_nestbydim)
    disp(' ');
    for imodel=1:length(model_types)
        model_type=model_types{imodel};
        if opts.if_showsig>0
            text_string=sprintf('nshuff %5.0f p=%5.3f  normalization:',nshuff,opts.sig_level);
        else
            text_string=sprintf('nshuff %5.0f',nshuff);
        end
        disp(sprintf('comparing model %s and nested model with lower dimension',model_type))
        d_model=d_models(:,:,imodel);
        if_sig=zeros(length(ref_dim_list),length(adj_dim_list),n_calc_types); %d3 is normalization type
         %
         %collect values of d across all dimensions
         d_nestdim=NaN(length(ref_dim_list),length(adj_dim_list));
         d_shuffs_nestdim=NaN(length(ref_dim_list),length(adj_dim_list),nshuff,n_calc_types);
         for ref_dim_ptr=1:length(ref_dim_list)
             ref_dim=ref_dim_list(ref_dim_ptr);
             for adj_dim_ptr=1:length(adj_dim_list)
                adj_dim=adj_dim_list(adj_dim_ptr);
                if have_data(ref_dim,adj_dim)
                    rz=results{ref_dim,adj_dim};
                    if (adj_dim_ptr>1) 
                        for calc_type=1:n_calc_types
                             if adj_dim_list(adj_dim_ptr-1)>=model_types_def.(model_type).min_inputdims
                                if_sig(ref_dim_ptr,adj_dim_ptr,calc_type)=double(surrogate_count_nestdim(ref_dim_ptr,adj_dim_ptr,imodel,adj_dim_list(adj_dim_ptr-1),calc_type)<opts.sig_level*nshuff);
                             end
                        end %calc_type
                        d_shuffs_nestdim(ref_dim_ptr,adj_dim_ptr,:,:)=reshape(d_models_nestdim(ref_dim_ptr,adj_dim_ptr,imodel,:,adj_dim_list(adj_dim_ptr-1),:),[1 1 nshuff n_calc_types]);
                        d_nestdim(ref_dim_ptr,adj_dim_ptr)=d_model(ref_dim_ptr,adj_dim_ptr-1);
                    end %adj_dim_ptr>1
                end
             end %adj_dim_ptr
         end %ref_dim_ptr
        figure;
        figname_tstring=cat(2,'model-',model_type,'-nested-by-dim');
        figname_tstring=strrep(figname_tstring,' ','-');
        figname_tstring=strrep(figname_tstring,'_','-');
        figname_tstring=strrep(figname_tstring,'procrustes','proc');
        figname_tstring=strrep(figname_tstring,'affine','aff');
        figname_tstring=strrep(figname_tstring,'projective','proj');
        set(gcf,'Position',[100 100 1200 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',figname_tstring);
        hold on;
        h_model=surf_augvec(adj_dim_list,ref_dim_list,d_model);
        set(h_model,'FaceColor','none');
        set(h_model,'LineWidth',opts.lw_model);
        h_nestdim=surf_augvec(adj_dim_list,ref_dim_list,d_nestdim);
        set(h_nestdim,'FaceColor','none');
        set(h_nestdim,'LineWidth',opts.lw_nest);
        set(h_nestdim,'LineStyle',':');
        %
        if (opts.if_omnicolors) %use colors from ombnibus ;plots?
            set(h_model,'EdgeColor',opts.colors_models{1+mod(imodel-1,length(opts.colors_models))});
            set(h_nestdim,'EdgeColor',opts.colors_models{1+mod(imodel-1,length(opts.colors_models))});
        else
            set(h_model,'EdgeColor',opts.colors_mn{1});
            set(h_nestdim,'EdgeColor',opts.colors_mn{1});
        end
        legend_labels=strvcat(model_type,'lower dim');
        %plot quantile at significance level
        if opts.if_showquant
            qstring=sprintf('p=%6.3f, den: ',opts.sig_level);
            for calc_type=1:n_calc_types
                if showsigs(calc_type)
                    h_quant=surf_augvec(adj_dim_list,ref_dim_list,quantile(d_shuffs_nestdim(:,:,:,calc_type),opts.sig_level,3));
                    set(h_quant,'FaceColor','none');
                    set(h_quant,'EdgeColor',get(h_model,'EdgeColor'));
                    set(h_quant,'LineStyle',opts.quant_lines{calc_type});
                    set(h_quant,'LineWidth',opts.lw_quant);
                    legend_labels=strvcat(legend_labels,cat(2,qstring,norm_labels{calc_type}));
                end
            end
        end
         %plot significance flags
         for calc_type=1:n_calc_types
             if showsigs(calc_type)
                 for ref_dim_ptr=1:length(ref_dim_list)
                     for adj_dim_ptr=1:length(adj_dim_list)
                         if if_sig(ref_dim_ptr,adj_dim_ptr,calc_type)
                             hsig=plot3(adj_dim_list(adj_dim_ptr),ref_dim_list(ref_dim_ptr),d_model(ref_dim_ptr,adj_dim_ptr),opts.sig_symbols{calc_type});
                             set(hsig,'Color',get(h_model,'EdgeColor')); %color of model
                             set(hsig,'MarkerSize',opts.sig_symsize);
                         end
                     end %adj_dim_ptr
                 end %ref_dim_ptr
                 text_string=cat(2,text_string,sprintf('%s denom (%s) ',norm_labels{calc_type},opts.sig_symbols{calc_type}));
             end
         end
        xlabel('adj dim');
        ylabel('ref dim');
        set(gca,'XLim',xlim_adj);
        set(gca,'YLim',ylim_ref);
        set(gca,'XTick',adj_dim_list);
        set(gca,'YTick',ref_dim_list);
        zlabel('d');
        set(gca,'ZLim',[0 1]);
        grid on;
        box on;
        view(3);
        hl=legend(legend_labels);
        set(hl,'Interpreter','none');
        %
        axes('Position',[0.01,0.05,0.01,0.01]); %for text
        text(0,0,text_string,'Interpreter','none');
        axis off;
        %
        fig_counter=fig_counter+1;
        opts.fig_handles{1,fig_counter}=gcf;
        opts.fig_names{1,fig_counter}=figname_tstring;
    end
end
opts_used=opts;
return
end


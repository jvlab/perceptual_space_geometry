%psg_geomodels_summ: summarize analysis of geometric models
% omnibus plot, and comparison of models with nested models
%
% The plots are 3d, and surf_augvec is used if adj model or reference model are
% only examined on one dimension.
%
%runs on "results" variable created by psg_geomodels_run
%
%by default, only considers critical comparisons:  
% a comparison with a nested model is not critical if that nested model
% is already contained in an intermediate nested model
%
% 15Jun24: take into account minimal number of dimensins required for a model
% 18Jun24: plot nest-by-dimension data, if present in results
% 08Dec25: a mode in which a single file contains a variable 'r', and r{j,k}.results is to be analyzed
%    (e.g., output of hlid_mds_coords_geomodels)
% 12Dec25: automate saving plots
%
%   See also:   PSG_GEOMODELS_RUN, PSG_GEOMODELS_DEFINE, SURF_AUGVEC, HLID_MDS_COORDS_GEOMODELS.
%
if ~exist('sig_level') sig_level=0.05; end
if ~exist('sig_symbols') sig_symbols={'+','x'}; end
if ~exist('sig_symsize') sig_symsize=14; end
if ~exist('quant_lines') quant_lines={'--','-.'}; end
if ~exist('colors_mn') colors_mn={'k','b'}; end
if ~exist('colors_models') colors_models={'k','b','c','m','r',[1 0.5 0],[0.7 0.7 0],'g'}; end
if ~exist('lw_model') lw_model=2; end %line width for a model
if ~exist('lw_nest') lw_nest=2; end %line width for a nested model
if ~exist('lw_quant') lw_quant=1.5; end %line width for quantiles
norm_labels={'orig','shuff'}; %denominator used for normalization of d in significance calculations
%
mode=getinp('mode: 1 for a single set of results, 2 for multiple ''results'' structures in a single file','d',[1 2],1);
if (mode==1)
    fn=getinp('mat-file from psg_geomodels_run with results to summarize','s',[]);
    results_rows=1;
    results_cols=1;
    rows_anal=1;
    cols_anal=1;
else
    fn=getinp('mat-file containing multiple sets of results to summarize','s',[]);
    results_all=getfield(load(fn),'r');
    results_rows=size(results_all,1);
    results_cols=size(results_all,2);
    disp(sprintf('file contains %2.0f x %2.0f cell array of results',results_rows,results_cols));
    rows_anal=getinp('rows to analyze','d',[1 results_rows],[1:results_rows]);
    cols_anal=getinp('cols to analyze','d',[1 results_cols],[1:results_cols]);
end
if_savefig_close=getinp('1 to save figures and close','d',[0 1]);
if if_savefig_close
    figname_base=getinp('figure name base','s',[],strrep(fn,'.mat',''));
end
if_log=getinp('1 for detailed log','d',[0 1],0);
if_allcompare=getinp('1 for all comparisons (0 just for critical ones)','d',[0 1],0);
if_omnicolors=getinp('1 to use colors from omnibus plots in comparison plots','d',[0 1],1); %changed from 0, 12Dec25
sig_level=getinp('significance level','f',[0 1],sig_level);
if_showsig=getinp('show significance flags (0: none, 1: orig, 2: shuff, 3: both','d',[0 3],3);
showsigs(1)=mod(if_showsig,2);
showsigs(2)=double(if_showsig>=2);
if_showquant=getinp(sprintf('1 to show quantile at p=%6.4f',sig_level),'d',[0 1],0);
%
embed_meth_short=cell(results_rows,results_cols);
for row_ptr=1:length(rows_anal)
    row_anal=rows_anal(row_ptr);
    for col_ptr=1:length(cols_anal)
        col_anal=cols_anal(col_ptr);
        switch mode
            case 1
                results=getfield(load(fn),'results');
                embed_meth_desc='unspecified';
                embed_meth_short{row_anal,col_anal}='std';
            case 2
                results=results_all{row_anal,col_anal}.results;
                embed_meth_desc=results_all{row_anal,col_anal}.method_desc;
                embed_meth_short{row_anal,col_anal}=embed_meth_desc; %create a short string for a file name
                embed_meth_short{row_anal,col_anal}=strrep(embed_meth_short{row_anal,col_anal},'distance ','');
                embed_meth_short{row_anal,col_anal}=strrep(embed_meth_short{row_anal,col_anal},'Euclidean ','euc_');
                embed_meth_short{row_anal,col_anal}=strrep(embed_meth_short{row_anal,col_anal},'via ','');
                embed_meth_short{row_anal,col_anal}=strrep(embed_meth_short{row_anal,col_anal},'SVD','svd_');
                embed_meth_short{row_anal,col_anal}=strrep(embed_meth_short{row_anal,col_anal},'MDS','mds_');
                embed_meth_short{row_anal,col_anal}=strrep(embed_meth_short{row_anal,col_anal},'similarity ','');        
                embed_meth_short{row_anal,col_anal}=strrep(embed_meth_short{row_anal,col_anal},'similarity, ','');        
                embed_meth_short{row_anal,col_anal}=strrep(embed_meth_short{row_anal,col_anal},'cosine ','cos_');
                embed_meth_short{row_anal,col_anal}=strrep(embed_meth_short{row_anal,col_anal},'Pearson ','pea_');
                embed_meth_short{row_anal,col_anal}=strrep(embed_meth_short{row_anal,col_anal},'as ','');
                embed_meth_short{row_anal,col_anal}=strrep(embed_meth_short{row_anal,col_anal},'angle','angle_');
                embed_meth_short{row_anal,col_anal}=strrep(embed_meth_short{row_anal,col_anal},'chord','chord_');
                embed_meth_short{row_anal,col_anal}=strrep(embed_meth_short{row_anal,col_anal},', nosub mean','nosub');
                embed_meth_short{row_anal,col_anal}=strrep(embed_meth_short{row_anal,col_anal},'_nosub mean','_nosub');
                embed_meth_short{row_anal,col_anal}=strrep(embed_meth_short{row_anal,col_anal},', sub mean','sm');
                embed_meth_short{row_anal,col_anal}=strrep(embed_meth_short{row_anal,col_anal},'_sub mean','_sm');
        end
        disp(' ');
        disp('******');
        disp(sprintf(' examining results {%2.0f,%2.0f}, embedding method: %s (tag: %s)',row_anal,col_anal,embed_meth_desc,embed_meth_short{row_anal,col_anal}));
        %
        have_data=zeros(size(results));
        for ref_dim=1:size(results,1)
            for adj_dim=1:size(results,2);
                have_data(ref_dim,adj_dim)=isfield(results{ref_dim,adj_dim},'model_types_def');
            end
        end
        ref_dim_list=find(any(have_data,2)'>0);
        adj_dim_list=find(any(have_data,1)>0);
        disp('reference set dimension list:')
        disp(ref_dim_list);
        disp('adjusted  set dimension list')
        disp(adj_dim_list);
        %
        r=results{ref_dim_list(1),adj_dim_list(1)};
        s=struct;
        switch mode
            case 1
                s.ref_file=r.ref_file;
                s.adj_file=r.adj_file;
                disp(s);
                nshuff=r.nshuff;
                clear s
            case 2
                nshuff=r.opts_geofit.nshuffs;
        end
        model_types_def=r.model_types_def;
        n_calc_types=size(r.surrogate_count,3); %typically 2       
        if_nestbydim=double(isfield(r,'nestdim_list'));
        clear r
        %
        model_types=model_types_def.model_types;
        disp('model types');
        disp(model_types);
        %
        disp(sprintf('nest by dimension analysis present: %1.0f',if_nestbydim));
        compares=zeros(0,3); %list of nested comparisons (as pointers into model_types) and 1 if critical
        %build up complete list of comparisons
        for imodel=1:length(model_types)
            model_type=model_types{imodel};
            nested_list=model_types_def.(model_type).nested;
            if if_log
                disp(' ');
                disp(sprintf('for model %s, nested models are',model_type));
                disp(nested_list);
            end
            for inest_ptr=1:length(nested_list)
                nested_type=nested_list{inest_ptr};
                inest=strmatch(nested_type,model_types,'exact');
                compares=[compares;[imodel inest]];
                if if_log
                    disp(sprintf('model %s contains %s as a nested model',model_type,model_types{inest})); %do this to verify proper match
                end
            end %inset_ptr
        end %imodel
        compares=sortrows(compares,[1 2],{'ascend' 'descend'}); %sort nested models downward
        %find critical comparisons
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
            if (if_log)
                disp(sprintf('for model %s and nested model %s: critical=%1.0f',model_type,nested_type,crit))
            end
        end %icompare
        %
        %collect values of d across all models
        %
        d_models=nan(length(ref_dim_list),length(adj_dim_list),length(model_types)); %d1: ref_dim_ptr, d2: adj_dim_ptr, dim3: model
        d_models_nestdim=       nan(length(ref_dim_list),length(adj_dim_list),length(model_types),nshuff,length(adj_dim_list)-1,n_calc_types); %d1: ref_dim_ptr, d2: adj_dim_ptr, dim3: model, d4: shuf, d5: nestdim, d6: n_calc_types
        surrogate_count_nestdim=zeros(length(ref_dim_list),length(adj_dim_list),length(model_types)       ,length(adj_dim_list)-1,n_calc_types); %d1: ref_dim_ptr, d2: adj_dim_ptr, dim3: model,           d4: nestdim, d5: n_calc_types
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
                    if adj_dim>=model_types_def.(model_type).min_inputdims
                        d_models(ref_dim_ptr,adj_dim_ptr,imodel)=results{ref_dim,adj_dim}.d(imodel);
                    end
                    if if_nestbydim & adj_dim_ptr>1
                        if adj_dim_list(adj_dim_ptr-1)>=model_types_def.(model_type).min_inputdims
                            d_models_nestdim(ref_dim_ptr,adj_dim_ptr,:,:,1:adj_dim_ptr-1,:)=...
                                reshape(results{ref_dim,adj_dim}.d_shuff_nestdim,[1 1 length(model_types),nshuff,adj_dim_ptr-1,n_calc_types]);
                            surrogate_count_nestdim(ref_dim_ptr,adj_dim_ptr,:,1:adj_dim_ptr-1,:)=...
                                reshape(results{ref_dim,adj_dim}.surrogate_count_nestdim,[1 1 length(model_types),adj_dim_ptr-1,n_calc_types]);
                        end
                    end %if_nestbydim
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
        figure;
        tstring_omni='all models';
        figname_tstring='allmodels';
        set(gcf,'Position',[100 100 1200 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',tstring_omni);
        hold on;
        for imodel=1:length(model_types)
            h_model=surf_augvec(adj_dim_list,ref_dim_list,d_models(:,:,imodel));
            set(h_model,'FaceColor','none');
            set(h_model,'EdgeColor',colors_models{1+mod(imodel-1,length(colors_models))});
            set(h_model,'LineWidth',lw_model);
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
        axes('Position',[0.01,0.02,0.01,0.01]); %for text
        text(0,0,cat(2,tstring_omni,' ',fn),'Interpreter','none');
        axis off;
        axes('Position',[0.5,0.02,0.01,0.01]); %for text
        text(0,0,sprintf('embed meth: %s',embed_meth_desc),'Interpreter','none');
        axis off;
        if if_savefig_close
            figname=cat(2,figname_base,'_',figname_tstring,'_',embed_meth_short{row_anal,col_anal});
            disp(sprintf('saving figure as %s and then closing',figname));
            savefig(figname);
            close (gcf);
        end
        %
        %plot nested model comparisons
        %
        disp(' ');
        for icompare=1:size(compares,1)
            imodel=compares(icompare,1);
            inest=compares(icompare,2);
            model_type=model_types{imodel};
            nested_type=model_types{inest};
            if if_showsig>0
                text_string=sprintf('nshuff %5.0f p=%5.3f  normalization:',nshuff,sig_level);
            else
                text_string=sprintf('nshuff %5.0f',nshuff);
            end
            if compares(icompare,3)==1 | if_allcompare==1
                disp(sprintf('comparing model %s and nested model %s',model_type,nested_type))
                %
                d_model=d_models(:,:,imodel);
                d_nest=d_models(:,:,inest);
                if_sig=zeros(length(ref_dim_list),length(adj_dim_list),n_calc_types); %d3 is normalization type
                d_shuffs=zeros(length(ref_dim_list),length(adj_dim_list),nshuff,n_calc_types);
                %
                %collect values of d across all dimensions
                for ref_dim_ptr=1:length(ref_dim_list)
                    ref_dim=ref_dim_list(ref_dim_ptr);
                    for adj_dim_ptr=1:length(adj_dim_list)
                        adj_dim=adj_dim_list(adj_dim_ptr);
                        for calc_type=1:n_calc_types
                            if adj_dim>=model_types_def.(model_type).min_inputdims
                                if_sig(ref_dim_ptr,adj_dim_ptr,calc_type)=double(results{ref_dim,adj_dim}.surrogate_count(imodel,inest,calc_type)<sig_level*nshuff);
                            end
                            d_shuffs(ref_dim_ptr,adj_dim_ptr,:,calc_type)=reshape(results{ref_dim,adj_dim}.d_shuff(imodel,:,inest,calc_type),[1 1 nshuff]);
                        end
                    end
                end
                figure;
                tstring=cat(2,'model type: ',model_type,', nested model:',nested_type);
                figname_tstring=cat(2,'model-',model_type,'-nested-',nested_type);
                figname_tstring=strrep(figname_tstring,' ','-');
                figname_tstring=strrep(figname_tstring,'_','-');
                figname_tstring=strrep(figname_tstring,'procrustes','proc');
                figname_tstring=strrep(figname_tstring,'affine','aff');
                figname_tstring=strrep(figname_tstring,'projective','proj');
                %
                set(gcf,'Position',[100 100 1200 800]);
                set(gcf,'NumberTitle','off');
                set(gcf,'Name',tstring);
                hold on;
                h_model=surf_augvec(adj_dim_list,ref_dim_list,d_model);
                set(h_model,'FaceColor','none');
                set(h_model,'LineWidth',lw_model);
                h_nest=surf_augvec(adj_dim_list,ref_dim_list,d_nest);
                set(h_nest,'FaceColor','none');
                set(h_nest,'LineWidth',lw_nest);
                %
                if (if_omnicolors) %use colors from ombnibus ;plots?
                    set(h_model,'EdgeColor',colors_models{1+mod(imodel-1,length(colors_models))});
                    set(h_nest,'EdgeColor',colors_models{1+mod(inest-1,length(colors_models))});
                else
                    set(h_model,'EdgeColor',colors_mn{1});
                    set(h_nest,'EdgeColor',colors_mn{2});
                end
                %
                legend_labels=strvcat(model_type,nested_type);
                %plot quantile at significance level
                if if_showquant
                    qstring=sprintf('p=%6.3f, den: ',sig_level);
                    for calc_type=1:n_calc_types
                        if showsigs(calc_type)
                            h_quant=surf_augvec(adj_dim_list,ref_dim_list,quantile(d_shuffs(:,:,:,calc_type),sig_level,3));
                            set(h_quant,'FaceColor','none');
                            set(h_quant,'EdgeColor',get(h_model,'EdgeColor'));
                            set(h_quant,'LineStyle',quant_lines{calc_type});
                            set(h_quant,'LineWidth',lw_quant);
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
                                    hsig=plot3(adj_dim_list(adj_dim_ptr),ref_dim_list(ref_dim_ptr),d_model(ref_dim_ptr,adj_dim_ptr),sig_symbols{calc_type});
                                    set(hsig,'MarkerSize',sig_symsize);
                                    set(hsig,'Color',get(h_model,'EdgeColor')); %color of model
                                end
                            end %adj_dim_ptr
                        end %ref_dim_ptr
                        text_string=cat(2,text_string,sprintf('%s denom (%s) ',norm_labels{calc_type},sig_symbols{calc_type}));
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
                axes('Position',[0.01,0.02,0.01,0.01]); %for text
                text(0,0,cat(2,tstring,' ',fn),'Interpreter','none');
                axis off;
                axes('Position',[0.01,0.05,0.01,0.01]); %for text
                text(0,0,text_string,'Interpreter','none');
                axis off;
                axes('Position',[0.5,0.05,0.01,0.01]); %for text
                text(0,0,sprintf('embed meth: %s',embed_meth_desc),'Interpreter','none');
                axis off;
                %
                if if_savefig_close
                    figname=cat(2,figname_base,'_',figname_tstring,'_',embed_meth_short{row_anal,col_anal});
                    disp(sprintf('saving figure as %s and then closing',figname));
                    savefig(figname);
                    close (gcf);
                end
            end
        end %icompare
        %
        % plot nested comparisons by dimension
        %
        if (if_nestbydim)
            disp(' ');
            for imodel=1:length(model_types)
                model_type=model_types{imodel};
                if if_showsig>0
                    text_string=sprintf('nshuff %5.0f p=%5.3f  normalization:',nshuff,sig_level);
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
                        if (adj_dim_ptr>1) 
                            for calc_type=1:n_calc_types
                                 if adj_dim_list(adj_dim_ptr-1)>=model_types_def.(model_type).min_inputdims
                                    if_sig(ref_dim_ptr,adj_dim_ptr,calc_type)=double(results{ref_dim,adj_dim}.surrogate_count_nestdim(imodel,adj_dim_ptr-1,calc_type)<sig_level*nshuff);
                                 end
                            end %calc_type
                            d_shuffs_nestdim(ref_dim_ptr,adj_dim_ptr,:,:)=reshape(d_models_nestdim(ref_dim_ptr,adj_dim_ptr,imodel,:,adj_dim_ptr-1,:),[1 1 nshuff n_calc_types]);
                            d_nestdim(ref_dim_ptr,adj_dim_ptr)=d_model(ref_dim_ptr,adj_dim_ptr-1);
                        end %adj_dim_ptr>1
                     end %adj_dim_ptr
                 end %ref_dim_ptr
                figure;
                tstring=cat(2,'model type: ',model_type,', nested by dimension');
                figname_tstring=cat(2,'model-',model_type,'-nested-by-dim');
                figname_tstring=strrep(figname_tstring,' ','-');
                figname_tstring=strrep(figname_tstring,'_','-');
                figname_tstring=strrep(figname_tstring,'procrustes','proc');
                figname_tstring=strrep(figname_tstring,'affine','aff');
                figname_tstring=strrep(figname_tstring,'projective','proj');
                set(gcf,'Position',[100 100 1200 800]);
                set(gcf,'NumberTitle','off');
                set(gcf,'Name',tstring);
                hold on;
                h_model=surf_augvec(adj_dim_list,ref_dim_list,d_model);
                set(h_model,'FaceColor','none');
                set(h_model,'LineWidth',lw_model);
                h_nestdim=surf_augvec(adj_dim_list,ref_dim_list,d_nestdim);
                set(h_nestdim,'FaceColor','none');
                set(h_nestdim,'LineWidth',lw_nest);
                set(h_nestdim,'LineStyle',':');
                %
                if (if_omnicolors) %use colors from ombnibus ;plots?
                    set(h_model,'EdgeColor',colors_models{1+mod(imodel-1,length(colors_models))});
                    set(h_nestdim,'EdgeColor',colors_models{1+mod(imodel-1,length(colors_models))});
                else
                    set(h_model,'EdgeColor',colors_mn{1});
                    set(h_nestdim,'EdgeColor',colors_mn{1});
                end
                legend_labels=strvcat(model_type,'lower dim');
                %plot quantile at significance level
                if if_showquant
                    qstring=sprintf('p=%6.3f, den: ',sig_level);
                    for calc_type=1:n_calc_types
                        if showsigs(calc_type)
                            h_quant=surf_augvec(adj_dim_list,ref_dim_list,quantile(d_shuffs_nestdim(:,:,:,calc_type),sig_level,3));
                            set(h_quant,'FaceColor','none');
                            set(h_quant,'EdgeColor',get(h_model,'EdgeColor'));
                            set(h_quant,'LineStyle',quant_lines{calc_type});
                            set(h_quant,'LineWidth',lw_quant);
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
                                     hsig=plot3(adj_dim_list(adj_dim_ptr),ref_dim_list(ref_dim_ptr),d_model(ref_dim_ptr,adj_dim_ptr),sig_symbols{calc_type});
                                     set(hsig,'Color',get(h_model,'EdgeColor')); %color of model
                                     set(hsig,'MarkerSize',sig_symsize);
                                 end
                             end %adj_dim_ptr
                         end %ref_dim_ptr
                         text_string=cat(2,text_string,sprintf('%s denom (%s) ',norm_labels{calc_type},sig_symbols{calc_type}));
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
                axes('Position',[0.01,0.02,0.01,0.01]); %for text
                text(0,0,cat(2,tstring,' ',fn),'Interpreter','none');
                axis off;
                axes('Position',[0.01,0.05,0.01,0.01]); %for text
                text(0,0,text_string,'Interpreter','none');
                axis off;
                axes('Position',[0.5,0.05,0.01,0.01]); %for text
                text(0,0,sprintf('embed meth: %s',embed_meth_desc),'Interpreter','none');
                axis off;
                if if_savefig_close
                    figname=cat(2,figname_base,'_',figname_tstring,'_',embed_meth_short{row_anal,col_anal});
                    disp(sprintf('saving figure as %s and then closing',figname));
                    savefig(figname);
                    close (gcf);
                end
            end
        end
    end %col
end %row

%psg_geomodels_summ: summarize analysis of geometric models
% omnibus plot, and comparisoin of models with nested models
%
% The plots are 3d, and assume that both adj model and reference model are
% examined for more than one dimension.
%   to do: make it work for one-dimension analyses, with 
%   hs=surf([4 4],[1:5],repmat(rand(5,1),1,2));set(hs,'FaceColor','none') or similar
%
%runs on "results" variable created by psg_geomodels_run
%
%optionally only considers critical comparisons:  
% a comparison with a nested model is not critical if that nested model
% is already contained in an intermediate nested model
%
%   See also:   PSG_GEOMODELS_RUN,  PSG_GEOMODELS_DEFINE.
%
if ~exist('sig_level') sig_level=0.05; end
if ~exist('sig_symbols') sig_symbols={'+','x'}; end
if ~exist('colors_mn') colors_mn={'k','b'}; end
if ~exist('colors_models') colors_models={'k','b','c','m','r',[1 0.5 0],[0.7 0.7 0],'g'}; end
%
fn=getinp('mat-file from psg_geomodels_run with results to summarize','s',[]);
if_log=getinp('1 for detailed log','d',[0 1],0);
if_allcompare=getinp('1 for all comparisons (0 just for critical ones)','d',[0 1],0);
sig_level=getinp('significance level','f',[0 1],sig_level);
%
results=getfield(load(fn),'results');
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
s.ref_file=r.ref_file;
s.adj_file=r.adj_file;
disp(s);
nshuff=r.nshuff;
model_types_def=r.model_types_def;
clear r s
model_types=model_types_def.model_types;
disp('model types');
disp(model_types);
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
d_models=zeros(length(ref_dim_list),length(adj_dim_list),length(model_types));
for imodel=1:length(model_types)
    model_type=model_types{imodel};
    %
    %collect values of d across all dimensions and models
    %
    for ref_dim_ptr=1:length(ref_dim_list)
        ref_dim=ref_dim_list(ref_dim_ptr);
        for adj_dim_ptr=1:length(adj_dim_list)
            adj_dim=adj_dim_list(adj_dim_ptr);
            d_models(ref_dim_ptr,adj_dim_ptr,imodel)=results{ref_dim,adj_dim}.d(imodel);
        end
    end
end
%
%make omnibus plot
%
figure;
tstring_omni='all models';
set(gcf,'Position',[100 100 1200 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name',tstring_omni);
hold on;
for imodel=1:length(model_types)
    h_model=surf(adj_dim_list,ref_dim_list,d_models(:,:,imodel));
    set(h_model,'FaceColor','none');
    set(h_model,'EdgeColor',colors_models{1+mod(imodel-1,length(colors_models))});
end
xlabel('adj dim');
ylabel('ref dim');
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
%
%do comparisons
%
for icompare=1:size(compares,1)
    imodel=compares(icompare,1);
    inest=compares(icompare,2);
    model_type=model_types{imodel};
    nested_type=model_types{inest};
    if compares(icompare,3)==1 | if_allcompare==1
        disp(sprintf('comparing model %s and nested model %s',model_type,nested_type))
        %
        d_model=d_models(:,:,imodel);
        d_nest=zeros(length(ref_dim_list),length(adj_dim_list));
        if_sig=zeros(length(ref_dim_list),length(adj_dim_list),2); %d3 is normalization type      
        %
        %collect values of d across all dimensions
        for ref_dim_ptr=1:length(ref_dim_list)
            ref_dim=ref_dim_list(ref_dim_ptr);
            for adj_dim_ptr=1:length(adj_dim_list)
                adj_dim=adj_dim_list(adj_dim_ptr);
%                d_model(ref_dim_ptr,adj_dim_ptr)=results{ref_dim,adj_dim}.d(imodel);
                d_nest(ref_dim_ptr,adj_dim_ptr)=results{ref_dim,adj_dim}.d(inest);
                for norm_type=1:2
                    if_sig(ref_dim_ptr,adj_dim_ptr,norm_type)=double(results{ref_dim,adj_dim}.surrogate_count(imodel,inest,norm_type)<sig_level*nshuff);
                end
            end
        end
        figure;
        tstring=cat(2,'model type: ',model_type,', nested model:',nested_type);
        set(gcf,'Position',[100 100 1200 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',tstring);
        hold on;
        h_model=surf(adj_dim_list,ref_dim_list,d_model);
        set(h_model,'FaceColor','none');
        set(h_model,'EdgeColor',colors_mn{1});
        h_nest=surf(adj_dim_list,ref_dim_list,d_nest);
        set(h_nest,'FaceColor','none');
        set(h_nest,'EdgeColor',colors_mn{2});
        %plot significance flags
        for ref_dim_ptr=1:length(ref_dim_list)
            for adj_dim_ptr=1:length(adj_dim_list)
                for norm_type=1:2
                    if if_sig(ref_dim_ptr,adj_dim_ptr,norm_type)
                        hsig=plot3(adj_dim_list(adj_dim_ptr),ref_dim_list(ref_dim_ptr),d_model(ref_dim_ptr,adj_dim_ptr),sig_symbols{norm_type});
                        set(hsig,'Color',colors_mn{1}); %color of model
                    end
                end
            end
        end
        xlabel('adj dim');
        ylabel('ref dim');
        set(gca,'XTick',adj_dim_list);
        set(gca,'YTick',ref_dim_list);
        zlabel('d');
        set(gca,'ZLim',[0 1]);
        grid on;
        box on;
        view(3);
        hl=legend(strvcat(model_type,nested_type));
        set(hl,'Interpreter','none');
        %
        axes('Position',[0.01,0.02,0.01,0.01]); %for text
        text(0,0,cat(2,tstring,' ',fn),'Interpreter','none');
        axis off;
        axes('Position',[0.01,0.05,0.01,0.01]); %for text
        text(0,0,sprintf('nshuff %5.0f p=%5.3f (normalization: orig denom (%s) shuff denom(%s)',...
            nshuff,sig_level,sig_symbols{1},sig_symbols{2}),...
            'Interpreter','none');
        axis off;

    end
end %icompare

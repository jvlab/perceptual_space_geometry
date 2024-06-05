%psg_geomodels_summ: summarize analysis of geometric models
%
%runs on "results" variable created by psg_geomodels_run
%
%optionally only considers critical comparisons:  
% a comparison with a nested model is not critical if that nested model
% is already contained in an intermediate nested model
%
%   See also:   PSG_GEOMODELS_RUN,  PSG_GEOMODELS_DEFINE.
fn=getinp('mat-file from psg_geomodels_run with results to summarize','s',[]);
if_log=getinp('1 for detailed log','d',[0 1],0);
if_allcompare=getinp('1 for all comparisons (0 just for critical ones)','d',[0 1],0);
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
%do comparisons
for icompare=1:size(compares,1)
    imodel=compares(icompare,1);
    inest=compares(icompare,2);
    model_type=model_types{imodel};
    nested_type=model_types{inest};
    if compares(icompare,3)==1 | if_allcompare==1
        disp(sprintf('comparing model %s and nested model %s',model_type,nested_type))
    end
end %icompare



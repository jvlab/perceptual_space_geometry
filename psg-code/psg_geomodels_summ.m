%psg_geomodels_summ: summarize analysis of geometric models
%
%runs on "results" variable created by psg_geomodels_run
%
fn=getinp('mat-file from psg_geomodels_run with results to summarize','s',[]);
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
for imodel=1:length(model_types)
    model_type=model_types{imodel};
    nested_list=model_types_def.(model_type).nested;
    disp(sprintf(' for model %s, nested models are',model_type))
    disp(nested_list)
end

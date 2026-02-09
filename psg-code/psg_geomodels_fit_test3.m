%psg_geomodels_fit_test3: run psg_geomodels_fit with hlid files, comparing psg and rs
%
%this is a test dataset in which the aug file is nested, so both methods of testing nest-by-dimension should give very similar results
if ~exist('ref_file') ref_file='../../OlfNav/HongLab/data/kc_soma_nls/hlid_odor17_coords_kc_soma_nls_pooled.mat'; end
if ~exist('adj_file') adj_file='../../OlfNav/HongLab/data/orn_terminals/hlid_odor17_coords_orn_terminals7sets_pooled.mat'; end
hlid_setup;
opts_read.input_type=1;
opts_read.if_auto=1;
[sets_ref,ds_ref,sas_ref,rayss_ref,opts_read_used_ref]=psg_get_coordsets(setfield(opts_read,'data_fullname_def',ref_file),[],[],1);
[sets_adj,ds_adj,sas_adj,rayss_adj,opts_read_used_adj]=psg_get_coordsets(setfield(opts_read,'data_fullname_def',adj_file),[],[],1);
% 
if ~exist('opts_geofit') opts_geofit=struct; end
%
disp('for standard run, do not excluse any model types.');
opts_geofit.model_types_def=psg_geomodels_define(1);
%
opts_geofit=filldefault(opts_geofit,'ref_dim_list',[1:4]);
opts_geofit=filldefault(opts_geofit,'adj_dim_list',[1:4]);
opts_geofit=filldefault(opts_geofit,'if_center',0); %non-default
opts_geofit=filldefault(opts_geofit,'if_frozen',1);
opts_geofit=filldefault(opts_geofit,'if_log',1);
opts_geofit=filldefault(opts_geofit,'if_summary',1);
opts_geofit=filldefault(opts_geofit,'nshuffs',5);
opts_geofit=filldefault(opts_geofit,'if_nestbydim',1);
%
%
data_in=struct;
data_in.ds=ds_adj;
data_in.sas=sas_adj;
data_in.sets=sets_adj;
%
data_out=struct;
data_out.ds=ds_ref;
data_out.sas=sas_ref;
data_out.sets=sets_ref;
%
aux=struct;
aux.opts_geof=rmfield(rmfield(rmfield(opts_geofit,'model_types_def'),'ref_dim_list'),'adj_dim_list');
aux.opts_geof.dim_max_in=max(opts_geofit.adj_dim_list);
aux.opts_geof.dim_max_out=max(opts_geofit.ref_dim_list);
aux.opts_geof.model_list=opts_geofit.model_types_def.model_types;
aux.opts_geof.dimpairs_method='all';
aux.opts_geof
[rs,xs,aux_out]=rs_geofit(data_in,data_out,aux);

%
[results_psg,opts_geofit_used_psg]=psg_geomodels_fit(ds_ref{1},ds_adj{1},opts_geofit);
%
% for km=1:length(opts_geofit.model_types_def.model_types)
%     disp(sprintf(' model type %2.0f: %s',km,opts_geofit.model_types_def.model_types{km}));
%     for r=opts_geofit.ref_dim_list
%         for a=opts_geofit.adj_dim_list
%             md=max(max(max(abs(results_nbdpos{r,a}.d_shuff(km,:,:,:)-results_nbdneg{r,a}.d_shuff(km,:,:,:)))));
%             mc=max(max(abs(results_nbdpos{r,a}.surrogate_count(km,:,:)-results_nbdneg{r,a}.surrogate_count(km,:,:))));
%             md_nest=max(max(max(abs(results_nbdpos{r,a}.d_shuff_nestdim(km,:,:,:)-results_nbdneg{r,a}.d_shuff_nestdim(km,:,:,:)))));
%             mc_nest=max(max(abs(results_nbdpos{r,a}.surrogate_count_nestdim(km,:,:)-results_nbdneg{r,a}.surrogate_count_nestdim(km,:,:))));
%             disp(sprintf(' d_ref %1.0f d_adj %1.0f max dev d_shuff %10.8f (counts: %2.0f) max dev d_shuff_nest %10.8f (counts,%2.0f)',r,a,md,mc,md_nest,mc_nest));
%         end
%     end
% end
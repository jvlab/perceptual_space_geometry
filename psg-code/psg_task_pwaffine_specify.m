%psg_task_pwaffine_specify: run psg_geomodels_fit on task consensus files, with fixed cutplanes for pwaffine model
%
if ~exist('ref_file') ref_file='psg_data/bgca3pt_coords_6subjs_consensus.mat'; end
if ~exist('adj_file') adj_file=ref_file; end
%
opts_read=struct();
disp('enter the reference coordinate set');
[sets_ref,ds_ref,sas_ref,rayss_ref,opts_read_used_ref]=psg_get_coordsets(setfield(opts_read,'data_fullname_def',ref_file),[],[],1);
%
disp(' ');
disp('enter the coordinate set to be adjusted (typically a qform model)');
opts_read.setup_fullname_def=opts_read_used_ref{1}.setup_fullname;
[sets_adj,ds_adj,sas_adj,rayss_adj,opts_read_used_adj]=psg_get_coordsets(setfield(opts_read,'data_fullname_def',adj_file),[],[],1);
% 
if ~exist('opts_geofit') opts_geofit=struct; end
%
disp('for standard run, include procrustes_noscale_offset,procrustes_scale_offset,affine_offset,projective,pwaffine.');
opts_geofit.model_types_def=psg_geomodels_define(1);
%
if ~exist('dim_list') dim_list=4; end %just dimension 4
opts_geofit=filldefault(opts_geofit,'ref_dim_list',[4]);
opts_geofit=filldefault(opts_geofit,'adj_dim_list',[4]);
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
aux.opts_geof=rmfield(opts_geofit,'model_types_def');
aux.opts_geof.dim_max_in=max(opts_geofit.adj_dim_list);
aux.opts_geof.dim_max_out=max(opts_geofit.ref_dim_list);
aux.opts_geof.model_list=opts_geofit.model_types_def.model_types;
aux.opts_geof.dimpairs_method='list';
aux.opts_geof.dimpairs_list=[dim_list(:) dim_list(:)];
aux.opts_geof
[rs,xs,aux_out]=rs_geofit(data_in,data_out,aux);
results_rs=rs{1}.results;
%
[results_psg,opts_geofit_used_psg]=psg_geomodels_fit(ds_ref{1},ds_adj{1},opts_geofit);
%

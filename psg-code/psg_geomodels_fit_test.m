%psg_geomodels_fit_test: run psg_geomodels_fit with the built-in files for psg_geomdels_run, and compare
%with outputs of psg_geomodels_run -- NOTE that psg_geomodels_run should be
%run with these same files, but NOT in built-in mode, as reading the files
%alphabetizes the stimuli but using the built-in mode does not.
%
if ~exist('ref_file') ref_file='./psg_data/bgca3pt_coords_MC_sess01_10.mat'; end
if ~exist('adj_file') adj_file='./psg_data/bgca3pt_coords_BL_sess01_10.mat'; end
if ~exist('opts_read') opts_read=struct; end
if ~exist('opts_rays') opts_rays=struct; end
if ~exist('opts_qpred') opts_qpred=struct; end
opts_read.input_type=1;
opts_read.if_auto=1;
[sets_ref,ds_ref,sas_ref,rayss_ref,opts_read_used_ref]=psg_get_coordsets(setfield(opts_read,'data_fullname_def',ref_file),opts_rays,opts_qpred,1);
[sets_adj,ds_adj,sas_adj,rayss_adj,opts_read_used_adj]=psg_get_coordsets(setfield(opts_read,'data_fullname_def',adj_file),opts_rays,opts_qpred,1);
% 
%options to match tests of psg_geomodels_run_test_02Dec25a.txt (and 13Jun24a.txt)
if ~exist('opts_geofit') opts_geofit=struct; end
%
disp('for standard run, do not excluse any model types.');
opts_geofit.model_types_def=psg_geomodels_define(1);
%
opts_geofit=filldefault(opts_geofit,'ref_dim_list',[1 3 4]);
opts_geofit=filldefault(opts_geofit,'adj_dim_list',[1 2 4]);
opts_geofit=filldefault(opts_geofit,'if_center',1);
opts_geofit=filldefault(opts_geofit,'if_frozen',1);
opts_geofit=filldefault(opts_geofit,'if_log',1);
opts_geofit=filldefault(opts_geofit,'if_summary',1);
opts_geofit=filldefault(opts_geofit,'nshuffs',5);
opts_geofit=filldefault(opts_geofit,'if_nestbydim',1);
%
[results,opts_geofit_used]=psg_geomodels_fit(ds_ref{1},ds_adj{1},opts_geofit);

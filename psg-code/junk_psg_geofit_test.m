% rs_geofit_test: test rs_geofit
%
%  See also:  RS_GEOFIT, RS_BENCHMARK_COMPARE, RS_SAVE_MAT.
%
test_descs='transforming three binary texture coordinate files to three others, mean, bogus,and procrustes models, maximal nesting, explicit dimensions';
filenames_in={'../rs/src/samples/bwtextures/bgca3pt_coords_BL_sess01_10.mat','../rs/src/samples/bwtextures/bgca3pt_coords_MC_sess01_10.mat','../rs/src/samples/bwtextures/bgca3pt_coords_SN_sess01_10.mat'};
aux_ins.opts_read=setfields(struct(),{'input_type','if_auto','if_log'},{1,1,1});
aux_ins.nsets=3;
filenames_out={'../rs/src/samples/bwtextures/bgca3pt_coords_BL-br_sess01_10.mat','../rs/src/samples/bwtextures/bgca3pt_coords_MC-br_sess01_10.mat','../rs/src/samples/bwtextures/bgca3pt_coords_SN-br_sess01_10.mat'};
aux_outs=aux_ins;
auxs.opts_geof.model_list={'mean','procrustes_noscale_nooffset','procrustes_scale_nooffset','procrustes_noscale_offset','procrustes_scale_offset'};
auxs.opts_geof.if_stats=1;
auxs.opts_geof.if_nestbymodel=-1;
auxs.opts_geof.if_nestbydim=-1; %nest by dimension, with pca
auxs.opts_geof.dimpairs_method='list';
auxs.opts_geof.dimpairs_list=[2 2;2 3;2 4;3 2;3 3;3 4;4 3;3 3;4 5];
auxs.opts_geof.nshuffs=3;
auxs.opts_geof.if_fit_log=1;
%
%read the input and output data structures
aux_ins.opts_read.if_log=0;
[data_ins,aux_ins]=rs_get_coordsets(filenames_in,aux_ins);
%
aux_outs.opts_read.if_log=1;
[data_outs,aux_outs]=rs_get_coordsets(filenames_out,aux_outs);
%
[gfs,xs,aux_geofits]=rs_geofit(data_ins,data_outs,auxs);
%

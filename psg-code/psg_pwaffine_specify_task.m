%psg_pwaffine_specify_task: run psg_geomodels_fit on task consensus files, with fixed cutplanes for pwaffine model
%looks at only the 4-dim model, steps through the four possible cutplanes
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
dim_model=4;
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
aux.opts_geof.dim_max_in=dim_model;
aux.opts_geof.dim_max_out=dim_model;
aux.opts_geof.dimpairs_method='list';
aux.opts_geof.dimpairs_list=[dim_model dim_model];
aux.opts_geof.model_list={'procrustes_noscale_offset','procrustes_scale_offset','affine_offset','projective','pwaffine'};
aux.opts_geof.nshuffs=20;
aux.opts_geof.vcut_specified={[1],[1 0],[1 0 0],[0 0 0 1]}; %change last entry to [0 1 0 0], etc, to cut on other axes
aux.opts_geof.acut_specified={0,0,0,0}; %all cuts at 0
%aux.opts_geof
disp('*****************');
results_cut=cell(1,dim_model);
for cut_axis=1:4
    disp(sprintf('examining cut axis %2.0f in %s',cut_axis,opts_read_used_ref{1}.setup_fullname));
    aux.opts_geof.vcut_specified{4}={find([1:4]==cut_axis)}; %set cut dimension orthogonal to the cut axis
    [rs,xs,aux_out]=rs_geofit(data_in,data_out,aux);
    results_cut{cut_axis}=rs{1}.gf;
end


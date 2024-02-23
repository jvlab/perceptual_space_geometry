% psg_vss24_task_fig_pool: create a figure showing the pooled datasets for the bc55 paradigm
% 
% Pooled datasets must have been previously created, e.g., by psg_align_knit_demo
% 
% See also:  PSG_GET_COORDSETS, PSG_PLOTCOORDS, PSG_ALIGN_COORDSETS, PROCRUSTES_CONSENSUS.
%
opts_read=struct;
opts_rays=struct;
opts_qpred=struct;
if getinp('1 for interactive mode','d',[0 1])
    disp('enter two datasets: first is experimental data, second is model')
    [sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred,2);
else
    opts_read.input_type=[1 2]; %first dataset is data, second is model
    opts_read.data_fullnames={'./psg_data/bc55pooled_coords_MC.mat',[]};
    opts_read.setup_fullnames={'./psg_data/bc55pooled9.mat','./psg_data/bc55pooled9.mat'};
    opts_read.if_auto=1; 
    nsets=length(opts_read.input_type);
    [sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred,nsets);
end
if ~exist('rings_select') rings_select=4; end %which rings to show
rays_mod=rmfield(rayss{1},'rings');
for iring=1:length(rings_select)
    rays_mod.rings{iring}=rayss{1}.rings{rings_select(iring)};
end
rays_mod.nrings=length(rings_select);
%
if ~exist('dim_align') dim_align=3; end
if ~exist('dim_plot') dim_plot=[1 2 3]; end
%
opts_procrustes=struct;
opts_procrustes.if_scale=0;
if ~exist('axis_width') axis_width=1; end
if ~exist('marker_size') marker_size=14; end
%
% align data to model
%
[d_procrustes,aligned]=psg_geo_procrustes(ds{2}{dim_align},ds{1}{dim_align},opts_procrustes);
%
fig_name=opts_read_used{1}.data_fullname;
figure;
set(gcf,'NumberTitle','off');
set(gcf,'Position',[50 50 1200 800]);
set(gcf,'Name',fig_name)
%
%plot data (aligned to model) in usual fashion
% 
opts_plot=struct;
opts_plot.axis_handle=subplot(2,2,1);
psg_plotcoords(ds{1}{dim_align},dim_plot,sas{1},rayss{1},opts_plot);
%
%plot model in usual fashion
%
opts_plot.axis_handle=subplot(2,2,2);
psg_plotcoords(ds{2}{dim_align},dim_plot,sas{1},rayss{1},opts_plot);
%
%plot data without points that are off rays, and keeping only selectd rings
%
opts_plot.axis_handle=subplot(2,2,3);
opts_plot.marker_noray='';
opts_plot.if_nearest_neighbor=0;
opts_plot.noray_connect=0;
opts_plot.if_rings=1;
opts_plot.line_type_ring='-';
psg_plotcoords(ds{1}{dim_align},dim_plot,sas{1},rays_mod,opts_plot);
%
%plot modelwithout points that are off rays
%
opts_plot.axis_handle=subplot(2,2,4);
psg_plotcoords(ds{2}{dim_align},dim_plot,sas{1},rays_mod,opts_plot);
%
axes('Position',[0.01,0.02,0.01,0.01]); %for text
text(0,0,fig_name,'Interpreter','none','FontSize',8);
axis off;

% psg_vss24_task_fig_meth: create a methods figure showing a dataset as a point set, and the distances
%
%  See also:  PSG_GET_COORDSETS, PSG_PLOTCOORDS, PSG_GEO_PROCRUSTES, PSG_ALIGN_KNIT. 
%
opts_read=struct;
opts_rays=struct;
opts_qpred=struct;
if getinp('1 for interactive mode','d',[0 1])
    disp('enter two datasets: first is experimental data, second is model')
    [sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred,2);
else
    opts_read.input_type=[1 2]; %first dataset is data, second is model
    opts_read.data_fullnames={'./psg_data/bgca3pt_coords_BL_sess01_10.mat',[]};
    opts_read.setup_fullnames={'./psg_data/bgca3pt9.mat','./psg_data/bgca3pt9.mat'};
    opts_read.if_auto=1; 
    nsets=length(opts_read.input_type);
    [sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred,nsets);
end
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
%plot data as points as aligned to model
%
opts_plot.axis_handle=subplot(2,2,3);
opts_plot.if_use_rays=0;
opts_plot.marker_origin='.';
opts_plot.marker_size=marker_size;
psg_plotcoords(ds{1}{dim_align},dim_plot,sas{1},rayss{1},opts_plot);
hold on;
set(gca,'XLim',[-1 1]*max(abs(get(gca,'XLim'))));
plot3(get(gca,'XLim'),[0 0],[0 0],'k','LineWidth',axis_width);
set(gca,'YLim',[-1 1]*max(abs(get(gca,'YLim'))));
plot3([0 0],get(gca,'YLim'),[0 0],'k','LineWidth',axis_width);
set(gca,'ZLim',[-1 1]*max(abs(get(gca,'ZLim'))));
plot3([0,0],[0,0],get(gca,'ZLim'),'k','LineWidth',axis_width);
legend off;
%
%plot distances as heatmap
%
subplot(2,2,4);
dots=ds{1}{dim_align}*ds{1}{dim_align}';
npts=size(dots,1);
dsq=repmat(diag(dots),1,npts)+repmat(diag(dots)',npts,1)-2*dots;
imagesc(sqrt(dsq));
axis equal;
axis tight;
set(gca,'XTick',[1:npts]);
set(gca,'XTickLabel',sas{1}.typenames,'FontSize',7);
set(gca,'XTickLabelRotation',90.);
set(gca,'YTick',[1:npts]);
set(gca,'YTickLabel',sas{1}.typenames,'FontSize',7);
colorbar;
%
axes('Position',[0.01,0.02,0.01,0.01]); %for text
text(0,0,fig_name,'Interpreter','none','FontSize',8);
axis off;

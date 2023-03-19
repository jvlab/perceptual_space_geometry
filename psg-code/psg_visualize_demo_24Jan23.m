%psg_visualize_demo: demonstrate visualization of psg coordinates
%
% plotformats sets number of dimensions, number to plot at once
%
%  See also: PSG_READ_COORDDATA, PSG_FINDRAYS, PSG_DEFOPTS, BTC_DEFINE,
%  PSG_PLOTCOORDS, PSG_RAYFIT, PSG_RAYANGLES, PSG_SPOKES_SETUP, PSG_VISUALIZE, PSG_PLOTANGLES.
%  PSG_PROCRUSTES_DEMO, PSG_COLORS_LEGACY.
%
if ~exist('opts_plot') opts_plot=struct(); end %for psg_plotcoords
if ~exist('opts_vis') opts_vis=struct(); end %for psg_visualize
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_fit') opts_fit=struct(); end %for psg_rayfit
if ~exist('opts_ang') opts_ang=struct(); end %for psg_rayangles
if ~exist('opts_visang') opts_visang=struct(); end %for psg_plotangles
if ~exist('data_fullname') data_fullname=[]; end
if ~exist('setup_fullname') setup_fullname=[]; end
%
opts_plot=psg_colors_legacy(opts_plot);
if isfield(opts_plot,'colors')
    opts_visang.colors=opts_plot.colors;
end
if ~exist('plotformats')
    plotformats=[2 2;3 2;3 3;4 3;5 3]; %dim model, number of dims at a time
end
if ~exist('if_plotrays') if_plotrays=1; end
if ~exist('if_plotbids') if_plotbids=0; end
if ~exist('if_plotrings') if_plotrings=0; end
%
opts_read=filldefault(opts_read,'if_log',1);
[d,sa,opts_read_used]=psg_read_coorddata(data_fullname,setup_fullname,opts_read);
opts_rays.permute_raynums=opts_read_used.permute_raynums;
[rays,opts_rays_used]=psg_findrays(sa.btc_specoords,opts_rays);
%
disp(sprintf('stimulus coordinates group along %2.0f rays',rays.nrays));
origin_ptr=find(rays.whichray==0);
if length(origin_ptr)==1
    disp(sprintf('origin found at stimulus %1.0f (%s)',origin_ptr,sa.typenames{origin_ptr}));
    offset_ptr=getinp('pointer for stimulus to plot at origin or 0 for none or -1 for centroid)','d',[-1 sa.nstims],origin_ptr);
    if_origin=1;
else
    disp('origin not found');
    offset_ptr=0;
    if_origin=0;
end
opts_vis.offset_ptr=offset_ptr;
opts_vis.if_pcrot=getinp('1 to apply pca rotation','d',[0 1],0);
%
ray_counts=full(sparse(rays.whichray(rays.whichray~=0),ones(sa.nstims-if_origin,1),1,rays.nrays,1));
for iray=1:rays.nrays
    disp(sprintf('ray %2.0f: %2.0f points; endpoint: %s',iray,ray_counts(iray),sprintf('%5.2f',rays.endpt(iray,:))));
end
%
if_plotrays=getinp('1 to superimpose plots of (unidirectional, through origin) rays ','d',[0 1],if_plotrays);
if_plotbids=getinp('1 to superimpose plots of (bidirectional, not through origin) rays','d',[0 1],if_plotbids);
opts_plot.if_rings=getinp('1 to plot rings','d',[0 1],if_plotrings);
%
%for each dimension model, find best-fitting signed and unsigned rays, including the origin
%
model_dim_max=max(opts_read_used.dim_list);
d_rayfit=cell(1,model_dim_max); %coordinate structure of best-fitting rays
d_bidfit=cell(1,model_dim_max); %coordinate structure of best-fitting bidirectional rays
ray_ends=cell(1,model_dim_max); %unidirectional ray directions at max
bid_ends=cell(1,model_dim_max); %bidirectional  ray directions at max
%
angles_ray=cell(1,model_dim_max); %angle data, unidirectional
angles_bid=cell(1,model_dim_max); %angle data, bidirectional
%
opts_fit_used=cell(model_dim_max,2); %d1: model dim, d2: 1+if_bid
opts_ang_used=cell(model_dim_max,2);
%
for idimptr=1:length(opts_read_used.dim_list) %compute ray fits and angles
    idim=opts_read_used.dim_list(idimptr);
    [d_rayfit{idim},ray_ends{idim},opts_fit_used{idim,1}]=psg_rayfit(d{idim},rays,filldefault(opts_fit,'if_bid',0));
    [d_bidfit{idim},bid_ends{idim},opts_fit_used{idim,2}]=psg_rayfit(d{idim},rays,filldefault(opts_fit,'if_bid',1));
    [angles_ray{idim},opts_ang_used{idim,1}]=psg_rayangles(ray_ends{idim},sa,rays,opts_ang);
    [angles_bid{idim},opts_ang_used{idim,2}]=psg_rayangles(bid_ends{idim},sa,rays,opts_ang);
end
%
file_string=sprintf('%s %s',opts_read_used.setup_fullname,opts_read_used.data_fullname);
%
% simple plot: 2- and 3-way combinations of all axes
%
opts_vis.if_plotrays=if_plotrays;
opts_vis.if_plotbids=if_plotbids;
opts_vis.d_rayfit=d_rayfit;
opts_vis.d_bidfit=d_bidfit;
opts_vis.file_string=file_string;
[opts_vis_used,opts_plot_used]=psg_visualize(plotformats,d,sa,rays,opts_vis,opts_plot);
opts_visang_ray_used=psg_plotangles(angles_ray,sa,rays,opts_visang);
opts_visang_bid_used=psg_plotangles(angles_bid,sa,rays,opts_visang);

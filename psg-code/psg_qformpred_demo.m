%psg_qformpred_demo: demonstrate predictions from quadratic form model of thresholds
%
% This fits the threshold data two ways, but they should be equivalent:
%    see psg_qformpred_notes.doc.
%
%  d_qform:  Find the square root of the quadratic form.  Then since
%   x'Qy=[sqrt(Q)x]'[sqrt(Q)y], then sqrt(Q)x yields the coordinates that
%   whose Euclidean distances are given by the quadratic form x'Qy.  For an
%   n-dimensional model, take the first n PC's of that.
%  
% d_mds: compute the distances via the quadratic form.  Then do
%    muldimensional scaling, and take the n-dimensional model
%
% plotformats sets number of dimensions, number to plot at once
%
% 17Dec22: move logic for permute_raynums to psg_read_coorddata
%
%   See also: PSG_VISUALIZE_DEMO, PSG_READ_COORDDATA, PSG_FINDRAYS, PSG_DEFOPTS, BTC_DEFINE,
%  PSG_PLOTCOORDS, PSG_RAYFIT, PSG_RAYANGLES, PSG_SPOKES_SETUP, BTC_SOID_DEMO, PSG_QFORMPRED
%  PSG_VISUALIZE, PSG_PLOTANGLES, PSG_PROCRUSTES_DEMO.
%
if ~exist('opts_plot') opts_plot=struct(); end %for psg_plotcoords
if ~exist('opts_vis') opts_vis=struct(); end %for psg_visualize
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_fit') opts_fit=struct(); end %for psg_rayfit
if ~exist('opts_ang') opts_ang=struct(); end %for psg_rayangles
if ~exist('opts_qpred') opts_qpred=struct(); end %for psg_qformpred
%
if ~exist('data_fullname') data_fullname=[]; end
if ~exist('setup_fullname') setup_fullname=[]; end
%
if ~exist('qform_datafile') qform_datafile='../stim/btc_allraysfixedb_avg_100surrs_madj.mat';end
if ~exist('qform_modeltype') qform_modeltype=12; end %if_symm=1 (symmetrize around origin), if_axes=1 (symmetrize bc, de, tuvw); ifaug=1 (augmented coords)
%
if ~exist('plotformats')
    plotformats=[2 2;3 2;3 3;4 3;5 3]; %dim model, number of dims at a time
end
if ~exist('if_plotrays') if_plotrays=1; end
if ~exist('if_plotbids') if_plotbids=0; end
%
% read the psg setup
opts_read=filldefault(opts_read,'if_log',1);
[d,sa,opts_read_used]=psg_read_coorddata(data_fullname,setup_fullname,setfield(opts_read,'if_justsetup',1));
opts_rays.permute_raynums=opts_read_used.permute_raynums;
[rays,opts_rays_used]=psg_findrays(sa.btc_specoords,opts_rays);
%
disp(sprintf('stimulus coordinates group along %2.0f rays',rays.nrays));
if_origin=any(rays.whichray==0);
disp(sprintf('origin found: %1.0f',if_origin));
ray_counts=full(sparse(rays.whichray(rays.whichray~=0),ones(sa.nstims-if_origin,1),1,rays.nrays,1));
for iray=1:rays.nrays
    disp(sprintf('ray %2.0f: %2.0f points; endpoint: %s',iray,ray_counts(iray),sprintf('%5.2f',rays.endpt(iray,:))));
end
%
if_plotrays=getinp('1 to superimpose plots of (unidirectional, through origin) rays ','d',[0 1],if_plotrays);
if_plotbids=getinp('1 to superimpose plots of (bidirectional, not through origin) rays','d',[0 1],if_plotbids);
%
%get quadratic form
%
nbtc=size(sa.btc_augcoords,2);
qform_source_type=getinp(' quadratic form choice: 1->use threshold data file, 2->use identity','d',[1 2],1);
switch qform_source_type
    case 1
        btc_thresh_data=getfield(load(qform_datafile),'r');
        q=btc_thresh_data{qform_modeltype}.results.qfit;
        q_label=btc_thresh_data{qform_modeltype}.setup.label;
        disp(sprintf(' model %2.0f (%s) loaded from %s',...
            qform_modeltype,q_label,qform_datafile));
        qform_source=qform_datafile;
    case 2
        q=eye(nbtc);
        qform_source=sprintf('[%2.0f x %2.0f identity]',nbtc,nbtc);
        q_label=' ';
end
pred_string_pref=' ';
if_aug_spe=getinp('1 to use augmented coords, 2 to use spec coords','d',[1 2],1);
switch if_aug_spe
    case 1
        btc_coords=sa.btc_augcoords;
        aug_spe_string='aug coords';
    case 2
        btc_coords=sa.btc_specoords;
        btc_coords(isnan(btc_coords))=0;
        aug_spe_string='spec coords';
end  
%
%fit quadratic form two ways
%
[d_qform,d_mds,opts_qpred_used]=psg_qformpred(q,btc_coords,rays,opts_qpred);
%
disp('max deviations between coordinates of the two methods (sign flipped if needed)');
for idim=1:nbtc
    disp(sprintf(' %8.6f',min([max(abs(d_qform{idim}-d_mds{idim}),[],1);max(abs(d_qform{idim}+d_mds{idim}),[],1)],[],1)));
end
%
%for each dimension model, find best-fitting signed and unsigned rays, including the origin
%
dim_list=[1:nbtc];
opts_read_used.dim_list=dim_list;
model_dim_max=max(opts_read_used.dim_list);
%
opts_vis.if_plotrays=if_plotrays;
opts_vis.if_plotbids=if_plotbids;
for ipred=1:2
    if (ipred==1)
        d_pred=d_qform;
        pred_string_meth=sprintf('qform if_centroid_pca=%1.0f',opts_qpred_used.if_pca_centroid);
    else
        d_pred=d_mds;
        pred_string_meth='mds';
    end
    pred_string=cat(2,pred_string_pref,' method:',pred_string_meth);
    d_rayfit=cell(1,model_dim_max); %coordinate structure of best-fitting rays
    d_bidfit=cell(1,model_dim_max); %coordinate structure of best-fitting bidirectional rays
    ray_ends=cell(1,model_dim_max); %ray directions at max
    bid_ends=cell(1,model_dim_max); %ray directions at max
    %
    angles_ray=cell(1,model_dim_max); %angle data, unidirectional
    angles_bid=cell(1,model_dim_max); %angle data, bidirectional
    %
    opts_fit_used=cell(model_dim_max,2); %d1: model dim, d2: 1+if_bid
    opts_ang_used=cell(model_dim_max,2);
    %
    for idimptr=1:length(opts_read_used.dim_list) %compute ray fits and angles
        idim=opts_read_used.dim_list(idimptr);
        [d_rayfit{idim},ray_ends{idim},opts_fit_used{idim,1}]=psg_rayfit(d_pred{idim},rays,filldefault(opts_fit,'if_bid',0));
        [d_bidfit{idim},bid_ends{idim},opts_fit_used{idim,2}]=psg_rayfit(d_pred{idim},rays,filldefault(opts_fit,'if_bid',1));
        [angles_ray{idim},opts_ang_used{idim,1}]=psg_rayangles(ray_ends{idim},sa,rays,opts_ang);
        [angles_bid{idim},opts_ang_used{idim,2}]=psg_rayangles(bid_ends{idim},sa,rays,opts_ang);
    end
    if (ipred==1)
        angles_ray_qform=angles_ray;
        angles_bid_qform=angles_bid;
    else
        angles_ray_mds=angles_ray;
        angles_bid_mds=angles_bid;
    end
    %
    opts_vis.d_rayfit=d_rayfit;
    opts_vis.d_bidfit=d_bidfit;
    opts_vis.file_string=sprintf('%s %s %s %s method: %s',opts_read_used.setup_fullname,aug_spe_string,qform_source,q_label,pred_string_meth);
    opts_vis.vis_string_format=cat(2,'raw %1.0fd fit, ',pred_string);
    if (ipred==1)
        [opts_vis_used_qform,opts_plot_used_qform]=psg_visualize(plotformats,d_qform,sa,rays,opts_vis,opts_plot);
        opts_visang_ray_used_qform=psg_plotangles(angles_ray_qform,sa,rays,opts_vis);
        opts_visang_bid_used_qform=psg_plotangles(angles_bid_qform,sa,rays,opts_vis);
    else
        [opts_vis_used_mds,opts_plot_used_mds]=psg_visualize(plotformats,d_mds,sa,rays,opts_vis,opts_plot);
        opts_visang_ray_used_mds=psg_plotangles(angles_ray_mds,sa,rays,opts_vis);
        opts_visang_bid_used_mds=psg_plotangles(angles_bid_mds,sa,rays,opts_vis);
    end
end %ipred (prediction type)

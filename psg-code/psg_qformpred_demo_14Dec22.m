%psg_qformpred_demo: demonstrate predictions from quadratic form model of thresholds
%
% This fits the threshold data two ways:
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
% To do: this code, and psg_visualize_demo, would benefit from a common
% call to a plotting routine, modularizing the two fitting procedures
% with options including whether to center at rand or centroid, what to do about negative 
% eigenvalues, and having a pca routine for model coordinates
% This will also enable passing a unit qform for an information-theoretic
% (ideal observer) distance Also should be runnable on augmented or specified coordinates.
%
%  See also: PSG_VISUALIZE_DEMO, PSG_READ_COORDDATA, PSG_FINDRAYS, PSG_DEFOPTS, BTC_DEFINE,
%  PSG_PLOTCOORDS, PSG_RAYFIT, PSG_RAYANGLES, PSG_SPOKES_SETUP, NICESUBP
%  BTC_SOID_DEMO.
%
if ~exist('opts_vis') opts_vis=struct(); end %for psg_plotcoords
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_fit') opts_fit=struct(); end %for psg_rayfit
if ~exist('opts_ang') opts_ang=struct(); end %for psg_rayangles
if ~exist('data_fullname') data_fullname=[]; end
if ~exist('setup_fullname') setup_fullname=[]; end
if ~exist('lim_margin') lim_margin=1.1; end %margin for axis limits
%
if ~exist('qform_datafile') qform_datafile='../stim/btc_allraysfixedb_avg_100surrs_madj.mat';end
if ~exist('qform_modeltype') qform_modeltype=12; end %if_symm=1 (symmetrize around origin), if_axes=1 (symmetrize bc, de, tuvw); ifaug=1 (augmented coords)
%
permutes.bgca=[2 1 3 4]; % permute ray numbers to the order gbca 
%
if ~exist('plotformats')
    plotformats=[2 2;3 2;3 3;4 3;5 3]; %dim model, number of dims at a time
end
% read the psg setup
opts_read=filldefault(opts_read,'if_log',1);
[d,sa,opts_read_used]=psg_read_coorddata(data_fullname,setup_fullname,setfield(opts_read,'if_justsetup',1));
%
if strfind(opts_read_used.setup_fullname,'bgca')
    disp('suggested row permutation:')
    disp(permutes.bgca);
    if getinp('1 if ok','d',[0 1],1)
        opts_rays.permute_raynums=permutes.bgca;
    end
end
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
%get quadratic form
if ~exist('qform_datafile') qform_datafile='../stim/btc_allraysfixedb_avg_100surrs_madj.mat';end
btc_thresh_data=getfield(load(qform_datafile),'r');
q=btc_thresh_data{qform_modeltype}.results.qfit;
q_label=btc_thresh_data{qform_modeltype}.setup.label;
disp(sprintf(' model %2.0f (%s) loaded from %s',...
    qform_modeltype,q_label,qform_datafile));
%diagonalize
[qeivecs,qeivals]=eig(q);
[qeivals_sorted,inds]=sort(diag(qeivals),'descend');
qeivecs_sorted=qeivecs(:,inds);
nbtc=length(qeivals_sorted);
disp(sprintf('%3.0f eigenvalues, largest is %5.2f, smallest is %5.2f, %3.0f are less than 0 and set to 0',...
    nbtc,qeivals_sorted(1),qeivals_sorted(end),sum(qeivals_sorted<0)));
qeivals_sorted(qeivals_sorted<0)=0;
disp('eigenvectors as columns, sorted from largest to smallest eigenvalues')
disp(qeivecs_sorted)
%
%compute coordinates (square root of the quadratic form, computed "by hand" from eigenvectors
%and square roots of eigenvalues)
%
qeivecs_scaled=qeivecs_sorted.*repmat(sqrt(qeivals_sorted(:))',nbtc,1);
qform_predcoords=sa.btc_augcoords*qeivecs_scaled;
if (if_origin) %if origin is found, use it; otherwise, use centroid
    qform_origin=qform_predcoords(find(rays.whichray==0),:);
    disp('origin for qform coords taken from random stimulus.');
else
    qform_origin=mean(qform_predcoords,1);
    disp('origin for qform coords taken from centroid');
end
qform_predcoords=qform_predcoords-repmat(qform_origin,sa.nstims,1); %move origin to zero
%
opts_read_used.dim_list=[1:nbtc];
%
%First method:
%For each dimension, use the highest principal components of the full model obtained from the square root
%of the quadratic form. Note that this PCA step is necessary, since, depending on the
%texture params to be modeled, i.e., btc_augcoords, only the "higher" dimensions of the MDS fit 
%may be relevant
%
[qu,qs,qv]=svd(qform_predcoords,0); %qform_predcoords = qu*qs*qv';
d_qform=cell(1,nbtc);
var_tot=sum(qform_predcoords(:).^2);
for idim=1:nbtc
    d_qform{idim}=qu(:,[1:idim])*qs([1:idim],[1:idim]); %coordinates in pca space
    d_recon=d_qform{idim}*(qv(:,[1:idim]))';
    coord_maxdiff=max(max(abs(d_recon-qform_predcoords)));
    var_ex=sum(d_qform{idim}(:).^2);
    disp(sprintf(' %2.0f-dim model, variance explained =%7.3f of %7.3f (frac unex: %8.5f), max coord diff to %2.0f-dim MDS is %6.3f',...
        idim,var_ex,var_tot,1-var_ex/var_tot,nbtc,coord_maxdiff));
end
%
%Second method: do MDS on the predicted coordinates and then peel off
%successive dimensions
%
d_mds=cell(1,nbtc);
qf_dsq=zeros(sa.nstims,sa.nstims);
for istim=1:sa.nstims
    for jstim=1:sa.nstims
        qf_dsq(istim,jstim)=sum((qform_predcoords(istim,:)-qform_predcoords(jstim,:)).^2);
    end
end
qf_centered=repmat(mean(qf_dsq),sa.nstims,1);
qf_mds=qf_dsq-qf_centered-qf_centered'+mean(qf_dsq(:));
[qf_vec,qf_val]=eig(-qf_mds/2);
qform_mdscoords=real(qf_vec)*diag(sqrt(max(0,real(diag(qf_val)))));
if (if_origin) %if origin is found, use it; otherwise, use centroid
    qform_origin=qform_mdscoords(find(rays.whichray==0),:);
    disp('origin for mds coords taken from random stimulus.');
else
    qform_origin=mean(qform_mdscoords,1);
    disp('origin for mds coords taken from centroid');
end
qform_mdscoords=qform_mdscoords-repmat(qform_origin,sa.nstims,1);
for idim=1:nbtc
    d_mds{idim}=qform_mdscoords(:,[1:idim]);
end
%
if_plotrays=getinp('1 to superimpose plots of (unidirectional) rays','d',[0 1]);
if_plotbids=getinp('1 to superimpose plots of (bidirectional)  rays','d',[0 1]);
%
%for each dimension model, find best-fitting signed and unsigned rays, including the origin
%
model_dim_max=max(opts_read_used.dim_list);
for ipred=1:2
    if (ipred==1)
        d_pred=d_qform;
        pred_string='qform';
    else
        d_pred=d_mds;
        pred_string='mds';
    end
    d_rayfit=cell(1,model_dim_max); %coordinate structure of best-fitting rays
    d_bidfit=cell(1,model_dim_max); %coordinate structure of best-fitting bidirectional rays
    ray_ends=cell(1,model_dim_max); %ray directions at max
    bid_ends=cell(1,model_dim_max); %ray directions at max
    %
    angles=cell(1,model_dim_max); %structure of angle data
    %
    opts_fit_used=cell(model_dim_max,2); %d1: model dim, d2: 1+if_bid
    opts_ang_used=cell(model_dim_max,1);
    %
    for idimptr=1:length(opts_read_used.dim_list) %compute ray fits and angles
        idim=opts_read_used.dim_list(idimptr);
        [d_rayfit{idim},ray_ends{idim},opts_fit_used{idim,1}]=psg_rayfit(d_pred{idim},rays,filldefault(opts_fit,'if_bid',0));
        [d_bidfit{idim},bid_ends{idim},opts_fit_used{idim,2}]=psg_rayfit(d_pred{idim},rays,filldefault(opts_fit,'if_bid',1));
        [angles{idim},opts_ang_used{idim}]=psg_rayangles(ray_ends{idim},sa,rays,opts_ang);
    end
    if (ipred==1)
        angles_qform=angles;
    else
        angles_mds=angles;
    end
    %
    file_string=sprintf('%s augcoords %s %s',opts_read_used.setup_fullname,qform_datafile,q_label);
    %
    % simple plot: 2- and 3-way combinations of all axes
    %
    opts_vis_raw=cell(1,size(plotformats,1));
    for iplot=1:size(plotformats,1)
        model_dim=plotformats(iplot,1);
        dims_together=plotformats(iplot,2);
        if ismember(model_dim,opts_read_used.dim_list)
            vis_string=sprintf('raw %1.0fd fit, prediction type: %s',model_dim,pred_string);
            figure;
            set(gcf,'Position',[50 50 1200 800]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',vis_string);
            dim_selects=nchoosek([1:model_dim],dims_together);
            ncombs=size(dim_selects,1);
            [nr,nc]=nicesubp(ncombs,0.7);
            opts_vis.lims=lim_margin*[-1 1]*max(abs(d_pred{model_dim}(:)));
            for icomb=1:ncombs
                ha=subplot(nr,nc,icomb);
                opts_vis_raw{iplot}=psg_plotcoords(d_pred{model_dim},dim_selects(icomb,:),sa,rays,setfield(opts_vis,'axis_handle',ha));
                ha=opts_vis_raw{iplot}.axis_handle;
                %
                if if_plotrays
                    psg_plotcoords(d_rayfit{model_dim},dim_selects(icomb,:),sa,rays,setfields(opts_vis,{'axis_handle','line_type','if_just_data'},{ha,':',1}));
                end
                if if_plotbids
                    psg_plotcoords(d_bidfit{model_dim},dim_selects(icomb,:),sa,rays,setfields(opts_vis,{'axis_handle','line_type','if_just_data','if_origin_on_rays'},{ha,'--',1,0}));
                end
            end
            axes('Position',[0.01,0.04,0.01,0.01]); %for text
            text(0,0,vis_string,'Interpreter','none');
            axis off;
            axes('Position',[0.01,0.02,0.01,0.01]); %for text
            text(0,0,file_string,'Interpreter','none');
            axis off;
        end %ismember
    end %iplot
 end %itype

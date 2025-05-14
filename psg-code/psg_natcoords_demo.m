%psg_natcoords_demo: find natural axes in a set of psg coordinates
%
% based on psg_visualize_demo, which uses psg_rayfit to calculate rays from origin
% this determines tangents to each stimulus-coordinate trajectory at the origin, 
% and related quantitites, which serve as natural coordinates,
% and examines a transformation into those natural coordinates
%
% to do:
% extract a transformation; find the major axes of the transformation and
% express their angles w.r.t. the natural coords; find expansion in the
% direction of the natural coordinates;  for transformations with cutpoints,
% then relate cutpoints to directions of the natural coordinates
%
%  See also: PSG_GET_COORDSETS, PSG_READ_COORDDATA, PSG_FINDRAYS, PSG_DEFOPTS, BTC_DEFINE,
%  PSG_PLOTCOORDS, PSG_RAYFIT, PSG_SPOKES_SETUP, PSG_VISUALIZE,
%  PSG_VISUALIZE_DEMO, PSG_NATCOORDS,
%  PSG_GET_GEOTRANSFORMS, PSG_GEOMODELS_DEFINE.
%
if ~exist('opts_plot') opts_plot=struct(); end %for psg_plotcoords
if ~exist('opts_vis') opts_vis=struct(); end %for psg_visualize
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_fit') opts_fit=struct(); end %for psg_rayfit
if ~exist('opts_visang')  opts_visang=struct(); end %for psg_plotangles
if ~exist('opts_vismult') opts_vismult=struct(); end %for psg_plotangles
if ~exist('opts_qpred') opts_qpred=struct(); end %for psg_qformpred
if ~exist('data_fullname') data_fullname=[]; end
if ~exist('setup_fullname') setup_fullname=[]; end
%
if ~exist('opts_natc') opts_natc=struct; end %for psg_natcoords
if ~exist('opts_geot') opts_geot=struct(); end % for psg_get_geotransforms
%
if ~exist('magnif_tol') magnif_tol=10^-5; end %tolerance for checking mangifications against eigenvalues
%
opts_geot=filldefault(opts_geot,'if_uigetfile',0); %use a console prompt
%
model_types_def=psg_geomodels_define();
%
opts_read=filldefault(opts_read,'if_spray',0); %default to not define rays by single points
if (opts_read.if_spray==0) %if ray mode is never used, color is irrelevant
    opts_plot=psg_colors_legacy(opts_plot);
    if isfield(opts_plot,'colors')
        opts_visang.colors=opts_plot.colors;
        opts_vismult.colors=opts_plot.colors;
    end
end
%
if ~exist('plotformats')
    plotformats=[2 2;3 2;3 3;4 3;5 3]; %dim model, number of dims at a time
end
if ~exist('if_plotrays') if_plotrays=1; end
if ~exist('if_plotbids') if_plotbids=0; end
if ~exist('if_plotrings') if_plotrings=0; end
if ~exist('if_nearest_neighbor') if_nearest_neighbor=-1; end
%
if_spray=getinp('1 to define rays by single points, suppress ray angle calculations, and suppress ray angle plots','d',[0 1],opts_read.if_spray);
if (if_spray)
    opts_rays.ray_minpts=1;
    if_plotrays=0;
end
%
nsets=1;
opts_read=filldefault(opts_read,'if_log',1);
[sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred,nsets);
d=ds{1};
sa=sas{1};
rays=rayss{1};
opts_read_used=opts_read_used{1};
opts_rays_used=opts_rays_used{1};
opts_rays_used.permute_raynums=opts_read_used.permute_raynums;
%
disp(sprintf('stimulus coordinates group along %2.0f rays',rays.nrays));
origin_ptr=find(rays.whichray==0);
if length(origin_ptr)==1
    disp(sprintf('origin found at stimulus %1.0f (%s)',origin_ptr,sa.typenames{origin_ptr}));
    offset_ptr=getinp('pointer for stimulus to plot at origin or 0 for none or -1 for centroid)','d',[-1 sa.nstims],origin_ptr);
else
    disp('origin not found');
    offset_ptr=0;
    opts_fit.if_origin=0;
end
opts_vis.offset_ptr=offset_ptr;
opts_vis.if_pcrot=getinp('1 to apply pca rotation when plotting','d',[0 1],0);
%
ray_counts=full(sparse(rays.whichray(rays.whichray>0),ones(sum(rays.whichray>0),1),1,rays.nrays,1));
for iray=1:rays.nrays
    disp(sprintf('ray %2.0f: %2.0f points; endpoint: %s',iray,ray_counts(iray),sprintf('%5.2f',rays.endpt(iray,:))));
end
%
if_plotrays=getinp('1 to superimpose plots of (unidirectional, through origin) rays ','d',[0 1],if_plotrays);
if_plotbids=getinp('1 to superimpose plots of (bidirectional, not through origin) rays','d',[0 1],if_plotbids);
opts_plot.if_rings=getinp('1 to plot rings','d',[0 1],if_plotrings);
opts_plot.if_nearest_neighbor=getinp('1 to connect nearest neighbors, 0 not, -1 if unassigned points','d',[-1 1],if_nearest_neighbor);
%
%for each dimension model, find best-fitting signed and unsigned rays, including the origin
%
dim_list=sets{1}.dim_list;
model_dim_max=max(dim_list);
d_rayfit=cell(1,model_dim_max); %coordinate structure of best-fitting rays
d_bidfit=cell(1,model_dim_max); %coordinate structure of best-fitting bidirectional raysopts_mult_used
ray_ends=cell(1,model_dim_max); %unidirectional ray directions at max
bid_ends=cell(1,model_dim_max); %bidirectional  ray directions at max
%
angles_ray=cell(1,model_dim_max); %angle data, unidirectional
angles_bid=cell(1,model_dim_max); %angle data, bidirectional
%
opts_fit_used=cell(model_dim_max,2); %d1: model dim, d2: 1+if_bid
%
opts_natc_used=cell(model_dim_max,1);
natcoords=cell(1,model_dim_max);
%
file_string=sprintf('%s %s',opts_read_used.setup_fullname,opts_read_used.data_fullname);
%
if ~if_spray
    for idimptr=1:length(dim_list)
        idim=dim_list(idimptr);
        %compute ray fits
        [d_rayfit{idim},ray_ends{idim},opts_fit_used{idim,1}]=psg_rayfit(d{idim},rays,filldefault(opts_fit,'if_bid',0));
        [d_bidfit{idim},bid_ends{idim},opts_fit_used{idim,2}]=psg_rayfit(d{idim},rays,filldefault(opts_fit,'if_bid',1));
        %compute tangents
        [natcoords{idim},opts_natc_used{idim}]=psg_natcoords(d{idim},sa,rays,opts_natc);
    end
end
%
if ~if_spray
    %get a geometrical transform
    [transforms_avail,dims_avail,geot_desc,opts_geot_used]=psg_get_geotransforms(opts_geot); %get a geometric transform from a file
    disp(sprintf('analyzing transformations in %s',opts_geot_used.fullname));
    disp(sprintf('%s used for domain natural coordinates',file_string));
    disp(geot_desc)
    for idimptr=1:min(length(dim_list),max(dims_avail))
        idim=dim_list(idimptr);
        if ~isempty(transforms_avail{idim}) & ~isempty(natcoords{idim})
            disp(' ');
            npieces=size(transforms_avail{idim}.T,3);
            for ipiece=1:npieces
                tr=transforms_avail{idim}.b*transforms_avail{idim}.T(:,:,ipiece); %include overall scale factor
                for natc=1:length(natcoords{idim}.avail)
                    natc_name=natcoords{idim}.avail{natc};
                    natc_vecs=natcoords{idim}.(natc_name);
                    natc_vecs_length=sqrt(sum(natc_vecs.^2,2)); %to turn dot-prdouct into cosines
                    if ~any(isnan(natc_vecs))
                        disp(sprintf('natural coordinate analysis using %12s, transformation dimension %2.0f, piece %2.0f',natc_name,idim,ipiece));
                        ray_string='';
                        for iray=1:rays.nrays
                            ray_string=cat(2,ray_string,sprintf('%12s',natcoords{idim}.labels{iray}),' ');
                        end
                        disp(sprintf('                                 length  magnif  cosines: %s',ray_string));
                        [eivecs,eivals]=eig(tr*tr');
                        eivals=real(diag(eivals)); %A is self-adjoint
                        [eivals,sort_inds]=sort(eivals,'descend'); %obtain eigenvalues in descending order
                        eivecs=real(eivecs(:,sort_inds));                       
                        for iv=1:rays.nrays+idim %look at natural coords and eigenvectors
                            if (iv<=rays.nrays)
                                vec_lab=sprintf('nat coord on ray %2.0f (%s)',iv,natcoords{idim}.labels{iv});
                                vec=natc_vecs(iv,:);
                                iev=0;
                            else
                                %
                                % Here we compare the major axes of the transformation to the natural coordinates. Logic is same as psg_majaxes, adjusted-dimension case.
                                % if tr is the transform matrix: then let [eivecs,eivals ]=eig(tr*tr').  Then, for each column r in eivecs, (r')tr*tr' is a multiple of r',
                                % and the eigenvalues d are the squares of the expansions (magnif_check). 
                                % 
                                % The eigenvectors are then compared sith the natural coordinates via cosines (dot produts normalized by the vector lengths)
                                %
                                iev=iv-rays.nrays;
                                vec_lab=sprintf('  transformation eigenvec %2.0f',iev);
                                vec=eivecs(:,iev)'; %eigenvectors are returned in columns but we need them in rows
                                magnif_check=sqrt(eivals(iev));
                            end
                            %
                            %  natural coord vectors are post-multiplied by tr to see how % much they expand
                            %
                            vec_length=sqrt(sum(vec.^2));
                            vec_transformed=vec*tr;
                            vec_length_transformed=sqrt(sum(vec_transformed.^2));
                            vec_magnif=vec_length_transformed/vec_length;
                            dots=vec*natc_vecs';
                            cosines=abs(dots./vec_length./natc_vecs_length');
                            cosine_string=sprintf('      %7.3f',cosines);
                            disp(sprintf('%s    %7.3f %7.3f           %s',vec_lab,vec_length,vec_magnif,cosine_string));
                            if iev>0
                                if abs(magnif_check-vec_magnif)>magnif_tol
                                    disp(sprintf(' disagreement in magnfication factor, sqrt(eiv) is %7.3f',magnif_check))
                                end
                            end %eiv check
                        end %iv
                    end %all present
                end %natc
            end %ipiece
        end %transformation present
    end %dim
end
% simple plot: 2- and 3-way combinations of all axes
%
opts_vis.if_plotrays=if_plotrays;
opts_vis.if_plotbids=if_plotbids;
opts_vis.d_rayfit=d_rayfit;
opts_vis.d_bidfit=d_bidfit;
opts_vis.file_string=file_string;
[opts_vis_used,opts_plot_used]=psg_visualize(plotformats,d,sa,rays,opts_vis,opts_plot);

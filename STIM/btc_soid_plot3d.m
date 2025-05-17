function [hplots,edirs_aug,ou]=btc_soid_plot3d(edirs,dirlist,opts)
% [hplots,edirs_aug,ou]=btc_soid_plot(edirs,dirlist,opts) plots directions, conditions,
% and data for a btc eperiment, as an ellipsoid; standard dictionary assumed
%
% for each coordinate set in dirlist, the isodiscrimination contours from the planes
% contained in this subspace are shown.  These "ellipses" are modified by adding 
% the nonzero coordinates from augcoords.  For example, when plotting the 'tg' data
% into the 'atg' space, a=g^4.  Most of the code is general, except for the
% final plotting (where dirlist must be a set of triplets)
%
% edirs:  a structure with fields for each plane ("xy")
%   indicating the experimental directions
%   key fields are edirs.yx.plane (a 2-char string)
%   and edirs.yx.allvecs_concat
% dirlist:  a cell array of sets of triplets
% opts: plot options
%    opts.ncircle: circle quality, default to 48
%    opts.radii: circle radii, defaults to [1:4]/4
%    opts.fontsize: font size, defaults to 8
%    opts.lims: limits, defaults to 1
%    opts.datafield: strvcat of data field(s) to plot; can be a strvcat to plot more than one variable
%       this must have two columns, first column is abscissa, second is ordinate
%    opts.marker: symbols to use, e.g., 'k.'; can be a strvcat to plot more than one variable
%    opts.aug_opts: options used for btc_augcoords, defaults to []
%    opts.dperm: the dirlist coords that correspond to x, y, and z; [3 2 1] so z is plotted up
%    opts.show_planes:  1 to show the coordinate planes in the plot
%         (default).  A fraction <1 sjhows the planes with an "alpha" equal to
%    show_planes (recommend 0.1)
%    opts.view:  arguments for 3-d view, defaults to [-50 15]
%    opts.use_sym:  1 to use symmetries to find coordinate planes, 0 not to do so (defaults to 1)
%    opts.verbose: 1 to provide a log of substitutions (defaults to 0)
%    opts.tstring: suffix for title string
%
% hplots: pointers to the individual axes of the plots
% ou: options used
% edirs_aug:  edirs, with new fields, one for each member of datafield,
%    listing all 10 coordinate values at the threshold
%
%  See also:  BTC_SOID_TEST, BTC_SOID_PLOT, BTC_AUGCOORD, BTC_ROTCODE, BTC_HVI.
%
if (nargin<=2) opts=[]; end
opts=filldefault(opts,'fontsize',8);
opts=filldefault(opts,'ncircle',48); %0 for no circle
opts=filldefault(opts,'radii',[1:4]/4); %empty for no circle
opts=filldefault(opts,'lims',1);
opts=filldefault(opts,'marker','k.');
opts=filldefault(opts,'cyclic',0);
opts=filldefault(opts,'aug_opts',setfield([],'ifstd',1));
opts=filldefault(opts,'dperm',[1:length(dirlist{1})]);
opts=filldefault(opts,'show_planes',1);
opts=filldefault(opts,'view',[-50 15]);
opts=filldefault(opts,'use_sym',1);
opts=filldefault(opts,'verbose',0);
opts=filldefault(opts,'tstring',' ');
%
opts.aug_opts.ifstd=1;
opts.aug_opts.nocheck=1;
dict_std=btc_define;
opts.dict_std=dict_std;
%
nsets=length(dirlist);
[nr,nc]=nicesubp(nsets,0.7);
%
ou=opts;
%
planes=char(fieldnames(edirs));
nplanes=size(planes,1);
%
nvars=size(opts.datafield,1);
%
for iplane=1:nplanes
    edir=getfield(edirs,planes(iplane,:));
    for ivar=1:nvars
        vecs_inplane=getfield(edir,deblank(opts.datafield(ivar,:)));
        datafield_aug=[];
        for icond=1:edir.ndirs
            spec=[];
            for ix=1:2
                %the "3-ix" is so that the SECOND coordinate of vecs_inplane 
                %is the FIRST btc coordinate (in planes)
                spec=setfield(spec,planes(iplane,ix),vecs_inplane(icond,3-ix));
            end
            augcoords=btc_augcoords(spec,opts.dict_std,opts.aug_opts);
            avec=augcoords.method{1}.vec;
            datafield_aug(icond,:)=avec;
        end %did each condition
        edir=setfield(edir,cat(2,deblank(opts.datafield(ivar,:)),'_allcoords'),datafield_aug);
    end %ivar
    edirs=setfield(edirs,planes(iplane,:),edir);
end %iplane
edirs_aug=edirs;
%
dperm=opts.dperm;
for iset=1:nsets
    hplots(iset)=subplot(nr,nc,iset,'FontSize',opts.fontsize);
    dirs=dirlist{iset};
    ndims=length(dirs);
    %consider each plot pair
    pairs=nchoosek([1:ndims],2);
    %plot locus of thresholds from the ellipse in each plane within specified coordinate set
    for ipair=1:size(pairs,1)
        req_plane=[dirs(pairs(ipair,1)),dirs(pairs(ipair,2))];
        use_plane=req_plane;
        use_dirs=dirs;
        have_plane=0;
        %look for a rotation or a hv flip that matches the required planes with one of the fields of edir
        for ihv=0:1*opts.use_sym
            for irot=0:3*opts.use_sym
                if (have_plane==0)
                    if (isfield(edirs,use_plane)) %is the field present?
                        have_plane=1;
                    elseif (isfield(edirs,fliplr(use_plane)))
                        have_plane=1;
                        use_plane=fliplr(use_plane);
                    else
                        use_plane=char(btc_rotcode(use_plane));
                        use_dirs=char(btc_rotcode(use_dirs));
                    end
                end %have_plane
            end %irot
            if (have_plane==0) %try an hv flip
                use_plane=char(btc_hvi(use_plane));
                use_dirs=char(btc_hvi(use_dirs));
            end
        end %ihv
        if ~(have_plane==0)
            if (opts.verbose==1)
                disp(sprintf('  using plane %s for data for %s to augment and project into %s',use_plane,req_plane,dirs));
            end
            for ivar=1:nvars
                inplane_data=getfield(edirs.(use_plane),cat(2,deblank(opts.datafield(ivar,:)),'_allcoords'));
                plot_data=[];
                for idim=1:ndims
                    plot_data(:,idim)=inplane_data(:,find(dict_std.codel==use_dirs(dperm(idim))));
                end
                %assume dimension =3; up to here, code has been general for
                %number of dimensions
                xyz=plot_data(:,1:3);
                marker=deblank(opts.marker(min(ivar,size(opts.marker,1)),:));
                cyclic=opts.cyclic(min(ivar,length(opts.cyclic)));
                indices=[1:size(xyz,1)]';
                if (cyclic==1)
                    indices=[indices;1];
                end
                plot3(xyz(indices,1),xyz(indices,2),xyz(indices,3),marker);
                hold on;
            end
        else %the plane is not present
            disp(sprintf(' no plane found with data for %s to augment and project into %s',req_plane,dirs));
        end % if have_plane
    end
    plot3([-1 1]*opts.lims,[0 0],[0,0],'k');
    plot3([0,0],[-1 1]*opts.lims,[0,0],'k');
    plot3([0,0],[0 0],[-1 1]*opts.lims,'k');
    xlabel(dirs(dperm(1)),'FontSize',opts.fontsize);
    ylabel(dirs(dperm(2)),'FontSize',opts.fontsize);
    zlabel(dirs(dperm(3)),'FontSize',opts.fontsize);
    title(deblank(cat(2,dirs,' ',opts.tstring)));
    set(gca,'XLim',[-1 1]*opts.lims);
    set(gca,'YLim',[-1 1]*opts.lims);
    set(gca,'ZLim',[-1 1]*opts.lims);
    set(gca,'XTick',[-1 1]*opts.lims);
    set(gca,'YTick',[-1 1]*opts.lims);
    set(gca,'ZTick',[-1 1]*opts.lims);
    tdisp=[0.9 -0.1 0.1]*opts.lims;
    text(tdisp(1),tdisp(2),tdisp(3),dirs(dperm(1)),'FontSize',opts.fontsize);
    text(tdisp(3),tdisp(1),tdisp(2),dirs(dperm(2)),'FontSize',opts.fontsize);
    text(tdisp(2),tdisp(3),tdisp(1),dirs(dperm(3)),'FontSize',opts.fontsize);
    %
    if (opts.show_planes>0)
        if (opts.show_planes==1)
            plot3([-1 1 1 -1 -1]*opts.lims,[1 1 -1 -1 1]*opts.lims,[0 0 0 0 0],'k');
            plot3([0 0 0 0 0],[-1 1 1 -1 -1]*opts.lims,[1 1 -1 -1 1]*opts.lims,'k');
            plot3([1 1 -1 -1 1]*opts.lims,[0 0 0 0 0],[-1 1 1 -1 -1]*opts.lims,'k');
        else
            hpatch=patch([-1 1 1 -1 -1]*opts.lims,[1 1 -1 -1 1]*opts.lims,[0 0 0 0 0],'k');
            set(hpatch,'FaceAlpha',opts.show_planes);
            hpatch=patch([0 0 0 0 0],[-1 1 1 -1 -1]*opts.lims,[1 1 -1 -1 1]*opts.lims,'k');
            set(hpatch,'FaceAlpha',opts.show_planes);
            hpatch=patch([1 1 -1 -1 1]*opts.lims,[0 0 0 0 0],[-1 1 1 -1 -1]*opts.lims,'k');
            set(hpatch,'FaceAlpha',opts.show_planes);
        end
    end
    %
    axis vis3d;
    set(gca,'View',opts.view);
end %the set
return

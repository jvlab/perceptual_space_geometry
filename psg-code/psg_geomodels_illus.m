%psg_geomodels_illus: demonstrate various geometric transforms with simple
%starting points, using standard visualization tools
%
%  See also:  PSG_MODELS_DEFINE, PSG_PWAFFINE_APPLY, PERSP_APPLY, 
%    PSG_FINDRAYS, PSG_VISUALIZE, PSG_PLOTCOORDS.
%
% multiple setups
% affine, projective, piecewise
% add noise
%
if ~exist('dim_max') dim_max=3; end %max dimension for layouts
%
if ~exist('dim_plot') dim_plot=min(dim_max,3); end %for psg_plotcoords
if ~exist('dim_select') dim_select=1:dim_plot; end %for psg_plotcoords
if ~exist('plotformats') plotformats=[dim_plot dim_plot]; end %plot format for psg_visualize: plot 3 dims, all together
%
if ~exist('opts_rays') opts_rays=struct; end %for psg_findrays
if ~exist('opts_plot') opts_plot=struct; end %for psg_visualize
if ~exist('opts_vis') opts_vis=struct; end %for psg_visualize
if ~exist('opts_mult') opts_mult=struct; end %for psg_visualize
%
if ~exist('transform_view') transform_view=[-50 21]; end
%
if ~exist('layouts')
    layouts=cell(0);
    nxy=[5 5];
    layouts{1}.name=sprintf('%1.0f x %1.0f grid',nxy);
    [xg,yg]=meshgrid([1:nxy(1)],[1:nxy(2)]);
    layouts{1}.coords=[xg(:),yg(:)];
    %
    nxy=[5 7];
    layouts{2}.name=sprintf('%1.0f x %1.0f grid',nxy);
    [xg,yg]=meshgrid([1:nxy(1)],[1:nxy(2)]);
    layouts{2}.coords=[xg(:),yg(:)];
    %
    nrays=8;
    neachray=3;
    layouts{3}.name=sprintf('%1.0f rays, %1.0f points each',nrays,neachray);
    angs=2*pi*floor([0:nrays*neachray-1]/neachray)/nrays;
    rads=1+mod([0:nrays*neachray-1],neachray);
    layouts{3}.coords=[[0; cos(angs(:)).*rads(:)],[0; sin(angs(:)).*rads(:)]];
    layouts{3}.keeprays=1;
    %
    nrays=4;
    neachray=6;
    layouts{4}.name=sprintf('%1.0f rays, %1.0f points each',nrays,neachray);
    angs=2*pi*floor([0:nrays*neachray-1]/neachray)/nrays;
    rads=1+mod([0:nrays*neachray-1],neachray);
    layouts{4}.coords=[[0; cos(angs(:)).*rads(:)],[0; sin(angs(:)).*rads(:)]];
    %
    ncirc=24;
    rads=3.;
    layouts{5}.name=sprintf('%2.0f point circle, radius %2.0f',ncirc,rads);
    angs=2.*pi*[0:ncirc-1]/ncirc;
    layouts{5}.coords=[[0; cos(angs(:)).*rads],[0; sin(angs(:)).*rads]];
    %
    nxy=[11 11];
    layouts{6}.name=sprintf('%1.0f x %1.0f grid (0.5 side length)',nxy);
    [xg,yg]=meshgrid([1:nxy(1)],[1:nxy(2)]);
    layouts{6}.coords=[xg(:),yg(:)]/2;
    %add extra dimensions, subtract mean, add needed fields
    for il=1:length(layouts)
        layouts{il}.npts=size(layouts{il}.coords,1);
        if size(layouts{il}.coords,2)<dim_max
            layouts{il}.coords=[layouts{il}.coords,zeros(layouts{il}.npts,dim_max-size(layouts{il}.coords,2))];
        end
        layouts{il}.coords=layouts{il}.coords-repmat(mean(layouts{il}.coords,1),layouts{il}.npts,1);
        layouts{il}=filldefault(layouts{il},'keeprays',0);
    end
end
nlayouts=length(layouts);
%
%set up transforms
%
if ~exist('transforms')
    transforms=cell(0);
    rotang12=pi/3;
    rotang13=pi/10;
    rot12=[cos(rotang12) sin(rotang12);-sin(rotang12) cos(rotang12)];
    rot13=[cos(rotang13) sin(rotang13);-sin(rotang13) cos(rotang13)];
    %
    rot=eye(dim_max);
    rot([1 2],[1 2])=rot12;
    if dim_max>=3
        rot3=eye(dim_max);
        rot3([1 3],[1 3])=rot13;
        rot=rot*rot3;
    end
    aff=eye(dim_max);
    proj=zeros(dim_max,1);
    vcut=zeros(1,dim_max);
    tdif=zeros(1,dim_max);
%   quantities specified up to 3d, extend beyond with zeros
    affratios=[0.6 1.4 1.2 0.7];
    projvec=[0.1 0.07 -0.3]';
    vcutvec=[.2 .8 .1];
    vcutvec=vcutvec./sqrt(sum(vcutvec.^2));
    tdifvec=[.5 0 1];
    for id=1:min(3,dim_max)
        aff(id,id)=affratios(id);
        proj(id)=projvec(id);
        vcut(id)=vcutvec(id);
        tdif(id)=tdifvec(id);
    end
    %
    transforms{1}.name='original';
    transforms{1}.class='procrustes';
    transforms{1}.params.T=eye(dim_max);
    transforms{1}.params.b=1;
    transforms{1}.params.c=zeros(1,dim_max);
    %
    transforms{2}=transforms{1};
    transforms{2}.name='procrustes_noscale';
    transforms{2}.params.T=rot;
    %
    transforms{3}=transforms{2};
    transforms{3}.name='procrustes_scale';
    transforms{3}.params.b=1.25;
    %
    transforms{4}=transforms{2};
    transforms{4}.name='affine_nooffset';
    transforms{4}.class='affine';
    transforms{4}.params.T=aff*rot;
    %
    transforms{5}=transforms{4};
    transforms{5}.name='projective';
    transforms{5}.class='projective';
    transforms{5}.params.p=proj;
    %
    transforms{6}=transforms{4};
    transforms{6}.name='piecewise affine';
    transforms{6}.class='pwaffine';
    transforms{6}.params.vcut=vcut;
    transforms{6}.params.acut=-0.5;
    transforms{6}.params.T=repmat(transforms{6}.params.T,[1 1 2]);
    transforms{6}.params.T(:,:,1)=transforms{6}.params.T(:,:,2)+vcut'*tdif; %will do nothing if orthog to vcut
    transforms{6}.params.c=zeros(2,dim_max);
end
ntransforms=length(transforms);
%set up structures needed for plotting un-transformed data
ds=cell(nlayouts,1);
sas=cell(nlayouts,1);
rayss=cell(nlayouts,1);
%
opts_rays_used=cell(nlayouts,1);
for il=1:nlayouts
    for id=1:dim_max
        ds{il}{id}=layouts{il}.coords(:,1:id);
    end
    sas{il}.typenames=cell(layouts{il}.npts,1);
    for ipt=1:layouts{il}.npts
        sas{il}.typenames{ipt}=sprintf('x=%5.2f y=%5.2f',layouts{il}.coords(ipt,[1 2]));
    end
    [rayss{il},opts_rays_used{il}]=psg_findrays(layouts{il}.coords,opts_rays);
    %remove rays (optionally) and rings
    if ~layouts{il}.keeprays
        rayss{il}.nrays=0;
        rayss{il}.whichray(rayss{il}.whichray>0)=NaN;
    end
    rayss{il}.nrings=0;
    rayss{il}.rings=cell(0);
    disp(sprintf(' layout %2.0f (%3.0f points):  %s',il,layouts{il}.npts,layouts{il}.name));
end
layout_list=getinp('layouts to use','d',[1 nlayouts],[1:nlayouts]);
ntransforms=length(transforms);
for it=1:ntransforms
    disp(sprintf(' %2.0f->%s',it,transforms{it}.name));
end
transform_list=getinp('transforms to show','d',[1 ntransforms],[1:ntransforms]);
%
%show the layouts, also to get useful defaults for opts_plot_used
opts_vis_used=cell(nlayouts,1);
opts_plot1_used=cell(nlayouts,1);
opts_mult_used=cell(nlayouts,1);
opts_plot2_used=cell(nlayouts,ntransforms);
for ilptr=1:length(layout_list)
    il=layout_list(ilptr);
    opts_vis.vis_string_format=layouts{il}.name;
    opts_plot.if_legend=0;
    [opts_vis_used{il},opts_plot_used_vis1{il},opts_mult_used{il}]=psg_visualize(plotformats,ds{il},sas{il},rayss{il},opts_vis,opts_plot,opts_mult);
end
%compute transforms
for ilptr=1:length(layout_list)
    il=layout_list(ilptr);
end
%plot transforms
[plot_rows,plot_cols]=nicesubp(length(transform_list),0.7);
plot_range=zeros(nlayouts,1);
for ilptr=1:length(layout_list)
    il=layout_list(ilptr);
    opts_plot_use=opts_plot_used_vis1{il}{1};
    opts_plot_use=rmfield(opts_plot_use,'xform_mult'); %since opts_plot_use might have been plotting a different number of dimensions
    opts_plot_use=rmfield(opts_plot_use,'xform_offset');
    opts_plot_use.plot_range=repmat(10,2,3);
    figure;
    set(gcf,'Position',[100 100 1200 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',layouts{il}.name);
    %
    for itptr=1:length(transform_list)
        it=transform_list(itptr);
        subplot(plot_rows,plot_cols,itptr);
        opts_plot_use=setfield(opts_plot_use,'axis_handle',gca);
        d_orig=ds{il}{dim_max};
        switch transforms{it}.class
            case {'mean','procrustes','affine'}
                d=transforms{it}.params.b*d_orig*transforms{it}.params.T+repmat(transforms{it}.params.c,layouts{il}.npts,1);
            case 'projective'
                d=persp_apply(transforms{it}.params.T,transforms{it}.params.c,transforms{it}.params.p,d_orig);
            case 'pwaffine'
                d=psg_pwaffine_apply(transforms{it}.params,d_orig);
            otherwise
                d=orig;
        end %model class
        opts_plot2_used{il,it}=psg_plotcoords(d,dim_select,sas{il},rayss{il},opts_plot_use);
        plot_range(il)=max([plot_range(il),max(abs(d(:))),max(opts_plot2_used{il,it}.plot_range(:))]);
        title(strrep(transforms{it}.name,'_',' '),'Interpreter','none');
    end %itptr
    %equal scale on all transforms
    for itptr=1:length(transform_list)
        it=transform_list(itptr);
        ha=opts_plot2_used{il,it}.axis_handle;
        axis equal;
        set(ha,'XLim',plot_range(il)*[-1 1]);
        set(ha,'YLim',plot_range(il)*[-1 1]);
        if length(dim_select)>2
            set(ha,'ZLim',plot_range(il)*[-1 1]);
            axis vis3d
            view(transform_view);
        end
    end
    axes('Position',[0.01,0.02,0.01,0.01]); %for text
    text(0,0,layouts{il}.name,'Interpreter','none','FontSize',8);
    axis off;
end %il

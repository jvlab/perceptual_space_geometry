%psg_geomodels_illus: demonstrate various geometric transforms with simple
%starting points, using standard visualization tools, and also shows fits
%
% could add noise
%
%  See also:  PSG_MODELS_DEFINE, PSG_PWAFFINE_APPLY, PERSP_APPLY, 
%    PSG_FINDRAYS, PSG_VISUALIZE, PSG_PLOTCOORDS, PSG_GEO_GENERAL,
%    PSG_GEOMODELS_DEFINE, PSG_GEO_LAYOUTS_SETUP, PSG_GEO_TRANSFORMS_SETUP.
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
%model fitting and display params
model_types_def=psg_geomodels_define();
model_types=model_types_def.model_types;
nmodels=length(model_types);
if ~exist('fit_colors') fit_colors='rmbcg';end
if ~exist('fit_offset') fit_offset=1; end %offset for plotting fits
%
if ~exist('layouts')
    layouts=psg_geo_layouts_setup(dim_max);
end
nlayouts=length(layouts);
%
%set up transforms
%
if ~exist('transforms')
    transforms=psg_geo_transforms_setup(dim_max);
end
ntransforms=length(transforms);
%
%set up structures needed for plotting un-transformed data
%
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
    disp(sprintf(' layout %2.0f (%3.0f points):  %s',il,layouts{il}.npts,layouts{il}.label));
end
layout_list=getinp('layouts to use','d',[1 nlayouts],[1:nlayouts]);
for it=1:ntransforms
    disp(sprintf(' %2.0f->%s',it,transforms{it}.label));
end
transform_list=getinp('transforms to simulate and show','d',[1 ntransforms],[1:ntransforms]);
%
transform_fit_list=getinp('transforms to use for fitting (0 for none)','d',[0 ntransforms]);
%
%show the layouts, also to get useful defaults for opts_plot_used
opts_vis_used=cell(nlayouts,1);
opts_plot1_used=cell(nlayouts,1);
opts_mult_used=cell(nlayouts,1);
opts_plot2_used=cell(nlayouts,ntransforms);
for ilptr=1:length(layout_list)
    il=layout_list(ilptr);
    opts_vis.vis_string_format=layouts{il}.label;
    opts_plot.if_legend=0;
    [opts_vis_used{il},opts_plot_used_vis1{il},opts_mult_used{il}]=psg_visualize(plotformats,ds{il},sas{il},rayss{il},opts_vis,opts_plot,opts_mult);
end
%
%compute and plot transforms and fits
%
%plot transforms
if ~exist('plot_aspect') plot_aspect=0.7; end
[plot_rows,plot_cols]=nicesubp(length(transform_list),plot_aspect);
plot_range=zeros(nlayouts,1);
transform_fit=cell(ntransforms,ntransforms,nlayouts);
coords_fit=cell(ntransforms,ntransforms,nlayouts);
d_fit=zeros(ntransforms,ntransforms,nlayouts);
for ilptr=1:length(layout_list)
    il=layout_list(ilptr);
    disp(' ');
    disp(sprintf('layout %s',layouts{il}.label));
    opts_plot_use=opts_plot_used_vis1{il}{1};
    opts_plot_use=rmfield(opts_plot_use,'xform_mult'); %since opts_plot_use might have been plotting a different number of dimensions
    opts_plot_use=rmfield(opts_plot_use,'xform_offset');
    opts_plot_use.plot_range=repmat(10,2,3);
    figure;
    set(gcf,'Position',[100 100 1200 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'name',layouts{il}.label);
    %
    for itptr=1:length(transform_list)
        it=transform_list(itptr);
        subplot(plot_rows,plot_cols,itptr);
        opts_plot_use=setfield(opts_plot_use,'axis_handle',gca);
        coords_orig=ds{il}{dim_max};
        switch transforms{it}.class
            case {'mean','procrustes','affine'}
                coords_transformed=transforms{it}.params.b*coords_orig*transforms{it}.params.T+repmat(transforms{it}.params.c,layouts{il}.npts,1);
            case 'projective'
                coords_transformed=persp_apply(transforms{it}.params.T,transforms{it}.params.c,transforms{it}.params.p,coords_orig);
            case 'pwaffine'
                coords_transformed=psg_pwaffine_apply(transforms{it}.params,coords_orig);
            otherwise
                warning('transformation class not recognized, no transformation applied.');
                coords_transformed=coords_orig;
        end %model class
        opts_plot2_used{il,it}=psg_plotcoords(coords_transformed,dim_select,sas{il},rayss{il},opts_plot_use);
        plot_range(il)=max([plot_range(il),max(abs(coords_transformed(:))),max(opts_plot2_used{il,it}.plot_range(:))]);
        %
        if any(transform_fit_list>0)
            disp(' ');
        end
        for ifitptr=1:length(transform_fit_list)
            ifit=transform_fit_list(ifitptr);
            if ifit~=0
                model_type=transforms{ifit}.model_type;
                opts_model=model_types_def.(model_type).opts;
                model_class=model_types_def.(model_type).class;
                opts_model=model_types_def.(model_type).opts;
                [d_fit(it,ifit,il),coords_fit{it,ifit,il},transform_fit{it,ifit,il},opts_model_used]=psg_geo_general(coords_transformed,coords_orig,model_class,opts_model);
                disp(sprintf(' transforming with %30s, fitting with model type %30s: d=%8.5f',transforms{it}.label,model_type,d_fit(it,ifit,il)));
                mcolor=fit_colors(mod(ifitptr-1,length(fit_colors))+1);
                hold on;
                opts_plot_fit_use=opts_plot_use;
                opts_plot_fit_use=setfield(opts_plot_fit_use,'color_norays',mcolor);
                opts_plot_fit_use=setfield(opts_plot_fit_use,'color_origin',mcolor);
                opts_plot_fit_use=setfield(opts_plot_fit_use,'color_nearest_nbr',mcolor);
                opts_plot_fit_use=setfield(opts_plot_fit_use,'if_just_data',mcolor);
                coords_fit_plot=coords_fit{it,ifit,il};
                if length(dim_select)>2
                    coords_fit_plot(:,3)=coords_fit_plot(:,3)+fit_offset*ifitptr;
                end
                ou=psg_plotcoords(coords_fit_plot,dim_select,sas{il},rayss{il},opts_plot_fit_use);
                plot_range(il)=max(plot_range(il),max(ou.plot_range(:)));
            end
        end
        %
        title(strrep(transforms{it}.label,'_',' '),'Interpreter','none');
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
    text(0,0,layouts{il}.label,'Interpreter','none','FontSize',8);
    axis off;
end %il

function [opts_vis_used,opts_plot_used,opts_mult_used]=psg_visualize(plotformats,d,sa,rays,opts_vis,opts_plot,opts_mult)
% [opts_vis_used,opts_plot_used,opts_mult_used]=psg_visualize(plotformats,d,sa,rays,opts_vis,opts_plot) plots several pages
% of visualizations of psg results, one for each dimension
%
% plotformats: one or more rows of [dim dims_together]: which dimension model to plot, and how many dims to plot together
% d: d{idim} is a model with idim dimensions
% sa: setup structure, from psg_read_coorddata
% rays: a ray structure, typically from psg_findrays
% opts_vis: options
%   opts_vis.if_plotrays=if_plotrays; 1 to plot fitted (unidirectional) rays, defaults to ~isempty(opts_vis.d_rayfit)
%   opts_vis.if_plotbids=if_plotbids; 1 to plot fitted bidirectional rays, defaults to ~isempty(opts_vis.d_bidfit)
%   opts_vis.d_rayfit: ray data in same format as d, from psg_rayfit
%   opts_vis.d_bidfit: bidirectional ray data, in same format as d, from psg_rayfit
%   opts_vis.vis_string_format='raw %1.0fd fit'; how dimension label is formatted
%   opts_vis.file_string: file name string
%   opts_vis.which_dimcombs
%       all': (default) plot all combinations
%      'keeplow': keep all but one dimensions low and only step the highest; plotformat=[5 3] yields [1 2 3],[1 2 4],[1 2 5]
%      'keepone': keep one dimension and step the rest; plotformat=[5 3] yields [1 2 3],[1 2 4],[1 2 5],[1 3 4],[1 3 5],[1 4 5]
%      'rolling': rolling contiguous subsets; plotformat=[5 3] yields [1 2 3],[2 3 4],[3 4 5],[4 5 1],[5 1 2]
%      'onlylowest': only the lowest dimensions, plotformat=[5 3] yields [1 2 3]
%   opts_vis.if_pcrot: use the principal components of each mode d{idim} as axes
%       note that principal components are determined for all of the d{idim} model, not just the combination of dimensions being shown
%       this prefixes 'pc ' onto opts_plot.axis_label_prefix
%   opts_vis.offset: a vector of length length(d), value to plot at origin
% opts_plot: options for psg_plotcoords, can be omitted
%       note that xform_offset is determined by xform_mult are computed from
% opts_mult: options for plotting multiple sets
%   opts_mult.line_widths: list of line widths
%
%  d, sa, rays, and opts_vis may also be cell arrays (1,nmult) of the structures described above. 
%    The setups (sa) should be consistent with each other, but this is not checked.
%    opts_vis{1} takes priority for which_dimcombs, maxcomb_legend, and prepending pca to dimension label;
%     other elements of opts_vis{*} set individually and combined intelligently (e.g., file name strings are concatenated)
%
%   opts_vis_used: options used, and handle to figure
%   opts_plot_used: options used from main call to psg_plotcoords
%   opts_mult_used: options usef for opts_mult
%
%  27Dec22: allow d, sa,rays, opts_vis to each be cell arrays whose subsidiary srtructures are as above; add opts_mult
%
%   See also: PSG_FINDRAYS, PSG_RAYFIT, PSG_PLOTCOORDS, PSG_VISUALIZE_DEMO, PSG_QFORMPRED_DEMO, PSG_PLOTANGLES, ISEMPTYSTRUCT.
%

%determine if multiple datasets will be superimposed,
%and set up dm, sam, raysm to be cell arrays 
if iscell(d{1})
    nmult=length(d);
    dm=d;
    sam=sa;
    raysm=rays;
    opts_vism=opts_vis;
else
    nmult=1;
    dm{1}=d;
    sam{1}=sa;
    raysm{1}=rays;
    opts_vism{1}=opts_vis;
end
if (nargin<5)
    opts_vism=cell(1,nmult);
end
for im=1:nmult
    opts_vism{im}=filldefault(opts_vism{im},'d_rayfit',struct());
    opts_vism{im}=filldefault(opts_vism{im},'d_bidfit',struct());
    opts_vism{im}=filldefault(opts_vism{im},'if_plotrays',~isemptystruct(opts_vism{im}.d_rayfit));
    opts_vism{im}=filldefault(opts_vism{im},'if_plotbids',~isemptystruct(opts_vism{im}.d_bidfit));
    opts_vism{im}=filldefault(opts_vism{im},'file_string',[]);
    opts_vism{im}=filldefault(opts_vism{im},'vis_string_format','raw %1.0fd fit');
    opts_vism{im}=filldefault(opts_vism{im},'lim_margin',1.1);
    opts_vism{im}=filldefault(opts_vism{im},'maxcomb_legend',6); %maximum number of subplots for which a legend is present
    opts_vism{im}=filldefault(opts_vism{im},'which_dimcombs','all');
    opts_vism{im}=filldefault(opts_vism{im},'if_pcrot',0);
end
file_string_cat=[];
for im=1:nmult
    file_string_cat=cat(2,file_string_cat,' ',opts_vism{im}.file_string);
end
file_string_cat=deblank(file_string_cat);
%
if (nargin<6)
    opts_plot=struct();
end
if (nargin<7)
    opts_mult=struct();
end
opts_mult=filldefault(opts_mult,'line_widths',[1:nmult]);
opts_mult_used=opts_mult;
%
model_dims_each=cell(1,nmult);
for im=1:nmult
    model_dims_each{im}=[];
    for idim=1:length(dm{im})
        if ~isempty(dm{im}{idim})
            model_dims_each{im}=[model_dims_each{im},size(dm{im}{idim},2)];
        end
    end
    if (im==1)
        model_dims=model_dims_each{im};
    else
        model_dims=intersect(model_dims,model_dims_each{im});
    end
end
for im=1:nmult
    if ~isfield(opts_vism{im},'offset')
        opts_vism{im}.offset=zeros(1,max(model_dims));
    end
end
rot=cell(1,nmult);
%
opts_plot=filldefault(opts_plot,'axis_label_prefix','dim');
axis_label_prefix_orig=opts_plot.axis_label_prefix;
opts_plot_used=cell(0);
for iplot=1:size(plotformats,1)
    model_dim=plotformats(iplot,1);
    dims_together=plotformats(iplot,2);
    if ismember(model_dim,model_dims) & model_dim>=dims_together
        for im=1:nmult
            model_dim_ptr(im)=find(model_dims_each{im}==model_dim); %needed since not all model dimensions are present
        end
        vis_string=sprintf(opts_vism{1}.vis_string_format,model_dim);
        for im=1:nmult
            if (opts_vism{im}.if_pcrot) %do a rotation into principal components before choosing a set of dims to plot
                opts_plot.axis_label_prefix=cat(2,'pc ',axis_label_prefix_orig);
                pca_offset=mean(dm{im}{model_dim_ptr(im)},1);
                [recon_pcaxes,recon_coords,var_ex,var_tot,coord_maxdiff,opts_pca_used]=psg_pcaoffset(dm{im}{model_dim_ptr(im)},pca_offset);
                rot{im}=opts_pca_used.qv;
            else
                rot{im}=eye(model_dim);
            end
        end
        opts_vism{1}.fig_handle{iplot}=figure;
        for im=1:nmult
            opts_vism{im}.fig_handle(iplot)=opts_vism{1}.fig_handle(iplot);
        end
        switch opts_vism{1}.which_dimcombs
            case 'all'
                dim_combs=nchoosek([1:model_dim],dims_together);
            case 'keeplow'
                dim_combs=[repmat([1:dims_together-1],model_dim-dims_together+1,1),[dims_together:model_dim]'];
            case 'keepone'
                dim_combs_hi=nchoosek([2:model_dim],dims_together-1);
                dim_combs=[ones(size(dim_combs_hi,1),1) dim_combs_hi];
            case 'rolling'
                dim_combs=mod(repmat([0:dims_together-1],model_dim,1)+repmat([0:model_dim-1]',1,dims_together),model_dim)+1;
            case 'onlylowest'
                dim_combs=[1:dims_together];
        end
        if size(dim_combs,1)>0
            ncombs=size(dim_combs,1);
            set(gcf,'Position',[50 50 1200 800]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',cat(2,file_string_cat,' ',vis_string));
            [nr,nc]=nicesubp(ncombs,0.7);
            dm_max=0;
            lim_max=0;
            for im=1:nmult
                dm_max=max(dm_max,max(abs(dm{im}{model_dim_ptr(im)}(:))));
                lim_max=max(lim_max,opts_vism{im}.lim_margin);
            end
            opts_plot.lims=lim_max*[-1 1]*dm_max;
            if isfield(opts_plot,'if_legend')
                if_legend=opts_plot.if_legend;
            else
               if_legend=double(ncombs<=opts_vism{1}.maxcomb_legend);
            end
            opts_plot_use=opts_plot;
            opts_plot_use.if_legend=if_legend;
            for icomb=1:ncombs
                ha=subplot(nr,nc,icomb);
                for im=1:nmult
                    if (nmult>1)
                        opts_plot_use=setfield(opts_plot_use,'line_width',opts_mult.line_widths(im));
                    end
                    opts_plot_use.xform_offset=opts_vism{im}.offset(1:model_dim);
                    opts_plot_used{iplot,icomb}=psg_plotcoords(dm{im}{model_dim_ptr(im)}*rot{im},dim_combs(icomb,:),sam{im},raysm{im},setfield(opts_plot_use,'axis_handle',ha));
                    ha=opts_plot_used{iplot,icomb}.axis_handle;
                    %
                    if opts_vism{im}.if_plotrays
                        psg_plotcoords(opts_vism{im}.d_rayfit{model_dim}*rot{im},dim_combs(icomb,:),sam{im},raysm{im},setfields(opts_plot_use,{'axis_handle','line_type','if_just_data'},{ha,':',1}));
                    end
                    if opts_vism{im}.if_plotbids
                        psg_plotcoords(opts_vism{im}.d_bidfit{model_dim}*rot{im},dim_combs(icomb,:),sam{im},raysm{im},setfields(opts_plot_use,{'axis_handle','line_type','if_just_data','if_origin_on_rays'},{ha,'--',1,0}));
                    end
                    opts_plot_use.if_legend=0; %after im=1, turn off lsegend
                end %im
            end %icomb
            axes('Position',[0.01,0.04,0.01,0.01]); %for text
            text(0,0,vis_string,'Interpreter','none','FontSize',8);
            axis off;
            axes('Position',[0.01,0.02,0.01,0.01]); %for text
            text(0,0,file_string_cat,'Interpreter','none','FontSize',8);
            axis off;
        end
    end %if ismember
end
opts_vis_used=opts_vism;
return


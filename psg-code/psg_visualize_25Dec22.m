function [opts_vis_used,opts_plot_used]=psg_visualize(plotformats,d,sa,rays,opts_vis,opts_plot)
% [opts_vis_used,opts_plot_used]=psg_visualize(plotformats,d,sa,rays,opts_vis,opts_plot) plots several pages
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
%   opts_vis.if_pcrot: use the principal components of each mode d{idim} as axes
%       note that principal components are determined for all of the d{idim} model, not just the combination of dimensions being shown
%       this prefixes 'pc ' onto opts_plot.axis_label_prefix
%   opts_vis.offset: a vector of length length(d), value to plot at origin
% opts_plot: options for psg_plotcoords, can be omitted
%       note that xform_offset is determined by xform_mult are computed from 
%
%   opts_vis_used: options used, and handle to figure
%   opts_plot_used: options used from main call to psg_plotcoords
%
%   See also: PSG_FINDRAYS, PSG_RAYFIT, PSG_PLOTCOORDS, PSG_VISUALIZE_DEMO, PSG_QFORMPRED_DEMO, PSG_PLOTANGLES.
%
if (nargin<5)
    opts_vis=struct();
end
opts_vis=filldefault(opts_vis,'d_rayfit',struct());
opts_vis=filldefault(opts_vis,'d_bidfit',struct());
opts_vis=filldefault(opts_vis,'if_plotrays',~isempty(opts_vis.d_rayfit));
opts_vis=filldefault(opts_vis,'if_plotbids',~isempty(opts_vis.d_bidfit));
opts_vis=filldefault(opts_vis,'file_string',[]);
opts_vis=filldefault(opts_vis,'vis_string_format','raw %1.0fd fit');
opts_vis=filldefault(opts_vis,'lim_margin',1.1);
opts_vis=filldefault(opts_vis,'maxcomb_legend',6); %maximum number of subplots for which a legend is present
opts_vis=filldefault(opts_vis,'which_dimcombs','all');
opts_vis=filldefault(opts_vis,'if_pcrot',0);
%%
if (nargin<6)
    opts_plot=struct();
end
model_dims=[];
for idim=1:length(d)
    if ~isempty(d{idim})
        model_dims=[model_dims,idim];
    end
end
if ~isfield(opts_vis,'offset')
    opts_vis.offset=zeros(1,max(model_dims));
end
opts_plot=filldefault(opts_plot,'axis_label_prefix','dim');
axis_label_prefix_orig=opts_plot.axis_label_prefix;
opts_plot_used=cell(0);
for iplot=1:size(plotformats,1)
    model_dim=plotformats(iplot,1);
    dims_together=plotformats(iplot,2);
    if ismember(model_dim,model_dims) & model_dim>=dims_together
        vis_string=sprintf(opts_vis.vis_string_format,model_dim);
        if (opts_vis.if_pcrot) %do a rotation into principal components before choosing a set of dims to plot
            opts_plot.axis_label_prefix=cat(2,'pc ',axis_label_prefix_orig);
            pca_offset=mean(d{model_dim},1);
            [recon_pcaxes,recon_coords,var_ex,var_tot,coord_maxdiff,opts_pca_used]=psg_pcaoffset(d{model_dim},pca_offset);
            rot=opts_pca_used.qv;
        else
            rot=eye(model_dim);
        end
        opts_vis.fig_handle{iplot}=figure;
        switch opts_vis.which_dimcombs
            case 'all'
                dim_combs=nchoosek([1:model_dim],dims_together);
            case 'keeplow'
                dim_combs=[repmat([1:dims_together-1],model_dim-dims_together+1,1),[dims_together:model_dim]'];
            case 'keepone'
                dim_combs_hi=nchoosek([2:model_dim],dims_together-1);
                dim_combs=[ones(size(dim_combs_hi,1),1) dim_combs_hi];
            case 'rolling'
                dim_combs=mod(repmat([0:dims_together-1],model_dim,1)+repmat([0:model_dim-1]',1,dims_together),model_dim)+1;
        end
        if size(dim_combs,1)>0
            ncombs=size(dim_combs,1);
            set(gcf,'Position',[50 50 1200 800]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',cat(2,opts_vis.file_string,' ',vis_string));
            [nr,nc]=nicesubp(ncombs,0.7);
            opts_plot.lims=opts_vis.lim_margin*[-1 1]*max(abs(d{model_dim}(:)));
            if isfield(opts_plot,'if_legend')
                if_legend=opts_plot.if_legend;
            else
               if_legend=double(ncombs<=opts_vis.maxcomb_legend);
            end
            opts_plot_use=opts_plot;
            opts_plot_use.if_legend=if_legend;
            for icomb=1:ncombs
                ha=subplot(nr,nc,icomb);
                opts_plot_used{iplot,icomb}=psg_plotcoords(d{model_dim}*rot,dim_combs(icomb,:),sa,rays,setfield(opts_plot_use,'axis_handle',ha));
                ha=opts_plot_used{iplot,icomb}.axis_handle;
                %
                if opts_vis.if_plotrays
                    psg_plotcoords(opts_vis.d_rayfit{model_dim}*rot,dim_combs(icomb,:),sa,rays,setfields(opts_plot_use,{'axis_handle','line_type','if_just_data'},{ha,':',1}));
                end
                if opts_vis.if_plotbids
                    psg_plotcoords(opts_vis.d_bidfit{model_dim}*rot,dim_combs(icomb,:),sa,rays,setfields(opts_plot_use,{'axis_handle','line_type','if_just_data','if_origin_on_rays'},{ha,'--',1,0}));
                end
            end
            axes('Position',[0.01,0.04,0.01,0.01]); %for text
            text(0,0,vis_string,'Interpreter','none');
            axis off;
            axes('Position',[0.01,0.02,0.01,0.01]); %for text
            text(0,0,opts_vis.file_string,'Interpreter','none');
            axis off;
        end
    end %if ismember
end
opts_vis_used=opts_vis;
return

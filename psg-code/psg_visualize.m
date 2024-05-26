function [opts_vis_used,opts_plot_used,opts_mult_used]=psg_visualize(plotformats,d,sa,rays,opts_vis,opts_plot,opts_mult)
% [opts_vis_used,opts_plot_used,opts_mult_used]=psg_visualize(plotformats,d,sa,rays,opts_vis,opts_plot,opts_mult) plots several pages
% of visualizations of psg results, one for each dimension.
%
% See psg_visualize_demo for examples that plot best-fitting rays, and psg_consensus_demo for examples that plot connections
%   between points of corresponding datasets.
%
% plotformats: one or more rows of [dim dims_together]: which dimension model to plot, and how many dims to plot together
% d: d{idim} is a model with idim dimensions
% sa: setup structure, from psg_read_coorddata
% rays: a ray structure, typically from psg_findrays
% opts_vis: options
%   opts_vis.if_plotrays: 1 to plot fitted (unidirectional) rays, defaults to ~isempty(opts_vis.d_rayfit)
%   opts_vis.if_plotbids: 1 to plot fitted bidirectional rays, defaults to ~isempty(opts_vis.d_bidfit)
%   opts_vis.d_rayfit: ray data in same format as d, from psg_rayfit
%   opts_vis.d_bidfit: bidirectional ray data, in same format as d, from psg_rayfit
%   opts_vis.vis_string_format='raw %1.0fd fit'; how dimension label is formatted
%   opts_vis.file_string: file name string
%   opts_vis.which_dimcombs
%       'all': (default) plot all combinations
%      'keeplow': keep all but one dimensions low and only step the highest; plotformat=[5 3] yields [1 2 3],[1 2 4],[1 2 5]
%      'keepone': keep one dimension and step the rest; plotformat=[5 3] yields [1 2 3],[1 2 4],[1 2 5],[1 3 4],[1 3 5],[1 4 5]
%      'rolling': rolling contiguous subsets; plotformat=[5 3] yields [1 2 3],[2 3 4],[3 4 5],[4 5 1],[5 1 2]
%      'onlylowest': only the lowest dimensions, plotformat=[5 3] yields [1 2 3]
%   opts_vis.if_pcrot: use the principal components of each mode d{idim} as axes
%       note that principal components are determined for all of the d{idim} model, not just the combination of dimensions being shown
%       this prefixes 'pc ' onto opts_plot.axis_label_prefix
%   opts_vis.offset: a vector of length length(d), value to plot at origin
%   opts_vis.offset_ptr: integer, if > 0, points to the condition to be plotted at origin: d{idim}(opts_vis.offset_ptr,:)
%                        if -1, use centroid
%   opts_vis.offset_norot: 0: offset is affected by pca rotation
%                          1: offset is independent of pca rotation
%                         -1: offset to data point (from offset_ptr) is affected by pca rotation but abaolute offset (opts_vis.offset) is not
%   opts_vis.tet_signs: 4-column array of sign choices for 4-d tetrahedral plots, defaults to [1 1 1 1], ignored if not 4-d
% opts_plot: options for psg_plotcoords, can be omitted
%   opts_plot.xform_offset is ignored, as it is determined by opts_vis.offset or opts_vis.offset_ptr and therefore can be separate for each datatset
%   opts_plot.marker_size: marker size
% opts_mult: options for plotting multiple sets
%   opts_mult.line_widths: list of line widths for individual datasets
%   opts_mult.connect_specs: list [nconnect 2] of corresponding datasets to
%      connect_specs=[1 2;1 3;1 4] or 'star' connects dataset 1 with each of 2,3,4
%      connect_specs=[1 2;2 3;3 1] or 'circuit' connects dataset 1 with 2, 2 with 3, 3 with 1
%      defaults to empty; used to generate connect_list, numeric list
%   opts_mult.connect_line_width: line width used for connection, defaults to 1
%   opts_mult.connect_line_type: line type used for connection, defaults to '-'
%   opts_mult.connect_line_type_neg: line type used for connection on negative coords, defaults to connect_line_type
%   opts_mult.connect_only: 1 if only the connections are drawn, defaults to 0
%   opts_mult.if_pcrot_whichuse: 0 (default) for each dataset to use its
%      own pca (if opts_vis{im}.if_pcrot=1; otherwise, indicates which pc to use for a common rotation
%   opts_mult.if_fit_range: 1: sets axis limits to range plotted (defaults to 0, which chooses equal round numbers)
%   opts_mult.color_norays_list: if present and opts_plot.if_use_rays=0, and connect_list is empty, 
%      a list of colors for each dataset; a cell array, e.g., {'b',[0.3 0.4 1],'k'};
%   opts_mult.color_norays_connect_mode: how segments between datasets are colored, if no rays
%     0->use origin color, 1->use color of dataset in connect_specs(:,1), 2[default]->use color of dataset in connect_specs(:,2)
%
%  d, sa, rays, and opts_vis may also be cell arrays (1,nmult) of the structures described above. 
%    The setups (sa) should be consistent with each other, but only typenames and spec_labels of sa are checked.
%    opts_vis{1} takes priority for which_dimcombs, maxcomb_legend, and prepending pca to dimension label;
%    other elements of opts_vis{*} set individually and combined intelligently (e.g., file name strings are concatenated)
%
%   opts_vis_used{im}: im in [1 nmult]: options used, and handle to figure
%   opts_plot_used{iplot,icomb}: options used from main call to psg_plotcoords, iplot is row of plotformats, icomb is subplot
%   opts_mult_used: options used for opts_mult, including connect_list
%
%  27Dec22: allow d, sa,rays, opts_vis to each be cell arrays whose subsidiary structures are as above; add opts_mult
%  14Jan23: add opts_mult.connect, to draw lines connecting corresponding points in multiple datasests
%           add opts_vis.offset_ptr, opts_vis.offset_norot, opts_mult.if_pcrot_whichuse
%  23Jan23: offset_use combines offset_ptr and an explicit offset, rather than mutually exclusive
%  03Feb23: add opts_vis.tet_signs
%  24Mar23: clarify comment on relation of offset and pca rotation, add -1 option to offset_norot
%  31Oct23: add opts_plot_used{iplot,icomb}.plot_range (array of size [2 3]) to indicate range plotted
%  31Oct23: add opts_mult.if_fit_range 
%  14Nov23: modify psg_visualize_check to avoid checking fields that don't exist
%  13Dec23: bug fix in legend cleanup
%  27Apr24: bug fix in legend cleanup when connect_list is present, add opts_mult.connect_line_type_neg
%  10May24: additional fill-ins for empty opts_vis input
%  14May24: add opts_mult.color_norays_list
%  24May24: further fixes for empty opts_vis input, pcaoffset now calculated with omitnan
%  25May24: data centroid now calculated with omitnan
%  26May24: add opts_mult.color_noays_connect_mode
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
    warn_msg=psg_visualize_check(sa);
    if ~isempty(warn_msg)
        disp(warn_msg);
    end
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
if isempty(opts_vis) %in case opts_vis is expliclty input, but as empty or struct();
    opts_vism=cell(1,nmult);
elseif isstruct(opts_vis)
    if isempty(fieldnames(opts_vis))
        opts_vism=cell(1,nmult);
    end
end
if length(opts_vism)<nmult
    for im=length(opts_vism)+1:nmult
        opts_vism{im}=struct;
    end
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
    opts_vism{im}=filldefault(opts_vism{im},'tet_signs',[1 1 1 1]);
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
opts_mult=filldefault(opts_mult,'connect_specs',[]);
opts_mult=filldefault(opts_mult,'connect_line_width',1);
opts_mult=filldefault(opts_mult,'connect_line_type','-');
opts_mult=filldefault(opts_mult,'connect_line_type_neg',opts_mult.connect_line_type);
opts_mult=filldefault(opts_mult,'connect_only',0);
opts_mult=filldefault(opts_mult,'if_pcrot_whichuse',0);
opts_mult=filldefault(opts_mult,'if_fit_range',0);
opts_mult=filldefault(opts_mult,'color_norays_list',[]);
opts_mult=filldefault(opts_mult,'color_norays_connect_mode',2); %use second color of connect list
%
[connect_list,warn_msg]=psg_visualize_getconnect(opts_mult.connect_specs,nmult);
if ~isempty(warn_msg)
    disp(warn_msg);
    disp('connect_specs');
    disp(opts_mult.connect_specs);
end
opts_mult.connect_list=connect_list;
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
    opts_vism{im}=filldefault(opts_vism{im},'offset',zeros(1,max(model_dims)));
    opts_vism{im}=filldefault(opts_vism{im},'offset_ptr',0);
    opts_vism{im}=filldefault(opts_vism{im},'offset_norot',0);
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
                if (opts_mult.if_pcrot_whichuse==0) %which dataset to use for pca
                    impc=im;
                else
                    impc=opts_mult.if_pcrot_whichuse;
                end
                pca_offset=mean(dm{impc}{model_dim_ptr(impc)},1,'omitnan');
                [recon_pcaxes,recon_coords,var_ex,var_tot,coord_maxdiff,opts_pca_used]=psg_pcaoffset(dm{impc}{model_dim_ptr(impc)},pca_offset);
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
        ncombs=size(dim_combs,1);
        if (dims_together==4)
            ntet_signs=size(opts_vism{1}.tet_signs,1);
        else
            ntet_signs=1;
        end
        if ncombs>0
            set(gcf,'Position',[50 50 1200 800]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',cat(2,file_string_cat,' ',vis_string));
            [nr,nc]=nicesubp(ncombs*ntet_signs,0.7);
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
            %
            %table of coords for plotting
            %offsets taken into account here, as this may depend on number of dimensions if offset_ptr>0
            %or if rot is nontrival because of pca
            %
            nconds=size(dm{1}{model_dim_ptr(1)},1);
            coords_all=zeros(nconds,model_dim,nmult); %without offset
            coords_all_offset=zeros(nconds,model_dim,nmult); %with offset
            offsets=zeros(1,model_dim,nmult);
            for im=1:nmult
                coords_orig=dm{im}{model_dim_ptr(im)};
                coords_all(:,:,im)=coords_orig*rot{im};
                if opts_vism{im}.offset_ptr>0
                    offset_use_data=coords_orig(opts_vism{im}.offset_ptr,:); %if needed, pca rotation is applied to offset later
                elseif opts_vism{1}.offset_ptr==-1
                    offset_use_data=mean(coords_orig,1,'omitnan'); %centroid, omitnan added 25May24
                else
                    offset_use_data=zeros(1,model_dim);
                end
                switch opts_vism{im}.offset_norot
                    case 1 %pca rotation has no effect on offset
                        offset_use=offset_use_data+opts_vism{im}.offset(1:model_dim);
                    case 0 %pca rotation affects both parts
                        offset_use=(offset_use_data+opts_vism{im}.offset(1:model_dim))*rot{im};
                    case -1 %pca rotation just affects data
                        offset_use=offset_use_data*rot{im}+opts_vism{im}.offset(1:model_dim);
                end
                offsets(1,:,im)=offset_use;
                coords_all_offset(:,:,im)=coords_all(:,:,im)-repmat(offset_use,nconds,1);
            end
            %
            for icomb_signs=1:ncombs*ntet_signs
                icomb=ceil(icomb_signs/ntet_signs); %outer loop
                tet_signs_ptr=mod(icomb_signs-1,ntet_signs)+1; %inner loop
                if (dims_together==4)
                    opts_plot_use.tet_signs=opts_vism{1}.tet_signs(tet_signs_ptr,:);
                end
                ha=subplot(nr,nc,icomb_signs);
                opts_plot_used{iplot,icomb}=struct;
                if (opts_mult.connect_only==0)
                    for im=1:nmult
                        if (nmult>1)
                            opts_plot_use=setfield(opts_plot_use,'line_width',opts_mult.line_widths(im));
                        end
                        if ~isempty(opts_mult.color_norays_list)
                            ncolors=length(opts_mult.color_norays_list);
                            opts_plot_use=setfield(opts_plot_use,'color_norays',opts_mult.color_norays_list{1+mod(im-1,ncolors)});
                        end
                        %only show tetrahedron after all data have been plotted
                        if_tet_show=[];
                        if (im==nmult)
                            if_tet_show=Inf;
                        end
                        opts_plot_use.xform_offset=zeros(1,model_dim); %offset handled above
                        opts_plot_use.tag_text=sprintf('ds %1.0f',im);
                        %
                        if ~isfield(opts_plot_use,'label_sets')
                            label_sets_use=0;
                        else
                            label_sets_use=double(ismember(im,opts_plot_use.label_sets));
                        end
                        opu=psg_plotcoords(coords_all_offset(:,:,im),dim_combs(icomb,:),sam{im},raysm{im},...
                            setfields(opts_plot_use,{'axis_handle','if_just_data','if_tet_show','label_sets'},{ha,double(im>1),if_tet_show,label_sets_use}));
                        %
                        if isempty(fieldnames(opts_plot_used{iplot,icomb}));
                            opts_plot_used{iplot,icomb}=opu;
                            ha=opts_plot_used{iplot,icomb}.axis_handle;
                        end
                        opts_plot_used{iplot,icomb}=psg_visualize_range(opts_plot_used{iplot,icomb},opu);
                        %for plotting connections, include rotation and offset
                        %
                        if opts_vism{im}.if_plotrays
                            opts_plot_use.tag_text=sprintf('fit_ray %1.0f',im);
                            opu=psg_plotcoords(opts_vism{im}.d_rayfit{model_dim}*rot{im},dim_combs(icomb,:),sam{im},raysm{im},...
                                setfields(opts_plot_use,{'axis_handle','line_type','if_just_data','xform_offset','if_rings'},{ha,':',1,offsets(1,:,im),0}));
                            opts_plot_used{iplot,icomb}=psg_visualize_range(opts_plot_used{iplot,icomb},opu);
                        end
                        if opts_vism{im}.if_plotbids
                            opts_plot_use.tag_text=sprintf('fit_bid %1.0f',im);
                            opu=psg_plotcoords(opts_vism{im}.d_bidfit{model_dim}*rot{im},dim_combs(icomb,:),sam{im},raysm{im},...
                                setfields(opts_plot_use,{'axis_handle','line_type','if_just_data','if_origin_on_rays','xform_offset','if_rings'},{ha,'--',1,0,offsets(1,:,im),0}));
                            opts_plot_used{iplot,icomb}=psg_visualize_range(opts_plot_used{iplot,icomb},opu);
                        end
                    end %im
                end %if_connect_only==0
                if size(connect_list,1)>0
                    opts_plot_use.tag_text='connection';
                    opts_plot_connect=setfields(opts_plot_use,{'axis_handle','if_just_data','connect_list','line_width','line_type','line_type_connect_neg'},...
                        {ha,1-opts_mult.connect_only,connect_list,opts_mult.connect_line_width,opts_mult.connect_line_type,opts_mult.connect_line_type_neg});
                    if ~isempty(opts_mult.color_norays_list) & opts_plot.if_use_rays==0 %if no rays and colors for each dataset, use the colors for the connections too
                        opts_plot_connect.color_connect_sets_norays=opts_mult.color_norays_list;
                        opts_plot.color_norays_connect_mode=opts_mult.color_norays_connect_mode;
                    end
                    opu=psg_plotcoords(coords_all_offset,dim_combs(icomb,:),sam{im},raysm{im},opts_plot_connect);
                    if isempty(fieldnames(opts_plot_used{iplot,icomb}))
                        opts_plot_used{iplot,icomb}=opu;
                    end
                    opts_plot_used{iplot,icomb}=psg_visualize_range(opts_plot_used{iplot,icomb},opu);
                    if  ~isempty(opts_mult.color_norays_list) & opts_plot.if_use_rays==0 
                        for im=1:nmult
                            opts_plot_use.tag_text=sprintf('replot ds %1.0f',im); %replot data points if each has its own color
                            ncolors=length(opts_mult.color_norays_list);
                            opts_plot_use=setfield(opts_plot_use,'color_norays',opts_mult.color_norays_list{1+mod(im-1,ncolors)});
                            opu=psg_plotcoords(coords_all_offset(:,:,im),dim_combs(icomb,:),sam{im},raysm{im},...
                                setfields(opts_plot_use,{'axis_handle','if_just_data','if_tet_show','label_sets'},{ha,1,0,0}));
                        end
                    end
                end
                opts_plot_use.if_legend=0; %after icomb=1, turn off legend
                %clean up legend
                if (icomb_signs==1) %only put legend in first panel of each figure
                    hc=get(ha,'Children');
                    tags=cell(length(hc),1);
                    for ich=1:length(hc)
                        tags{ich}=get(hc(ich),'Tag');
                    end
                    hc_set1=find(contains(tags,'ds 1'));
                    hc_rays=find(contains(tags,'ray'));
                    hc_sign=find(contains(tags,'signed'));
                    hc_conn=find(contains(tags,'connection'));
                    hc_c1=find(contains(tags,'connect   1'));
                    hc_p1=find(contains(tags,'point   1'));
                    %
                    hc_keep=intersect(intersect(hc_set1,hc_rays),hc_sign);
                    if opts_mult.connect_only
                        hc_keep_conn=intersect(intersect(hc_c1,hc_p1),hc_conn);
                        hc_keep=intersect(hc_keep_conn,hc_conn);
                    end
                    if opts_plot_use.if_legend | size(connect_list,1)>0 %test for if_legend added 13Dec23, connect_list added 27Apr24
                        legend(hc(hc_keep));
                    end
                    if isfield(opts_plot,'if_use_rays')
                        if opts_plot.if_use_rays==0
                            legend off;
                        end
                    end
                end
                if opts_mult.if_fit_range
                    plot_range=opts_plot_used{iplot,icomb}.plot_range;
                    set(gca,'XLim',plot_range(:,1));
                    set(gca,'YLim',plot_range(:,2));
                    axis equal;
                    if size(plot_range,2)==3
                        set(gca,'ZLim',plot_range(:,3));
                    end
                end
            end %icomb
            %
            axes('Position',[0.01,0.04,0.01,0.01]); %for text
            text(0,0,vis_string,'Interpreter','none','FontSize',8);
            axis off;
            axes('Position',[0.01,0.02,0.01,0.01]); %for text
            text(0,0,file_string_cat,'Interpreter','none','FontSize',8);
            axis off;
        end
    end %if ismember
end %iplot
opts_vis_used=opts_vism;
return

function [connect_list,warn_msg]=psg_visualize_getconnect(connect_specs,nmult)
%check format of connect_specs and make a list of datasets to be connected
connect_list=zeros(0,2);
warn_msg=[];
if isempty(connect_specs)
    return
end
if (nmult==1)
    warn_msg=strvcat(warn_msg,'connect_specs ignored, only one dataset plotted');
    return
end
if_ok=1;
if isnumeric(connect_specs)
    if size(connect_list,2)~=2
        warn_msg=strvcat(warn_msg,'bad size for connect_specs; ignored');
        if_ok=0;        
    end
    if any(connect_specs(:)>nmult) | any(connect_specs(:)<1)
        warn_msg=strvcat(warn_msg,'bad entries in connect_specs; ignored');
        if_ok=0;
    end
    if (if_ok)
        connect_list=connect_specs;
    end
else
    switch connect_specs
        case 'star'
            connect_list=[ones(1,nmult-1);[2:nmult]]';
        case 'circuit'
            if nmult==2
                connect_list=[1 2];
            else
                connect_list=[[1:nmult-1 nmult];[2:nmult 1]]';
            end
        otherwise
           warn_msg=strvcat(warn_msg,'bad non-numeric specifier for connect_specs; ignored');
    end
end
return

function warn_msg=psg_visualize_check(sa)
warn_msg=[];
nmult=length(sa);
sa1=sa{1};
for im=1:nmult
    if size(sa{1}.typenames,1)~=size(sa{im}.typenames,1)
        warn_msg=strvcat(warn_msg,sprintf('sa structures 1 and %2.0f differ in length of typenames',im));
    else
        if any(strcmp(sa{1}.typenames,sa{im}.typenames)==0)
            warn_msg=strvcat(warn_msg,sprintf('sa structures 1 and %2.0f differ in typenames',im));
        end
    end
    if isfield(sa{1},'spec_labels') & isfield(sa{im},'spec_labels')
        if size(sa{1}.spec_labels,1) ~=size(sa{im}.spec_labels,1)
            warn_msg=strvcat(warn_msg,sprintf('sa structures 1 and %2.0f in differ length of spec_labels',im));
        else
            if any(strcmp(sa{1}.spec_labels,sa{im}.spec_labels)==0)
                warn_msg=strvcat(warn_msg,sprintf('sa structures 1 and %2.0f differ in spec_labels',im));
            end
        end
    end
end
return

function opts_new=psg_visualize_range(opts,opu) %combine plot_range
if ~isempty(opu.plot_range)
    opts_new=opts;
    opts_new.plot_range(1,:)=min([opts.plot_range(1,:);opu.plot_range(1,:)],[],1);
    opts_new.plot_range(2,:)=max([opts.plot_range(2,:);opu.plot_range(2,:)],[],1);
return
end


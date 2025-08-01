function opts_used=psg_plotcoords(coords,dim_select,sa,rays,opts)
% opts_used=psg_plotcoords(coords,dim_select,sa,rays,opts) plots psg coordinates with nice coloring
%
% * "origin" refers to the coordinates of the random texture, which may not be zero
% * can use tags in individual plots to control legends
%
%  behavior with if_use_rays=1 (default) connects points along lines, and provides a legend with the ray endpoints
%   if_use_rays changes behavior:
%      if_use_rays=0, ray lines are not plotted; data points are not connected on ray lines; ray colors are not used, and legend does not refer to rays
%    BUT
%   individual points can be labeled, (see opts.label*), and points can be colored according to datasets (opts.color_norays), and segments connecting across datasets
%   can be colored (opts.color_connect_sets_norays)
% 
% coords: If 2d: [nstims nd] are the coordinates for each stimulus type in a single dataset
%    If 3d: [nstims nd nsets] are corresponding coordinates across
%    datasets; points plotted and then connected as specified in opts.connect_list
% dim_select: a vector of coordinates to plot, each must be in [1 nd]
%   coordinates are selected after the optional transformation
% sa: the setup structure returned from psg_readcoord_data
% rays: the ray structure returned by psg_findrays
% opts: options
%   opts.axis_handle: handle to the axis, new axis opened if empty or omitted
%   opts.line_width: line width; defaults to 1; 0 to omit lines
%   opts.line_type: line type (style)
%   opts.line_type_connect_neg: line type (style) to connect points between datasets on negative rays
%   opts.xform_offset: [1 nd], vector to translate coordinates, defaults to zeros(1,nd)
%   opts.xform_mult: [nd nd], vector to multiply coordinates
%     coordinates are plotted as (coords-xform_offset)*xform_mult
%   opts.connect_list: two-column array of datasets to be paired
%     defaults to empty, entries must be in range [1 nsets] typically generated by psg_visualize
%     if non-empty, only connecting points will be plotted.
%   opts.color_norays_connect_mode: how colors are assigned to a connection between datasets
%     0: use color_connect_sets_norays, 1: use color of first set, 2: use color of second set, 3: split the connection and use both to halfway point
%   opts.color_norays: an rgb triple, or a color letter, for points, when rays are not used
%   opts.tag_text: optional tag text to put in the 'tag' field of plots (for later choosing the legends)
%   opts.if_rings: plot rings between points at equal distances from origin
%     (see psg_findrays); characteristics controlled by opts.line_width_ring,
%     opts.line_type_ring, opts.color_ring
%   opts.tet_vertices: size is (4,3); coordinates of vertices of tetrahedron in 3-space, should all be at unit distance from origin
%     (this and all other opts.tet* ignored unless nd=4)
%   opts.tet_signs; 4-vector of +/- 1, indicating sign-flips for coordinates that point to vertices of tetrahedra
%   opts.tet_show: if a list, tetrahedra are drawn with vertices at these distances from origin; Inf (default) plots at largest value within box 
%   opts.tet_view: view orientation (defaults to [10 68], shows all vertices well)
%   opts.tet_color: color of edges of tetrahedron
%   opts.tet_line_width: line width of edges of tetrahedron
%   opts.tet_line_type_side: line type of edges of tetrahedron
%   opts.tet_line_type_axis[_neg]: line type of axes pointing to vertices of tetrahedron [_neg: line type used if entry in tet_signs is -1]
%   opts.label_sets: which datasets (can be empty) to label individual points, used only if if_use_rays=0, defaults to 0
%   opts.label_list: cell array of labels or 'typenames' to use sa.typenames
%   opts.colors, opts.colors_anymatch, opts.symbs_anymatch: options for psg_typenames2colors
% opts_used: options used
%
%  19Dec22: Invoke psg_typenames2colors for ray colors and symbols
%  21Dec22: Add transformations
%  13Jan23: Allow for dim 3 of coords to be other sets to attach, add opts.tag_text
%  24Jan23: Add ring plotting
%  31Jan23: Add 4-d renderings
%  17Jun23: Add plotting of nearest neighbor pairs
%  28Jun23: Add failsafe if legend is empty (could happen if no rays are identified)
%  01Jul23: Add failsafe if rgbs are NaN
%  04Jul23: Use psg_spec2legend for legend labels; fix a matlab issue in legends by adding DisplayName property
%  24Jul23: Use point with largest multiplier for the legend
%  31Oct23: Add opts_used.plot_range (array of size [2 3]) to indicate range plotted
%  14Nov23: Bug fix, typo (psg_plotcoords23->psg_plotcoords_23)
%  15Nov23: Add color_nearest_nbr,noray_connect
%  28Nov23: Plot rays in increasing order of multipliers
%  19Feb24: Add opts.colors_anymatch, opts.symbs_anymatch
%  27Apr24: Options for if_use_rays=0: add labels via opts.label_[sets|list|font_size], change marker to marker_noray
%  26May24: With norays:  color to connect datasets (opts.color_connect_sets_norays) can be specified separately from origin color
%           and opts.color_norays_connect_mode: how segments between datasets are colored
%  25Oct24: label_list can be 'typenames' to use sa.typenames for labels
%  24Dec24: add color_norays_connect_mode=3: use two colors, each halfway to midpoint
% 
%  See also: PSG_READ_COORDDATA, PSG_FINDRAYS, PSG_DEFOPTS, PSG_VISUALIZE_DEMO, FILLDEFAULT,
%    PSG_TYPENAMES2COLORS, PSG_VISUALIZE, PSG_SPEC2LEGEND, PSG_LEGEND_KEEP.
%
if (nargin<5)
    opts=struct;
end
opts=filldefault(opts,'axis_handle',[]);
opts=filldefault(opts,'line_width',1); %0 to omit lines
opts=filldefault(opts,'line_width_ring',1);
opts=filldefault(opts,'line_type',[]); %line type
opts=filldefault(opts,'line_type_connect_neg','--'); %line type for negative directions for connections
opts=filldefault(opts,'line_type_ring',':');
opts=filldefault(opts,'marker_sign','*+'); %symbols for negative and postive values on rays
opts=filldefault(opts,'marker_origin','o'); %symbol for origin
opts=filldefault(opts,'marker_noray','.'); %symbol if no ray
opts=filldefault(opts,'marker_size',8); %marker size
%opts=filldefault(opts,'color_rays',{[.3 .3 .3],[1 0 0],[0 .7 0],[0 0 1]}); %colors to cycle through for each ray, supplanted by psg_typenames2colors
opts=filldefault(opts,'color_origin',[0 0 0]); %color used for origin
opts=filldefault(opts,'color_nearest_nbr',[0 0 0]); %color for interconnections of nearest-neighbor points in same datset
opts=filldefault(opts,'color_ring',[0 0 0]);
opts=filldefault(opts,'noray_connect',1); %connect points not on rays (ray indicator=NaN) to each other
%
opts=filldefault(opts,'if_origin_on_rays',1); %1 to include origin on rays 
opts=filldefault(opts,'axis_label_prefix','dim'); % prefix for axis label
opts=filldefault(opts,'if_use_rays',1); %1 to use ray structure, 0 to ignore
opts=filldefault(opts,'lims',[]); %2-element row vector, or, [nd 2], indicating axis ranges
opts=filldefault(opts,'if_grid',1); %1 to display grid, 0 for no grid
opts=filldefault(opts,'if_just_data',0); %1 omits plotting of axis labels, setting limits
opts=filldefault(opts,'if_legend',1);% 1 for a legend
opts=filldefault(opts,'connect_list',zeros(0,2)); %which pairs of datasets to connect
opts=filldefault(opts,'tag_text','');
opts=filldefault(opts,'if_rings',0);
opts=filldefault(opts,'if_nearest_neighbor',-1); %connections within a dataset: 0 not to connect nearest neighbor, -1 to plot if any points unassigned, 1 to plot always
%
opts=filldefault(opts,'tet_signs',[1 1 1 1]);
opts=filldefault(opts,'tet_vertices',[1  1  1;1 -1 -1;-1 1  -1;-1 -1 1]/sqrt(3)); 
opts=filldefault(opts,'tet_show',Inf);
opts=filldefault(opts,'tet_color','k');
opts=filldefault(opts,'tet_line_width',1);
opts=filldefault(opts,'tet_line_type_side','-');
opts=filldefault(opts,'tet_line_type_axis','-');
opts=filldefault(opts,'tet_line_type_axis_neg','--');
opts=filldefault(opts,'tet_view',[10 68]);
%
%these are only active if if_use_rays=0
%
opts=filldefault(opts,'label_sets',0);
opts=filldefault(opts,'label_list',{' '});
opts=filldefault(opts,'label_font_size',8);
opts=filldefault(opts,'connect_only',0); %set to only show connections
opts=filldefault(opts,'color_norays',[0 0 0]);
opts=filldefault(opts,'color_connect_sets_norays',{opts.color_origin}); %color to use to connect corresponding points in different datasets, if no rays, must be cell
opts=filldefault(opts,'color_norays_connect_mode',2); % opts_mult.color_norays_connect_mode: how segments between datasets are colored, if no rays
%
if strcmp(opts.label_list,'typenames')
    opts.label_list=sa.typenames;
end
%
% 0->use origin color, 1->use color of dataset in connect_specs(:,1), 2[default]->use color of dataset in connect_specs(:,2)
%
opts.plot_range=[];
%
nsets=size(coords,3);
nconnect=size(opts.connect_list,1); %if nconnect>0, just plot connections, otherwise just plot individual dataset
%
nd=size(coords,2);
opts=filldefault(opts,'xform_offset',zeros(1,nd));
opts=filldefault(opts,'xform_mult',eye(nd));
%
%set up choices for psg_typenames2colors
%
opts_tn2c=struct();
opts_tn2c.symbs.z=opts.marker_origin;
opts_tn2c.symbs.m=opts.marker_sign(1);
opts_tn2c.symbs.p=opts.marker_sign(2);
opts_tn2c.symbs_nomatch=opts.marker_noray;
opts_tn2c.colors_nomatch=opts.color_origin;
if isfield(opts,'colors')
    opts_tn2c.colors=opts.colors;
end
if isfield(opts,'colors_anymatch')
    if ischar('opts.colors_anymatch') 
        opts.colors_anymatch=get(line('color',opts.colors_anymatch,'Visible','off'),'color'); %idea from StackOverflow
    end
    opts_tn2c.colors_anymatch=opts.colors_anymatch;
end
if isfield(opts,'symbs_anymatch')
    opts_tn2c.symbs_anymatch=opts.symbs_anymatch;
end
%
if isempty(opts.axis_handle)
    opts.axis_handle=gca;
end
nstims=size(coords,1);
ndplot=length(dim_select);
hl=[];
ht=[];
if (ndplot==2) | (ndplot==3) | (ndplot==4)
    if (opts.if_use_rays==0)
        if (nconnect==0) | opts.connect_only==0
            [hp,hps,opts]=psg_plotcoords_23(coords,dim_select,opts.marker_noray,opts);
            if ~isempty(hp)
                for ih=1:length(hps)
                    set(hps{ih},'Color',opts.color_norays);
                    set(hps{ih},'MarkerSize',opts.marker_size);
                    set(hps{ih},'Tag',sprintf('%s set %2.0f data',opts.tag_text,ih));
                end
                hl=hp;
                ht='data';
            end
        end
        if (nconnect>0)
            for ic=1:size(opts.connect_list,1)
                for ip=1:size(coords,1)
                    coords_connect=[coords(ip,:,opts.connect_list(ic,1));coords(ip,:,opts.connect_list(ic,2))];
                    if opts.color_norays_connect_mode<3
                        [hc,hcs,opts]=psg_plotcoords_23(coords_connect,dim_select,[],setfield(opts,'label_sets',0)); %plot with no symbol and no label
                         if ~isempty(hc)
                             for ih=1:length(hc)
                                if opts.color_norays_connect_mode==0
                                    set(hcs{ih},'Color',opts.color_connect_sets_norays{1});
                                 else
                                    color_index=mod(opts.connect_list(ic,opts.color_norays_connect_mode)-1,length(opts.color_connect_sets_norays))+1; 
                                    set(hcs{ih},'Color',opts.color_connect_sets_norays{color_index});
                                end
                                set(hcs{ih},'LineWidth',opts.line_width);
                                set(hcs{ih},'Tag',sprintf('%s set %2.0f connect %3.0f point %3.0f',opts.tag_text,ih,ic,ip));
                             end
                         end
                    else
                        %split the segment
                        coords_midpoint=mean(coords_connect,1);
                        for im=1:2
                            [hc,hcs,opts]=psg_plotcoords_23([coords_connect(im,:);coords_midpoint],dim_select,[],setfield(opts,'label_sets',0)); %plot with no symbol and no label
                            if ~isempty(hc)
                                for ih=1:length(hc)
                                    color_index=mod(opts.connect_list(ic,im)-1,length(opts.color_connect_sets_norays))+1; 
                                    set(hcs{ih},'Color',opts.color_connect_sets_norays{color_index});
                                    set(hcs{ih},'LineWidth',opts.line_width);
                                    set(hcs{ih},'Tag',sprintf('%s set %2.0f connect %3.0f part %1.0f point %3.0f',opts.tag_text,ih,ic,im,ip));
                                end
                            end
                        end %im
                    end %color_norays_connect_mode=3?
                end %ip
            end %ic
        end %nconnect
    else %use ray structure
        %plot origin
        points=find(rays.whichray==0);
        if (nconnect==0) %plot each dataset
            [hp,hps,opts]=psg_plotcoords_23(coords(points,:,:),dim_select,opts.marker_origin,opts);
            if ~isempty(hp)
                for ih=1:length(hps)
                    set(hps{ih},'Color',opts.color_origin);
                    set(hps{ih},'MarkerSize',opts.marker_size);
                    set(hps{ih},'Tag',sprintf('%s set %2.0f data origin',opts.tag_text,ih));
                end
            end
        else %connect origins
            for ic=1:size(opts.connect_list,1)
                for ip=1:length(points)
                    coords_connect=[coords(points(ip),:,opts.connect_list(ic,1));coords(points(ip),:,opts.connect_list(ic,2))];
                    [hc,hcs,opts]=psg_plotcoords_23(coords_connect,dim_select,[],opts); %plot with no symbol
                    if ~isempty(hc)
                        for ih=1:length(hc)
                            set(hcs{ih},'Color',opts.color_origin);
                            set(hcs{ih},'Tag',sprintf('%s set %2.0f connect %3.0f origin',opts.tag_text,ih,ic));
                        end
                    end
                end %ip
            end %ic
        end %nconnect
        %plot points not on rays, without connecting datasets
        points=find(isnan(rays.whichray));
        if opts.noray_connect
            [hp,hps,opts]=psg_plotcoords_23(coords(points,:,:),dim_select,opts.marker_noray,opts);
            if ~isempty(hp)
                for ih=1:length(hps)
                    set(hps{ih},'Color',opts.color_origin);
                    set(hps{ih},'MarkerSize',opts.marker_size);
                    set(hps{ih},'Tag',sprintf('%s set %2.0f data noray',opts.tag_text,ih));
                end
            end
        end
        %if requested, or if some points are not assigned to rays, connect nearest-neighbor points within datasets
        if_nearest_neighbor=opts.if_nearest_neighbor;
        if if_nearest_neighbor==-1
            if_nearest_neighbor=any(isnan(rays.whichray));
        end
        if if_nearest_neighbor
            for ipair=1:rays.npairs
                points=rays.pairs(ipair,:);
                [hp,hps,opts]=psg_plotcoords_23(coords(points,:,:),dim_select,'-',opts);
                if ~isempty(hp)
                    for ih=1:length(hps)
                        set(hps{ih},'Color',opts.color_nearest_nbr);
                        set(hps{ih},'MarkerSize',opts.marker_size);
                        set(hps{ih},'Tag',sprintf('%s set %2.0f data pair',opts.tag_text,ih,points));
                    end
                end
            end %ipair
        end
        %set up legends
        hl=cell(0);
        ht=[];
        %plot each ray in each direction, first as points and then connect
        for iray=1:rays.nrays
            allpoints_no=find(rays.whichray==iray);
            if (opts.if_origin_on_rays)
                allpoints=union(allpoints_no,find(rays.whichray==0)); %all points on the ray, including origin
            else
                allpoints=allpoints_no;
            end
            allpoints=allpoints(:);
%           raycolor=opts.color_rays{1+mod(iray-1,length(opts.color_rays))};
            rgbs=zeros(2,3); %neg and pos, typically identical
            for isign=1:2 %1->neg, 2->pos
                if (isign==1)
                    sign_sel=find(rays.mult<0);
                    isign_val=-1;
                else
                    sign_sel=find(rays.mult>0);
                    isign_val=+1;
                end
                points=intersect(allpoints,sign_sel); %points to plot
                maxmult_index=min(find(abs(rays.mult(points))==max(abs(rays.mult(points))))); %24Jul23
                [rgbs(isign,:),symb]=psg_typenames2colors(sa.typenames(max(intersect(allpoints_no,sign_sel))),opts_tn2c);
                if (nconnect==0)
                    %need to plot the ray in order of increasing abs value of multipliers (28Nov23)
                    [mo_sorted,m]=sort(isign_val*rays.mult(points));
                    [hp,hps,opts]=psg_plotcoords_23(coords(points(m),:,:),dim_select,symb,opts);
                    if ~isempty(hp)
                        for ih=1:length(hps)
                            set(hps{ih},'Color',rgbs(isign,:));
                            set(hps{ih},'MarkerSize',opts.marker_size);
                            set(hps{ih},'Tag',sprintf('%s set %2.0f ray %2.0f signed %2.0f',opts.tag_text,ih,iray,isign));
                            set(hps{ih},'DisplayName',psg_spec2legend(sa,points(maxmult_index),[])); %04Jul23
                        end
                        hl=[hl;hp];
                        ht=strvcat(ht,psg_spec2legend(sa,points(maxmult_index),[])); %04Jul23
                    end
                else %connect corresponding points
                    for ic=1:size(opts.connect_list,1)
                        for ip=1:length(points)
                            coords_connect=[coords(points(ip),:,opts.connect_list(ic,1));coords(points(ip),:,opts.connect_list(ic,2))];
                            if (isign==1)
                                line_type=opts.line_type_connect_neg;
                            else
                                line_type=opts.line_type;
                            end
                            [hc,hcs,opts]=psg_plotcoords_23(coords_connect,dim_select,line_type,opts); %plot with no symbol
                            if ~isempty(hc)
                                for ih=1:length(hc)
                                    set(hcs{ih},'LineWidth',opts.line_width);
                                    set(hcs{ih},'Color',rgbs(isign,:));
                                    set(hcs{ih},'Tag',sprintf('%s set %2.0f ray %2.0f signed %2.0f connect %3.0f point %3.0f',opts.tag_text,ih,iray,isign,ic,ip));
                                    set(hps{ih},'DisplayName',psg_spec2legend(sa,points(maxmult_index),[])); %04Jul23
                                end
                                if (ip==maxmult_index)
                                    hl=[hl,hc];
                                    ht=strvcat(ht,psg_spec2legend(sa,points(maxmult_index),[])); %04Jul23
                                end
                            end
                        end %ip
                    end %ic
                end %nconnect
            end %isign
            if (opts.line_width>0) & (nconnect==0) %connect points within a single dataset only along rays
                mults=rays.mult(allpoints);
                [mults_sorted,sort_order]=sort(mults);
                [hp,hps,opts]=psg_plotcoords_23(coords(allpoints(sort_order),:,:),dim_select,opts.line_type,opts);
                if ~isempty(hp)
                    rgb_nonan=find(~any(isnan(rgbs),2));
                    if ~isempty(rgb_nonan)
                        mean_rgb=mean(rgbs(rgb_nonan,:),1);
                    else
                        mean_rgb=zeros(1,3);
                    end
                    for ih=1:length(hps)
                        set(hps{ih},'LineWidth',opts.line_width);
                        set(hps{ih},'Color',mean_rgb);
                        set(hps{ih},'Tag',sprintf('%s set %2.0f ray %2.0f nosign',opts.tag_text,ih,iray));
                    end
                end
            end %line width
        end %iray
     end %use rays
     if opts.if_rings
         nrings=length(rays.rings);
         for iring=1:nrings
             coord_ptrs=rays.rings{iring}.coord_ptrs;
             coord_cycle=[coord_ptrs coord_ptrs(1)]; %make it a cycle
             [hp,hps,opts]=psg_plotcoords_23(coords(coord_cycle,:,:),dim_select,opts.line_type_ring,opts);
             if ~isempty(hp)
                 set(hps{ih},'LineWidth',opts.line_width_ring);
                 set(hps{ih},'Color',opts.color_ring);
                 set(hps{ih},'Tag',sprintf('%s set %2.0f ring %2.0f',opts.tag_text,ih,iring));                
             end
         end
     end
     %
     %special formatting 4d plot
     %
     if (ndplot==4) %draw a tetrahedron
        if ~opts.if_just_data
            boxmax=max(abs([get(gca,'XLim'),get(gca,'YLim'),get(gca,'ZLim')])); 
            edges=nchoosek([1:ndplot],2);
            %draw and label axes
            opts.plot_range=zeros(2,3);
            for iaxis=1:ndplot
                axis_end=sqrt(3)*boxmax*opts.tet_vertices(iaxis,:);
                if opts.tet_signs(iaxis)<0
                    sign_char='-';
                    line_type=opts.tet_line_type_axis_neg;
                else
                    sign_char='+';
                    line_type=opts.tet_line_type_axis;
                end
                axis_label=sprintf('%s %1.0f (%s)',opts.axis_label_prefix,dim_select(iaxis),sign_char);
                htet=plot3([0 axis_end(1)],[0 axis_end(2)],[0 axis_end(3)],cat(2,'k',line_type));
                opts.plot_range(1,:)=min([opts.plot_range(1,:);axis_end(:)'],[],1);
                opts.plot_range(2,:)=max([opts.plot_range(2,:);axis_end(:)'],[],1);               
                set(htet,'Tag',sprintf('tet axis %1.0f %s',iaxis,axis_label));
                text(axis_end(1),axis_end(2),axis_end(3),axis_label);
            end %iaxis
            for itet=1:length(opts.tet_show)
                tsize=opts.tet_show(itet);
                if tsize==Inf
                    tsize=boxmax*sqrt(3);
                end
                for iedge=1:size(edges,1)
                    htet=plot3(tsize*opts.tet_vertices(edges(iedge,:),1),tsize*opts.tet_vertices(edges(iedge,:),2),tsize*opts.tet_vertices(edges(iedge,:),3),...
                        cat(2,'k',opts.tet_line_type_side));
                    set(htet,'Color',opts.tet_color);
                    set(htet,'LineWidth',opts.tet_line_width);
                    set(htet,'Tag',sprintf('tet size %5.2f vertices %1.0f %1.0f',tsize,edges(iedge,:)));
                end %iedge
            end %itet
        end %if_just_data
        set(gca,'View',opts.tet_view);
        axis off
    end %nd=4
    %
    if ~opts.if_just_data
        if (opts.if_legend)
            if ~isempty(hl) & ~isempty(ht)
                hleg=legend(hl,ht);
                set(hleg,'FontSize',7);
%               set(hleg,'String',ht);
            end
        end
        if (ndplot<=3)
            xlabel(sprintf('%s %1.0f',opts.axis_label_prefix,dim_select(1)));
            ylabel(sprintf('%s %1.0f',opts.axis_label_prefix,dim_select(2)));
        end
        if ~isempty(opts.lims)
            lims=opts.lims;
            if size(lims,1)==1
                lims=repmat(lims,ndplot,1);
            end
            set(gca,'XLim',lims(1,:));
            set(gca,'YLim',lims(2,:));
            if (ndplot==3)
                set(gca,'ZLim',lims(3,:));
            end
        end
        if (ndplot==3)
            zlabel(sprintf('%s %1.0f',opts.axis_label_prefix,dim_select(3)));
            axis vis3d;
        else
            if ~isempty(opts.lims)
                axis square;
            else
                axis equal;
            end
        end
        if opts.if_grid
            grid on;
        else
            grid off;
        end
    end %if_just_data
    hleg_exist=get(gca,'Legend'); %re-create the legend to prevent extra entries
    if ~isempty(hleg_exist) & (nconnect==0)
        if ~isempty(hl) & ~isempty(ht)
            legend(hl,ht);
        end
    end
end
opts_used=opts;
return

function [hp,hps,opts_used]=psg_plotcoords_23(coords_untrans,dim_select,symb,opts)
%plot points or lines in 2 or 3 dimensions. 
%line type taken from symb
%line width set to 1
%color set to black
%offset and multiplier taken from opts.xform_offset, opts.xform_mult
%
%for 4-dimensions, the coords are dotted with four corners of a
%tetrahedron (chosen as 4 non-adjacent corners of a cube):
%
%
nconds=size(coords_untrans,3);
nd=length(dim_select);
hp=[];
hps=cell(1,nconds);
if isempty(coords_untrans)
    opts_used=opts;
    return
end
if isempty(opts.plot_range)
    opts.plot_range=repmat([Inf;-Inf],1,min(nd,3));
end
for icond=1:nconds
    %transform as requested
    coords_trans=(coords_untrans(:,:,icond)-repmat(opts.xform_offset,size(coords_untrans,1),1))*opts.xform_mult;
    %plot selected coordinates
    coords=coords_trans(:,dim_select);
    %2 or 3d plot, with symbol specified
    %4d: project onto tetrahedron
    if (nd==4)
        coords=coords.*repmat(opts.tet_signs,size(coords,1),1)*opts.tet_vertices;
    end
    if size(coords,1)==0
        hp=[];
    elseif nd==2
        hp=plot(coords(:,1),coords(:,2),cat(2,'k',symb));
        hps{icond}=hp;
        hold on;
        if opts.if_use_rays==0 & ismember(icond,opts.label_sets)
            for ipt=1:size(coords,1)
                label_ptr=mod(ipt-1,length(opts.label_list))+1;
                text(coords(ipt,1),coords(ipt,2),opts.label_list{label_ptr},'FontSize',opts.label_font_size);
            end
        end
        opts.plot_range(1,:)=min([opts.plot_range(1,:);min(coords,[],1)],[],1);
        opts.plot_range(2,:)=max([opts.plot_range(2,:);max(coords,[],1)],[],1);
    elseif nd==3 | nd==4
        hp=plot3(coords(:,1),coords(:,2),coords(:,3),cat(2,'k',symb));
        hps{icond}=hp;
        hold on;
        if opts.if_use_rays==0 & ismember(icond,opts.label_sets)
            for ipt=1:size(coords,1)
                label_ptr=mod(ipt-1,length(opts.label_list))+1;
                text(coords(ipt,1),coords(ipt,2),coords(ipt,3),opts.label_list{label_ptr},'FontSize',opts.label_font_size);
            end
        end
        opts.plot_range(1,:)=min([opts.plot_range(1,:);min(coords,[],1)],[],1);
        opts.plot_range(2,:)=max([opts.plot_range(2,:);max(coords,[],1)],[],1);
    end
end
opts_used=opts;
return

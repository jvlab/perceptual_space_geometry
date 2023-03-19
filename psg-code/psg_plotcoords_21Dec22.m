function opts_used=psg_plotcoords(coords,dim_select,sa,rays,opts)
% opts_used=psg_plotcoords(coords,dim_select,sa,rays,opts) plots psg
% coordinates with nice coloring
%
% note that "origin" refers to the coordinates of the random texture, which
% may not be zero
%
% coords: [nstims nd] are the coordinates for each stimulus type
% dim_select: a vector of coordinates to plot, each must be in [1 nd]
%   coordinates are selected after the optional transformation
% sa: the setup structure returned from psg_readcoord_data
% rays: the ray structure returned by psg_findrays
% opts: options
%   opts.axis_handle: handle to the axis, new axis opened if empty or omitted
%   opts.xform_offset: [1 nd], vector to translate coordinates, defaults to zeros(1,nd)
%   opts.xform_mult: [nd nd], vector to multiply coordinates
%     coordinates are plotted as (coords-xform_offset)*xform_mult
%
% opts_used: options used
%
%  19Dec22: Invoke psg_typenames2colors for ray colors and symbols
%  21Dec22: Add transformations
%
%  See also: PSG_READ_COORDDATA, PSG_FINDRAYS, PSG_DEFOPTS, PSG_VISUALIZE_DEMO, FILLDEFAULT,
%    PSG_TYPENAMES2COLORS.
%
if (nargin<5)
    opts=struct;
end
opts=filldefault(opts,'axis_handle',[]);
opts=filldefault(opts,'line_width',1); %0 to omit lines
opts=filldefault(opts,'line_type',[]); %line type
opts=filldefault(opts,'marker_sign','*+'); %symbols for negative and postive values on rays
opts=filldefault(opts,'marker_origin','o'); %symbol for origin
opts=filldefault(opts,'marker_noray','.'); %symbol if no ray
opts=filldefault(opts,'marker_size',8); %marker size
%opts=filldefault(opts,'color_rays',{[.3 .3 .3],[1 0 0],[0 .7 0],[0 0 1]}); %colors to cycle through for each ray, supplanted by psg_typenames2colors
opts=filldefault(opts,'color_origin',[0 0 0]);
opts=filldefault(opts,'color_norays',[0 0 0]);
opts=filldefault(opts,'if_origin_on_rays',1); %1 to include origin on rays 
opts=filldefault(opts,'axis_label_prefix','dim'); % prefix for axis label
opts=filldefault(opts,'if_use_rays',1); %1 to use ray structure, 0 to ignore
opts=filldefault(opts,'lims',[]); %2-element row vector, or, [nd 2], indicating axis ranges
opts=filldefault(opts,'if_grid',1); %1 to display grid, 0 for no grid
opts=filldefault(opts,'if_just_data',0); %1 omits plotting of axis labels, setting limits
opts=filldefault(opts,'if_legend',1);% 1 for a legend
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
%
if isempty(opts.axis_handle)
    opts.axis_handle=gca;
end
nstims=size(coords,1);
ndplot=length(dim_select);
if (ndplot==2) | (ndplot==3)
    if (opts.if_use_rays==0)
        hp=psg_plotcoords_23(coords,dim_select,opts.marker_origin,opts);
        if ~isempty(hp)
            set(hp,'Color',opts.color_norays);
            set(hp,'MarkerSize',opts.marker_size);
            hl=hp;
            ht='data';
        end
    else %use ray structure
        %plot origin
        points=find(rays.whichray==0);
        hp=psg_plotcoords_23(coords(points,:),dim_select,opts.marker_origin,opts);
        if ~isempty(hp)
            set(hp,'Color',opts.color_origin);
            set(hp,'MarkerSize',opts.marker_size);
        end
        %set up legends
        hl=cell(0);
        ht=[];
        %plot each ray in each direction
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
                else
                    sign_sel=find(rays.mult>0);
                end
                points=intersect(allpoints,sign_sel); %points to plot
                [rgbs(isign,:),symb]=psg_typenames2colors(sa.typenames(max(intersect(allpoints_no,sign_sel))),opts_tn2c);
                hp=psg_plotcoords_23(coords(points,:),dim_select,symb,opts);
                if ~isempty(hp)
                    set(hp,'Color',rgbs(isign,:));
                    set(hp,'MarkerSize',opts.marker_size);
                    hl=[hl;hp];
                    ht=strvcat(ht,sa.typenames{points(end)});
                end               
            end
            if (opts.line_width>0)
                mults=rays.mult(allpoints);
                [mults_sorted,sort_order]=sort(mults);
                hp=psg_plotcoords_23(coords(allpoints(sort_order),:),dim_select,opts.line_type,opts);
                if ~isempty(hp)
                    set(hp,'LineWidth',opts.line_width);
                    set(hp,'Color',mean(rgbs));
                end
            end %isign
        end %iray
     end %use rays
    if ~opts.if_just_data
        if (opts.if_legend)
            hleg=legend(hl,ht);
            set(hleg,'FontSize',7);
%            set(hleg,'String',ht);
        end
        xlabel(sprintf('%s%1.0f',opts.axis_label_prefix,dim_select(1)));
        ylabel(sprintf('%s%1.0f',opts.axis_label_prefix,dim_select(2)));
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
            zlabel(sprintf('%s%1.0f',opts.axis_label_prefix,dim_select(3)));
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
    if ~isempty(hleg_exist)
        legend(hl,ht);
    end
end
opts_used=opts;
return

function hp=psg_plotcoords_23(coords_untrans,dim_select,symb,opts)
%transform as requested
coords_trans=(coords_untrans-repmat(opts.xform_offset,size(coords_untrans,1),1))*opts.xform_mult;
%plot selected coordinates
coords=coords_trans(:,dim_select);
%2 or 3d plot, with symbol specified

if size(coords,1)==0
    hp=[];
    return
elseif size(coords,2)==2
    hp=plot(coords(:,1),coords(:,2),cat(2,'k',symb));
    hold on;
elseif size(coords,2)==3
    hp=plot3(coords(:,1),coords(:,2),coords(:,3),cat(2,'k',symb));
    hold on;
end
return

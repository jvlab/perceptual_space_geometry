function [hplanes,ou]=btc_soid_plot(edirs,opts)
% [hplanes,ou]=btc_soid_plot(edirs,opts) plots directions, conditions,
% and data for a btc eperiment
%
%   March 2017:  added ability to plot mtc data
%    if edir.ifalloy is not present: mode=0, binary (btc) data assumed
%    if edir.ifalloy is present: 3 gray levels, assumed
%       mode 1: edir.ifalloy=0, "mixed" correlations, alloy transformation not applied
%       mode 2: edir.ifalloy=1, alloy transformation and labelling applied
%
% edirs:  a structure with fields for each plane ("xy")
%   indicating the experimental directions
%   key fields are edirs.yx.plane, where yx is
%      for binary textures: a 2-char string, such as 'ag'
%      for gray-level textures: for single planes (alloys plots), strings such as A1G, A1Gmerge, AB11B, ABC111T[full],ABCD1111A
%                               for mixed planes, a 3-char string such as YGN, or
% opts: plot options
%    opts.ncircle: circle quality, default to 48
%    opts.radii: circle radii, defaults to steps of 1/4 up to value specified in opts.lims opts.lims_plane
%    opts.fontsize: font size, defaults to 8
%    opts.lims: limits, defaults to 1, overridden by opts_lims_plane
%    opts.lims_plane: a structure with field names a subset of those of edirs, and values that override opts.lims in that plane
%        defaults to struct()
%    opts.datafield: strvcat of data field(s) to plot; can be a strvcat to plot more than one variable
%       this must have two columns, first column is abscissa, second is ordinate
%    opts.marker: symbols to use, e.g., 'k.'; can be a strvcat to plot more than one variable
%    opts.rgbcolor:  colors to use, 3 columns, can have many rows to plot more than one variable
%        if NaN (default), then the color associated with opts.marker ris used
%    opts.linewidth: line width for data, defaults to 1, can be vector to plot more than one variable (added Aug 2017)
%    opts.cyclic:  set to 1 to plot values as a closed contour, defaults to 0.  Overridden by edir.cyclic if present
%    opts.tstring: suffix for title string
%    opts.xalloy_triangle_color: color alloy in triangle plot, [] to omit
%    opts.xalloy_label: 1 to label vertices in alloy plot
%    opts.plane_select: cell array or strvcat of names of planes to plot
%       if empty or omitted, all planes (i.e., each substructure) in edirs are plotted
%       if a plane listed in plane_select is not present in edirs, then a subplot is left blank
%    opts.plane_layout:  [nr,nc] -- number of rows and columns in layout, if empty or omitted, then 
%       set by nicesubp
%
% hplanes: pointers to the individual axes of the plots
% ou: options used
%
% (version of 24Sep17 frozen as btc_soid_plot_24Sep17)
% 11Dec21: add (via opts.lims_plane) option to set limits for each plane
%          change behavior of radii if it is not supplied
%          update documentation
%
%  See also:  BTC_SOID_TEST, BTC_PAIRSNEEDED, BTC_EDIRS, BTC_DEFINE, BTC_SOID_PLOTQF, MTC_LET2GBC,
%     MTC_XALLOY_MTX, ,BTC_SOID_STATS, BTC_SOID_DEMO, MTC_SOID_PLOT_DEMO, MGM_PRED_SUMM.
%
if (nargin<=1) opts=[]; end
opts=filldefault(opts,'fontsize',8);
opts=filldefault(opts,'ncircle',48); %0 for no circle
opts=filldefault(opts,'radii',[]); %empty for no circle
opts=filldefault(opts,'lims',1);
opts=filldefault(opts,'lims_plane',struct());
opts=filldefault(opts,'marker','k.');
opts=filldefault(opts,'rgbcolor',[NaN NaN NaN]);
opts=filldefault(opts,'linewidth',1);
opts=filldefault(opts,'cyclic',0);
opts=filldefault(opts,'tstring',' ');
opts=filldefault(opts,'xalloy_triangle_color','g'); %color alloy in triangle plot, [] to omit
opts=filldefault(opts,'xalloy_label',1); %1 to label vertices in alloy plot
opts=filldefault(opts,'plane_select',[]);
opts=filldefault(opts,'plane_layout',[]);
%
ou=opts;
%
if isempty(opts.plane_select)
    planes=char(fieldnames(edirs));
else
    planes=char(opts.plane_select);
end
nplanes=size(planes,1);
if isempty(opts.plane_layout)
    [nr,nc]=nicesubp(nplanes,0.7);
else
    nr=opts.plane_layout(1);
    nc=opts.plane_layout(2);
end
if (opts.ncircle>0)
    circ_x=cos([0:opts.ncircle]*2*pi/opts.ncircle);
    circ_y=sin([0:opts.ncircle]*2*pi/opts.ncircle);
end
hplanes=[];
for iplane=1:nplanes
    plane_name=deblank(planes(iplane,:));
    if isfield(edirs,plane_name)
        hplanes(iplane)=subplot(nr,nc,iplane,'FontSize',opts.fontsize);
        edir=getfield(edirs,plane_name);
        lims=opts.lims;
        if isfield(opts.lims_plane,plane_name)
            lims=opts.lims_plane.(plane_name);
        end
        if ~isfield(edir,'ifalloy')
            mode=0; %btc (binary) mode
            tstring=deblank(cat(2,planes(iplane,:),' (',edir.opts.layout,edir.plane,')  ',opts.tstring));
            xstring=edir.plane(2); %the x-y swap because of how planes are named (e.g., 'ag')
            ystring=edir.plane(1);
        else
            xstring=[]; %need to fix for alloys
            ystring=[];
            if (edir.ifalloy==0)
                mode=1; %mtc, not alloy
                if (edir.plane(1)=='Y')
                    [cgnp,cgnb,dnp,dnb]=mtc_let2gbc(edir.plane(2)); %plane name like "YGM", etc,
                    xstring=cat(2,edir.plane(2),': ',cgnp,' ',dnb);
                    [cgnp,cgnb,dnp,dnb]=mtc_let2gbc(edir.plane(3));
                    ystring=cat(2,edir.plane(3),': ',cgnp,' ',dnb);
                else
                    warning(sprintf('plane name (%s) not recognized',edir.plane));
                end
            else
                mode=2; %mtc, alloy
                [xalloy_mtx,triangle_raw]=mtc_xalloy_mtx;
            end
            BigConds_string=[];
            for ibc=1:length(edir.opts.bigconds_names)
                BigConds_string=cat(2,BigConds_string,edir.opts.bigconds_names{ibc},' ');
            end
            tstring=deblank(cat(2,planes(iplane,:),' (',deblank(BigConds_string),') ',opts.tstring));
        end
        nvars=size(opts.datafield,1);
        for ivar=1:nvars
            if isfield(edir,deblank(opts.datafield(ivar,:)))
                xy=getfield(edir,deblank(opts.datafield(ivar,:)));
                marker=deblank(opts.marker(min(ivar,size(opts.marker,1)),:));
                if isfield(edir,'cyclic')
                   cyclic=edir.cyclic;
                else
                    cyclic=opts.cyclic(min(ivar,length(opts.cyclic)));
                end
                linewidth=opts.linewidth(min(ivar,length(opts.linewidth)));
                indices=[1:size(xy,1)]';
                if (cyclic==1)
                    indices=[indices;1];
                end
                xyp=xy(indices,[1:2]);
                if (mode==2)
                    xyp=xyp*xalloy_mtx; %alloy transformation
                end
                hpp=plot(xyp(:,1),xyp(:,2),marker,'LineWidth',linewidth);
                rgbcolor=opts.rgbcolor(min(ivar,size(opts.rgbcolor,1)),:);
                if ~isnan(rgbcolor)
                    set(hpp,'color',rgbcolor);
                end
                hold on
            end %isfield
        end
        %special labeling if alloy transformation
        nt=3;
        if (mode==2)
            triangle_plot=triangle_raw*xalloy_mtx;
            if ~isempty(opts.xalloy_triangle_color)
                hp=plot(triangle_plot([1:nt 1],1),triangle_plot([1:nt 1],2));
                set(hp,'Color',opts.xalloy_triangle_color);
            end
            if (opts.xalloy_label==1)
                for it=1:nt
                    ht=text(triangle_plot(it,1),triangle_plot(it,2),0,sprintf('p=(%1.0f,%1.0f,%1.0f)',double([1:nt]==it)));
                    set(ht,'FontSize',opts.fontsize);
                end
            end
        end %end of mode 2 (alloy) labeling
        plot([-1 1]*lims,[0 0],'k');
        plot([0 0],[-1 1]*lims,'k');
        xlabel(xstring,'FontSize',opts.fontsize);
        ylabel(ystring,'FontSize',opts.fontsize);
        if (opts.ncircle>0)
            radii_list=opts.radii;
            if isempty(radii_list)
                radii_list=[0:1/4:lims];
            end
            for k=1:length(radii_list)
                plot(circ_x*radii_list(k),circ_y*radii_list(k),'k:');
            end
        end
        ht=title(tstring);
        set(ht,'Interpreter','none');
        set(gca,'XLim',[-1 1]*lims);
        set(gca,'YLim',[-1 1]*lims);
        if (mode<2)
            set(gca,'XTick',[-1 1]*lims);
            set(gca,'YTick',[-1 1]*lims);
        else
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
        end
        axis square;
    end %isfield
end
return


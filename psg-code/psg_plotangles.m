function opts_visang_used=psg_plotangles(angles,sa,rays,opts_visang)
% opts_visang_used=psg_plotangles(angles,sa,rays,opts_visang) shows how angles behave as function of dimension
%
% angles: angles to plot, typically returned by psg_rayangles,
%   each field is [isign, jsign point to neg and pos if if_bid=0, or are =1 if_bid=1]
%   dotprods: dot products, (iray,jray,isign,jsign)
%   cosangs:  cosine of angles (iray,jray,isign,jsign)
%   ang_degs: anlges in degress (iray,jray, isign, jsign)
% sa: setup structure, from psg_read_coorddata
% rays: a ray structure, typically from psg_findrays
% opts_visang: options
%   opts_visang.vis_string: 'raw fit'; how dimension label is formatted
%   opts_visang.file_string: file name string
%   opts_visang.dim_range: dimension range, defaults to [1 10]
%   opts_visang.var_plot: variable to plot ('ang_degs','cosangs','dotprods', defaults to 'ang_degs')
%
%  to do: legends
%
%   opts_visang_used: options used, and handle to figure
%
%   See also: PSG_FINDRAYS, PSG_RAYFIT, PSG_PLOTCOORDS, PSG_RAYANGLES, PSG_TYPENAMES2COLORS.
%
if (nargin<4)
    opts_visang=struct();
end
opts_visang=filldefault(opts_visang,'dim_range',[1 10]);
opts_visang=filldefault(opts_visang,'file_string',[]);
opts_visang=filldefault(opts_visang,'vis_string','raw fit');
opts_visang=filldefault(opts_visang,'var_plot','ang_degs');
opts_visang=filldefault(opts_visang,'marker_size',8); %marker size
opts_visang=filldefault(opts_visang,'marker_sign','*+'); %symbols for negative and postive values on rays
opts_visang=filldefault(opts_visang,'marker_origin','o'); %symbol for origin
opts_visang=filldefault(opts_visang,'marker_noray','.'); %symbol if no ray
opts_visang=filldefault(opts_visang,'if_legend',1);
%
%set up choices for psg_typenames2colors
%
opts_tn2c=struct();
opts_tn2c.symbs.z=opts_visang.marker_origin;
opts_tn2c.symbs.m=opts_visang.marker_sign(1);
opts_tn2c.symbs.p=opts_visang.marker_sign(2);
opts_tn2c.symbs_nomatch=opts_visang.marker_noray;
if isfield(opts_visang,'colors')
    opts_tn2c.colors=opts_visang.colors;
end
%
plotfmt.ang_degs.range=[0 180];
plotfmt.ang_degs.lines=90;
plotfmt.dotprods.range=[NaN NaN];
plotfmt.dotprods.lines=0;
plotfmt.cosangs.range=[-1 1];
plotfmt.cosangs.lines=0;
%
dims_avail=[];
dplot=[];
for idim=1:length(angles)
    if ~isempty(angles{idim})
        dplot=cat(5,dplot,angles{idim}.(opts_visang.var_plot));
        dims_avail=[dims_avail,idim];
    end
end
plotfmt_use=plotfmt.(opts_visang.var_plot);
if any(isnan(plotfmt_use.range))
    rmax=max(abs(dplot(:)));
    rmax_round=ceil(rmax/10^floor(log10(rmax)))*(10^floor(log10(rmax)));
    plotfmt_use.range=[-1 1]*rmax_round;
end
%    
np=size(dplot,3);
nrays=size(dplot,1);
if (np==1)
    fig_title=cat(2,opts_visang.file_string,' ',opts_visang.vis_string,': bid angles');
else
    fig_title=cat(2,opts_visang.file_string,' ',opts_visang.vis_string,': ray angles');
end
%
opts_visang.fig_handle=figure;
set(gcf,'Position',[50 50 1200 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name',fig_title);
%
%each row of the plot is one ray ignoring sign
%in each row, have np^2 plots
%
sign_select=cell(np,1);
if np==1
    sign_select{1}=find(rays.mult~=0);
    sign_char={' '};
else
    sign_select{1}=find(rays.mult<0);
    sign_select{2}=find(rays.mult>0);
    sign_char={'m','p'};
end

for iray=1:nrays
    for ip=1:np
        for jp=1:np
            ipjp=jp+(ip-1)*np;
            subplot(nrays,4,ipjp+(iray-1)*4);
            ipoints=intersect(find(rays.whichray==iray),sign_select{ip});
            ilab=sa.typenames{ipoints(end)};
            %set up legends
            hl=cell(0);
            ht=[];
            for jray=1:nrays
                jpoints=intersect(find(rays.whichray==jray),sign_select{jp});
                jlab=sa.typenames{jpoints(end)};
                [rgb,symb]=psg_typenames2colors(sa.typenames(jpoints(end)),opts_tn2c);
                hp=plot(dims_avail,squeeze(dplot(iray,jray,ip,jp,:)),cat(2,'k',symb));
                set(hp,'Color',rgb);
                set(hp,'MarkerSize',opts_visang.marker_size);
                hold on;
                hl=[hl;hp];
                ht=strvcat(ht,cat(2,sa.typenames{jpoints(1)},'-',sa.typenames{jpoints(end)}));
                hpline=plot(dims_avail,squeeze(dplot(iray,jray,ip,jp,:)),'k');
                set(hpline,'Color',rgb);
            end
            for k=1:length(plotfmt_use.lines)
                plot([-0.5 0.5]+opts_visang.dim_range,repmat(plotfmt_use.lines(k),1,2),'k:');
            end
            if (opts_visang.if_legend)
                hleg=legend(hl,ht);
                set(hleg,'FontSize',7);
%                set(hleg,'String',ht);
            end
            xlabel('dim');
            ylabel(opts_visang.var_plot,'Interpreter','none');
            set(gca,'XLim',[-0.5 0.5]+opts_visang.dim_range);
            set(gca,'XTick',dims_avail);
            set(gca,'YTick',sort([plotfmt_use.range,plotfmt_use.lines]));
            set(gca,'YLim',plotfmt_use.range);
            title(cat(2,ilab,' ',sign_char{ip},' ',sign_char{jp}));
        end %jp
    end %ip
end %iray
axes('Position',[0.01,0.02,0.01,0.01]); %for text
text(0,0,fig_title,'Interpreter','none');
axis off;
%
opts_visang_used=opts_visang;
return


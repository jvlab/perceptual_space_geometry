function opts_vismult_used=psg_plotmults(mults,sa,rays,opts_vismult)
% opts_vismult_used=psg_plotmults(mults,sa,rays,opts_vismult) shows how angles behave as function of dimension
%
% mults: multiples to plot, typically returned by psg_raymults,
%   each field is [isign, jsign point to neg and pos if if_bid=0, or are =1 if_bid=1]
%   dotprods: dot products, (iray,jray,isign,jsign)
%   cosangs:  cosine of angles (iray,jray,isign,jsign)
%   ang_degs: anlges in degress (iray,jray, isign, jsign)
% sa: setup structure, from psg_read_coorddata
% rays: a ray structure, typically from psg_findrays
% opts_vismult: options
%   opts_vismult.vis_string: 'raw fit'; how dimension label is formatted
%   opts_vismult.file_string: file name string
%   opts_vismult.dim_range: dimension range, defaults to [1 10]
%   opts_vismult.var_plot: variable to plot ('dist_gain','dist_raw','dist_endpts', defaults to 'dist_gain')
%
%  to do: legends
%
%   opts_vismult_used: options used, and handle to figure
%
% 17Jun23: protect from empty rays
% 04Jul23: use psg_spec2legend for legends
%
%   See also: PSG_FINDRAYS, PSG_RAYFIT, PSG_PLOTCOORDS, PSG_RAYMULTS, PSG_TYPENAMES2COLORS, PSG_PLOTANGLES.
%
if (nargin<4)
    opts_vismult=struct();
end
opts_vismult=filldefault(opts_vismult,'dim_range',[1 10]);
opts_vismult=filldefault(opts_vismult,'file_string',[]);
opts_vismult=filldefault(opts_vismult,'vis_string','raw fit');
opts_vismult=filldefault(opts_vismult,'var_plot','dist_gain');
opts_vismult=filldefault(opts_vismult,'marker_size',8); %marker size
opts_vismult=filldefault(opts_vismult,'marker_sign','*+'); %symbols for negative and postive values on rays
opts_vismult=filldefault(opts_vismult,'marker_origin','o'); %symbol for origin
opts_vismult=filldefault(opts_vismult,'marker_noray','.'); %symbol if no ray
opts_vismult=filldefault(opts_vismult,'if_legend',1);
%
%set up choices for psg_typenames2colors
%
opts_tn2c=struct();
opts_tn2c.symbs.z=opts_vismult.marker_origin;
opts_tn2c.symbs.m=opts_vismult.marker_sign(1);
opts_tn2c.symbs.p=opts_vismult.marker_sign(2);
opts_tn2c.symbs_nomatch=opts_vismult.marker_noray;
if isfield(opts_vismult,'colors')
    opts_tn2c.colors=opts_vismult.colors;
end
%
plotfmt.dist_gain.range=[NaN NaN];
plotfmt.dist_gain.lines=0;
plotfmt.dist_raw.range=[NaN NaN];
plotfmt.dist_raw.lines=0;
plotfmt.dist_endpts.range=[NaN NaN];
plotfmt.dist_endpts.lines=0;
%
dims_avail=[];
dplot=[];
for idim=1:length(mults)
    if ~isempty(mults{idim})
        dplot=cat(3,dplot,mults{idim}.(opts_vismult.var_plot));
        dims_avail=[dims_avail,idim];
    end
end
plotfmt_use=plotfmt.(opts_vismult.var_plot);
if any(isnan(plotfmt_use.range))
    rmax=max(abs(dplot(:)));
    rmax_round=ceil(rmax/10^floor(log10(rmax)))*(10^floor(log10(rmax)));
    plotfmt_use.range=[-1 1]*rmax_round;
end
%    
np=size(dplot,2);
nrays=size(dplot,1);
if (np==1)
    fig_title=cat(2,opts_vismult.file_string,' ',opts_vismult.vis_string,': bid mults');
else
    fig_title=cat(2,opts_vismult.file_string,' ',opts_vismult.vis_string,': ray mults');
end
%
opts_vismult.fig_handle=figure;
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
        subplot(nrays,2,ip+(iray-1)*2);
        ipoints=intersect(find(rays.whichray==iray),sign_select{ip});
        if length(ipoints)>0
            %ilab=sa.typenames{ipoints(end)};
            ilab=psg_spec2legend(sa,ipoints(end),[]); %04Jul23
            %set up legends
            hl=cell(0);
            ht=[];
                    [rgb,symb]=psg_typenames2colors(sa.typenames(ipoints(end)),opts_tn2c);
                    hp=plot(dims_avail,squeeze(dplot(iray,ip,:)),cat(2,'k',symb));
                    set(hp,'Color',rgb);
                    set(hp,'MarkerSize',opts_vismult.marker_size);
                    hold on;
                    hl=[hl;hp];
                    ht=strvcat(ht,psg_spec2legend(sa,ipoints(1),[]));
                    hpline=plot(dims_avail,squeeze(dplot(iray,ip,:)),'k');
                    set(hpline,'Color',rgb);
            for k=1:length(plotfmt_use.lines)
                plot([-0.5 0.5]+opts_vismult.dim_range,repmat(plotfmt_use.lines(k),1,2),'k:');
            end
            if (opts_vismult.if_legend) & ~isempty(ht)
                hleg=legend(hl,ht);
                set(hleg,'FontSize',7);
%                set(hleg,'String',ht);
            end
            xlabel('dim');
            ylabel(opts_vismult.var_plot,'Interpreter','none');
            set(gca,'XLim',[-0.5 0.5]+opts_vismult.dim_range);
            set(gca,'XTick',dims_avail);
            set(gca,'YTick',sort([plotfmt_use.range,plotfmt_use.lines]));
            set(gca,'YLim',plotfmt_use.range);
            title(cat(2,ilab,' ',sign_char{ip}));
        end %ipoints>0
    end %ip
end %iray
axes('Position',[0.01,0.02,0.01,0.01]); %for text
text(0,0,fig_title,'Interpreter','none');
axis off;
%
opts_vismult_used=opts_vismult;
return


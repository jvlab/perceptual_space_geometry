% psg_raystats_dbplot_ticks
% utility fo rpsg_raystats_dbplot to set ticks and scales
%
% uses tick_posits,tick_labels,xlims,vars_avail,vars_sel,ivar_sel,if_angle
%
%  See also: PSG_RAYSTATS_DBPLOT.
%
set(gca,'XTick',tick_posits(:));
set(gca,'XTickLabel',tick_labels);
set(gca,'XLim',xlims);
%set vertical scales based on variable being plotted
if ~isempty(strfind(vars_avail{vars_sel(ivar_sel)},'dist_gain'))
    set(gca,'YLim',[0 max(get(gca,'YLim'))]);
end
if ~isempty(strfind(vars_avail{vars_sel(ivar_sel)},'cosang'))
    if (if_angle)
        set(gca,'YLim',[0 180]);
        set(gca,'YTick',[0:45:180]);
    else
        set(gca,'YLim',[-1 1]);
        set(gca,'YTick',[-1:0.5:1]);
    end
end

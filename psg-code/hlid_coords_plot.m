%hlid_coords_plot: script to plot coordinates derived by svd from signals
%in regions of interest
%
%30Oct24: add a plot with coords forced to have a positive average (can be eliminated by defining if_show_force_pos=0)
%06May25: fix variance plot so that it is s^2, not s
%07May25: forceed-positive plot also works with NaN's
%27May25: chang nstims to sum(nstims) for compatibility with hlid_orn_merge2
%
%   See also: HLID_RASTIM2COORDS_DEMO, HLID_RASTIM2COORDS_POOL, HLID_ORN_MERGE, HLID_ORN_MERGE2.
%
var_total=sum(s_diag_all.^2);
figname_raw=sprintf('raw responses: %s',dsid_show);
if if_submean
    figname_raw=cat(2,figname_raw,' mean subtracted');
else
    figname_raw=cat(2,figname_raw,' mean not subtracted');
end
figure;
set(gcf,'NumberTitle','off');
set(gcf,'Name',figname_raw);
set(gcf,'Position',[100 100 1200 800]);
%
subplot(3,1,1)
imagesc(resps);
xlabel('roi')
ylabel('stimulus labels');
set(gca,'YTick',[1:sum(nstims)]);
set(gca,'YTickLabel',stim_labels);
if exist('roi_names')
    set(gca,'XTick',[1:length(roi_names)]);
    set(gca,'XTickLabel',roi_names);
end
colorbar;
%
subplot(3,1,2)
imagesc(f.coord_opts.aux.v');
xlabel('roi')
ylabel('component');
set(gca,'YTick',[1:maxdim]);
if exist('roi_names')
    set(gca,'XTick',[1:length(roi_names)]);
    set(gca,'XTickLabel',roi_names);
end
colorbar;
%
subplot(3,3,7)
plot([1:length(s_diag_all)],s_diag_all.^2/var_total,'k');
hold on;
plot([1:length(s_diag_all)],s_diag_all.^2/var_total,'k*');
set(gca,'XLim',[0 maxdim_allowed+if_submean+0.5]);
set(gca,'XTick',[1:maxdim_allowed+if_submean]);
set(gca,'YLim',[0 1]);
xlabel('component');
ylabel('variance fraction');
title('variance explained by each component');
%
subplot(3,3,9);
imagesc(coords_all,[-1 1]*max(abs(coords_all(:))));
set(gca,'XTick',[1:maxdim]);
set(gca,'YTick',[1:sum(nstims)]);
set(gca,'YTickLabel',stim_labels);
title('coordinates')
colorbar;
%
if ~exist('if_show_force_pos')
    if_show_force_pos=1;
end
if (if_show_force_pos)
    % plot coordinates forcing avg to be positive
    subplot(3,3,8);
    imagesc(coords_all.*repmat(sign(sum(coords_all,1,'omitnan')),[size(coords_all,1) 1]),[-1 1]*max(abs(coords_all(:)))); %omitnan added 07May25
    set(gca,'XTick',[1:maxdim]);
    set(gca,'YTick',[1:sum(nstims)]);
    set(gca,'YTickLabel',stim_labels);
    title('coordinates, net pos pcs')
    colorbar;
end
    axes('Position',[0.01,0.02,0.01,0.01]); %for text
    text(0,0,figname_raw,'Interpreter','none');
    axis off;
%

function figh=psg_align_stats_plot(ra,ra_setup)
% figh=psg_align_stats_plot(ra,ra_setup) plots the results of psg_align_stats
%
% ra: results structure from psg_align_stats
% ra_setup: analysis setup, see psg_align_stats_demo
%   ra_setup.dim_list_in_max: maximum input dimension
%   ra_setup.dim_list_in: list of input dimensions, defaults to [1:ra_setup.dim_list_in_max]
%   ra_setup.shuff_quantiles: quantiles to plot
%   ra_setup.figh: figure handle to use; if empty, figure will be opened
%   ra_setup.nrows: number of rows
%   ra_setup.row: row to plot
%   ra_setup.range_equalize: 1 (default) to equalize the ranges across rows, active only if setup.row=setup.nrows
%
% figh: figure handle
%
%   See also: PSG_ALIGN_STATS, PSG_ALIGN_STATS_DEMO, PSG_KNIT_STATS_PLOT.
%
ra_setup=filldefault(ra_setup,'dim_list_in',[1:ra_setup.dim_list_in]);
ra_setup=filldefault(ra_setup,'figh',[]);
ra_setup=filldefault(ra_setup,'nrows',1);
ra_setup=filldefault(ra_setup,'row',1);
ra_setup=filldefault(ra_setup,'range_equalize',1);
ra_setup=filldefault(ra_setup,'shuff_quantiles',[0.01 0.05 0.5 0.95 0.99]);
if isempty(ra_setup.figh)
    figh=figure;
    set(gcf,'NumberTitle','off');
    set(gcf,'Name','consensus analysis');
    set(gcf,'Position',[100 100 1300 800]);
else
    figh=figure(ra_setup.figh);
end
ncols=4; %unexplained variance by dataset, unexplained variance by stimulus, unexplained variance and shuffles, explained variance and shuffles
rms_plot_max1=max([max(abs(ra.rmsdev_setwise(:))),max(abs(ra.rmsdev_stmwise(:)))]); %max possible rms variance per category (dataset or stim)
rms_plot_max2=max(abs(ra.rmsdev_overall(:)));
if ra_setup.nshuffs>0
    rms_plot_max2=max([rms_plot_max2,max(abs(ra.rmsdev_overall_shuff(:)))]);
end
rms_plot_max3=max(ra.rmsavail_overall);
%
allow_scale=ra.opts_pcon.allow_scale;
if_normscale=ra.opts_pcon.if_normscale;
% for allow_scale=0:1
%     ia=allow_scale+1;
if (allow_scale==0)
    scale_string='no scaling';
else
    scale_string='scaling';
    if if_normscale
        scale_string=cat(2,scale_string,'+norm');
    end
end
%compare rms devs across datasets
subplot(ra_setup.nrows,ncols,(ra_setup.row-1)*ncols+1);
imagesc(ra.rmsdev_setwise,[0 rms_plot_max1]);
set(gca,'TickLabelInterpreter','none');
xlabel('dataset');
set(gca,'XTick',1:ra_setup.nsets);
set(gca,'XTickLabel',ra_setup.dataset_labels);
ylabel('dim');
set(gca,'YTick',1:ra_setup.dim_list_in_max);
title(cat(2,'rms var unex, by set, ',scale_string));
colorbar;
%compare rms devs across stimuli
subplot(ra_setup.nrows,ncols,(ra_setup.row-1)*ncols+2);
imagesc(ra.rmsdev_stmwise,[0 rms_plot_max1]);
set(gca,'TickLabelInterpreter','none');
xlabel('stim');
set(gca,'XTick',1:ra_setup.nstims);
set(gca,'XTickLabel',ra_setup.stimulus_labels);
ylabel('dim');
set(gca,'YTick',1:ra_setup.dim_list_in_max);
title(cat(2,'rms var unex, by stim, ',scale_string));
colorbar;
for iue=1:2 %unexplained or explained
%overall rms as function of dimension, and compare with shuffles
    if (iue==1)
        var_string='unexplained'; 
        iue_sign=1; %logic to add unexplained variance
        iue_mult=0; %and not include available variance
        ylim=rms_plot_max2;
    else
        var_string='explained';
        iue_sign=-1; %logic to subtract unexplained variance
        iue_mult=1; %and add explained variance
        ylim=rms_plot_max3;
    end
    subplot(ra_setup.nrows,ncols,(ra_setup.row-1)*ncols+2+iue);
    hl=cell(0);
    hp=plot(ra_setup.dim_list_in,iue_mult*ra.rmsavail_overall(ra_setup.dim_list_in,1)+iue_sign*ra.rmsdev_overall(ra_setup.dim_list_in,1),'k');
    hl=[hl,hp];
    ht='consensus, data';
    hold on;
    if (iue==2)
        hp=plot(ra_setup.dim_list_in,ra.rmsavail_overall(ra_setup.dim_list_in,1),'b');
        hl=[hl,hp];
        ht=strvcat(ht,'avail');
    end
    if ra_setup.nshuffs>0
        nquantiles=length(ra_setup.shuff_quantiles);
        for iq=1:nquantiles
            switch sign(ra_setup.shuff_quantiles(iq)-0.5)
                case -1
                    linetype=':';
                case 0
                    linetype='';
                case 1
                    linetype='--';
            end
            hp_last=plot(ra_setup.dim_list_in,iue_mult*ra.rmsavail_overall(ra_setup.dim_list_in,1)+...
                iue_sign*quantile(ra.rmsdev_overall_shuff(ra_setup.dim_list_in,1,1,:,1),ra_setup.shuff_quantiles(iq),4),cat(2,'r',linetype));
            hp_all=plot(ra_setup.dim_list_in,iue_mult*ra.rmsavail_overall(ra_setup.dim_list_in,1)+...
                iue_sign*quantile(ra.rmsdev_overall_shuff(ra_setup.dim_list_in,1,1,:,2),ra_setup.shuff_quantiles(iq),4),cat(2,'m',linetype));
            if iq==round((1+nquantiles)/2)
                hl=[hl,hp_last,hp_all];
                ht=strvcat(ht,'cons, last','cons, all');
             end
        end
    end
    set(gca,'XTick',ra_setup.dim_list_in);
    set(gca,'XTickLabel',ra_setup.dim_list_in);
    set(gca,'XLim',[0 ra_setup.dim_list_in_max]);
    xlabel('dim');
    set(gca,'YLim',[0 ylim]);
    ylabel('rms dev');
    title(cat(2,'rms ',var_string,' overall, ',scale_string));
    legend(hl,ht,'Location','Best','FontSize',7);
end %iue
% end
if (ra_setup.row==ra_setup.nrows)
    axes('Position',[0.01,0.04,0.01,0.01]); %for text
    text(0,0,'consensus analysis','Interpreter','none','FontSize',8);
    axis off;
    if (ra_setup.nshuffs>0)
        axes('Position',[0.5,0.04,0.01,0.01]); %for text
        text(0,0,cat(2,sprintf('quantiles from %5.0f shuffles: ',ra_setup.nshuffs),sprintf('%6.4f ',ra_setup.shuff_quantiles)),...
            'FontSize',8);
        axis off;
    end
    if ra_setup.range_equalize
        %set equal scales on colorbars
        for col=1:2
            cl=0;
            for row=1:ra_setup.nrows
                cl=max(cl,max(get(subplot(ra_setup.nrows,ncols,col+(row-1)*ncols),'Clim')));
            end
            for row=1:ra_setup.nrows
                set(subplot(ra_setup.nrows,ncols,col+(row-1)*ncols),'Clim',[0 cl]);
            end
        end
        %set equal scales on line plots
        for col=3:4
            cl=0;
            for row=1:ra_setup.nrows
                cl=max(cl,max(get(subplot(ra_setup.nrows,ncols,col+(row-1)*ncols),'YLim')));
            end
            for row=1:ra_setup.nrows
                set(subplot(ra_setup.nrows,ncols,col+(row-1)*ncols),'YLim',[0 cl]);
            end
        end
    end
end

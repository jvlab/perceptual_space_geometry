%psg_align_vara_plot: plotting script for psg_align_vara_demo
%
%plotting: only uses quantities in results
% shuff_quantiles, if_debug, plot_max_factor can be redefined before running
%
% 11Dec24: add tags
%
%  See also: PSG_ALIGN_VARA_PLOT.
%
if ~exist('if_debug') if_debug=0; end
if ~exist('shuff_quantiles') shuff_quantiles=[0.01 0.05 0.5 0.95 0.99]; end %quantiles for showing shuffled data
if ~exist('plot_max_factor') plot_max_factor=1; end
%provide empty default for tags
if ~isfield(results,'opts_multi')
    results.opts_multi=struct;
end
if ~isfield(results.opts_multi,'tags')
    results.opts_multi.tags=[];
    results.opts_multi.if_tagged=0;
end
if results.opts_multi.if_tagged
    tags=results.opts_multi.tags;
else
    tags=[];
end
%
figure;
set(gcf,'NumberTitle','off');
set(gcf,'Name','variance analysis');
set(gcf,'Position',[100 100 1400 800]);
ncols=3+if_debug; %several kinds of plots
rms_plot_max1=max(abs(results.rmsdev_overall(:)));
rms_plot_max2=max(abs(results.rmsdev_grpwise(:)));
if results.nshuffs>0
    rms_plot_max2=max([rms_plot_max2,max(abs(results.rmsdev_grpwise_shuff(:)))]);
end
rms_plot_max=plot_max_factor*max(rms_plot_max1,rms_plot_max2);
if results.if_shuff_all
    shuff_all_string=' all';
else
    shuff_all_string='';
end
for allow_scale=0:1
    ia=allow_scale+1;
    if (allow_scale==0)
        scale_string='no scaling';
    else
        scale_string='scaling';
        if results.if_normscale
            scale_string=cat(2,scale_string,'+norm');
        end
    end
    %concatenate groups, and set up pointers for reordering
    rmsdev_setwise_gp_concat=zeros(results.dim_max,0);
    gp_reorder=[];
    for igp=1:results.ngps
        rmsdev_setwise_gp_concat=cat(2,rmsdev_setwise_gp_concat,results.rmsdev_setwise_gp(:,[1:results.nsets_gp(igp)],ia,igp));
        gp_reorder=[gp_reorder,results.gp_list{igp}];
    end
    %
    if (if_debug)
        %rms dev from global, by set, native order
        subplot(2,ncols,allow_scale*ncols+ncols);
        imagesc(results.rmsdev_setwise(:,:,ia),[0 rms_plot_max]);
        hold on;
        psg_align_vara_util(setfield(results,'nsets_gp',results.nsets),[1:results.nsets],tags);
        title(cat(2,'rms dev fr gbl, ',scale_string));
    end
    %
    %rms dev from global, by set, group order
    subplot(2,ncols,allow_scale*ncols+1);
    imagesc(results.rmsdev_setwise(:,gp_reorder,ia),[0 rms_plot_max]);
    hold on;
    psg_align_vara_util(results,gp_reorder,tags);
    title(cat(2,'rms dev fr gbl, ',scale_string));
    %
    %rms dev from its group, by set, group order
    subplot(2,ncols,allow_scale*ncols+2);
    imagesc(rmsdev_setwise_gp_concat,[0 rms_plot_max]);
    hold on;
    psg_align_vara_util(results,gp_reorder,tags);
    title(cat(2,'rms dev fr grp, ',scale_string));
    %
    %compare global and group-wise rms devs
    hl=cell(0);
    if (results.nshuffs>0)
        ht=strvcat('overall','within-group','shuff mean');
    else
        ht=strvcat('overall','within-group');
    end
    subplot(2,ncols,allow_scale*ncols+3);
    hp=plot([1:results.dim_max],results.rmsdev_overall(:,1,ia),'b');
    hl=[hl;hp];
    hold on;
    hp=plot([1:results.dim_max],results.rmsdev_grpwise(:,1,ia),'k');
    hl=[hl;hp];
    if results.nshuffs>0
        hp=plot([1:results.dim_max],mean(results.rmsdev_grpwise_shuff(:,1,ia,:,:),5),'r*');   
        hl=[hl;hp];
        quant_plot=quantile(reshape(results.rmsdev_grpwise_shuff(:,1,ia,:,:),[results.dim_max,results.nshuffs]),shuff_quantiles,2);
        for iq=1:length(shuff_quantiles)
            switch sign(shuff_quantiles(iq)-0.5)
                case -1
                    linetype=':';
                case 0
                    linetype='';
                case 1
                    linetype='--';
            end
            hp=plot([1:results.dim_max],quant_plot(:,iq),cat(2,'r',linetype));
            if iq==round((1+length(shuff_quantiles))/2)
                hl=[hl;hp];
                ht=strvcat(ht,'shuff quantile');
            end
        end
    end %shuff
    xlabel('dim');
    ylabel('rms dev')
    set(gca,'XTick',[1 results.dim_max])
    set(gca,'XLim',[0 results.dim_max]);
    set(gca,'XTick',[1:results.dim_max]);
    set(gca,'YLim',[0 rms_plot_max]);
    title(cat(2,'rms dev, ',scale_string));
    legend(hl,ht,'Location','Best','FontSize',7);
end
axes('Position',[0.01,0.04,0.01,0.01]); %for text
text(0,0,'variance analysis','Interpreter','none','FontSize',8);
axis off;
if (results.nshuffs>0)
    axes('Position',[0.01,0.01,0.01,0.01]); %for text
    text(0,0,cat(2,sprintf('quantiles from%s %5.0f shuffles: ',shuff_all_string,results.nshuffs),sprintf('%6.4f ',shuff_quantiles)),...
        'FontSize',8);
    axis off;
end

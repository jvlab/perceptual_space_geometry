function figh=psg_vara_stats_plot(vs,vs_setup)
% figh=psg_vara_stats_plot(vs,vs_setup) plots the results of rs_vara_coordsets
%
% vs: results structure from ra_vara_coordsets
% vs_setup: analysis setup, see ra_vara_coordsets
%   vs_setup.nshuffs: number of shuffles
%   vs_setup.dim_list_in_max: maximum input dimension
%     optional:
%   vs_setup.dataset_labels: labels for each dataset
%   vs_setup.stimulus_labels: labels for each stimulus
%   vs_setup.dim_list_in: list of input dimensions to plot, defaults to [1:vs_setup.dim_list_in_max]
%   vs_setup.dim_list_out: list of output dimensions, defaults to 1:vs_setup.dim_list_in
%   vs_setup.shuff_quantiles: quantiles to plot
%   vs_setup.figh: figure handle to use; if empty, figure will be opened
%   vs_setup.row: row to plot, defaults to 1
%   vs_setup.nrows: number of rows, defaults to max(1,vs_setup.row)
%   vs_setup.range_equalize: 1 (default) to equalize the ranges across rows, active only if setup.row=setup.nrows
%
%   See also: PSG_ALIGN_STATS, PSG_ALIGN_STATS_DEMO, PSG_KNIT_STATS_PLOT.
%
%
%  See also: PSG_ALIGN_VARA_PLOT.
%

vs_setup=filldefault(vs_setup,'dataset_labels',[]);
vs_setup=filldefault(vs_setup,'stimulus_labels',[]);
vs_setup=filldefault(vs_setup,'dim_list_in',[1:vs_setup.dim_list_in_max]);
vs_setup=filldefault(vs_setup,'dim_list_out',vs_setup.dim_list_in);
vs_setup=filldefault(vs_setup,'figh',[]);
vs_setup=filldefault(vs_setup,'row',1);
vs_setup=filldefault(vs_setup,'nrows',max(1,vs_setup.row));
vs_setup=filldefault(vs_setup,'range_equalize',1);
vs_setup=filldefault(vs_setup,'shuff_quantiles',[0.01 0.05 0.5 0.95 0.99]);
vs_setup=filldefault(vs_setup,'if_debug',0);
vs_setup=filldefault(vs_setup,'plot_max_factor',1);
%
if isempty(vs_setup.figh)
    figh=figure;
    set(gcf,'NumberTitle','off');
    set(gcf,'Name','variance analysis');
    set(gcf,'Position',[100 100 1400 800]);
else
    figh=figure(vs_setup.figh);
end
%
tags=vs.groupings.tags;
ncols=3+vs_setup.if_debug;
%
rms_plot_max1=max(abs(vs.rmsdev_overall(:)));
rms_plot_max2=max(abs(vs.rmsdev_grpwise(:)));
if vs.nshuffs>0
    rms_plot_max2=max([rms_plot_max2,max(abs(vs.rmsdev_grpwise_shuff(:)))]);
end
rms_plot_max=vs_setup.plot_max_factor*max(rms_plot_max1,rms_plot_max2);
if isfield(vs.opts_multi_used,'if_exhaust')
    if vs.opts_multi_used.if_exhaust==1
        shuff_desc_string='exhaustive';
    else
        shuff_desc_string='random';
    end
else
    shuff_desc_string='custom';
end
if vs_setup.opts_vara.allow_scale==0
    scale_string='no scaling';
else
    scale_string='scaling';
end
if vs_setup.opts_vara.if_normscale
    scale_string=cat(2,scale_string,'+norm');
end
dim_max=max(vs_setup.dim_list_out);
%concatenate groups, and set up pointers for reordering
rmsdev_setwise_gp_concat=zeros(dim_max,0);
gp_reorder=[];
for igp=1:vs.groupings.ngps
    rmsdev_setwise_gp_concat=cat(2,rmsdev_setwise_gp_concat,vs.rmsdev_setwise_gp(:,[1:vs.groupings.nsets_gp(igp)],igp));
    gp_reorder=[gp_reorder,vs.groupings.gp_list{igp}];
end
%
vsu=vs; %versionfor psg_align_vara_util
vsu=setfield(vsu,'nsets',vs_setup.nsets);
vsu=setfield(vsu,'gps',vs.groupings.gps);
vsu=setfield(vsu,'dataset_labels',vs_setup.dataset_labels);
vsu=setfield(vsu,'dim_max',dim_max);
if vs_setup.if_debug %rms dev from global, by set, native order
    subplot(vs_setup.nrows,ncols,ncols*(vs_setup.row-1)+ncols);
    imagesc(vs.rmsdev_setwise,[0 rms_plot_max]);
    hold on;
 %       psg_align_vara_util(setfield(results,'nsets_gp',results.nsets),[1:results.nsets],tags);
    psg_align_vara_util(setfield(vsu,'nsets_gp',vs_setup.nsets),[1:vs_setup.nsets],tags);
    title(cat(2,'rms dev fr gbl, ',scale_string));
end
return
end
%     %
%     %rms dev from global, by set, group order
%     subplot(2,ncols,allow_scale*ncols+1);
%     imagesc(results.rmsdev_setwise(:,gp_reorder,ia),[0 rms_plot_max]);
%     hold on;
%     psg_align_vara_util(results,gp_reorder,tags);
%     title(cat(2,'rms dev fr gbl, ',scale_string));
%     %
%     %rms dev from its group, by set, group order
%     subplot(2,ncols,allow_scale*ncols+2);
%     imagesc(rmsdev_setwise_gp_concat,[0 rms_plot_max]);
%     hold on;
%     psg_align_vara_util(results,gp_reorder,tags);
%     title(cat(2,'rms dev fr grp, ',scale_string));
%     %
%     %compare global and group-wise rms devs
%     hl=cell(0);
%     if (results.nshuffs>0)
%         ht=strvcat('overall','within-group','shuff mean');
%     else
%         ht=strvcat('overall','within-group');
%     end
%     subplot(2,ncols,allow_scale*ncols+3);
%     hp=plot([1:dim_max],results.rmsdev_overall(:,1,ia),'b');
%     hl=[hl;hp];
%     hold on;
%     hp=plot([1:dim_max],results.rmsdev_grpwise(:,1,ia),'k');
%     hl=[hl;hp];
%     if results.nshuffs>0
%         hp=plot([1:dim_max],mean(results.rmsdev_grpwise_shuff(:,1,ia,:,:),5),'r*');   
%         hl=[hl;hp];
%         quant_plot=quantile(reshape(results.rmsdev_grpwise_shuff(:,1,ia,:,:),[dim_max,results.nshuffs]),shuff_quantiles,2);
%         for iq=1:length(shuff_quantiles)
%             switch sign(shuff_quantiles(iq)-0.5)
%                 case -1
%                     linetype=':';
%                 case 0
%                     linetype='';
%                 case 1
%                     linetype='--';
%             end
%             hp=plot([1:dim_max],quant_plot(:,iq),cat(2,'r',linetype));
%             if iq==round((1+length(shuff_quantiles))/2)
%                 hl=[hl;hp];
%                 ht=strvcat(ht,'shuff quantile');
%             end
%         end
%     end %shuff
%     xlabel('dim');
%     ylabel('rms dev')
%     set(gca,'XTick',[1 dim_max])
%     set(gca,'XLim',[0 dim_max]);
%     set(gca,'XTick',[1:dim_max]);
%     set(gca,'YLim',[0 rms_plot_max]);
%     title(cat(2,'rms dev, ',scale_string));
%     legend(hl,ht,'Location','Best','FontSize',7);
% end
% axes('Position',[0.01,0.04,0.01,0.01]); %for text
% text(0,0,'variance analysis','Interpreter','none','FontSize',8);
% axis off;
% if (results.nshuffs>0)
%     axes('Position',[0.01,0.01,0.01,0.01]); %for text
%     text(0,0,cat(2,sprintf('quantiles from%s %5.0f shuffles: ',shuff_all_string,results.nshuffs),sprintf('%6.4f ',shuff_quantiles)),...
%         'FontSize',8);
%     axis off;
% end
%%%%%%%%%%%%%%%%%
%from psg_align_stats_plot
% 
% 
% vs_setup=filldefault(vs_setup,'dataset_labels',[]);
% vs_setup=filldefault(vs_setup,'stimulus_labels',[]);
% vs_setup=filldefault(vs_setup,'dim_list_in',[1:vs_setup.dim_list_in_max]);
% vs_setup=filldefault(vs_setup,'dim_list_out',vs_setup.dim_list_in);
% vs_setup=filldefault(vs_setup,'figh',[]);
% vs_setup=filldefault(vs_setup,'row',1);
% vs_setup=filldefault(vs_setup,'nrows',max(1,vs_setup.row));
% vs_setup=filldefault(vs_setup,'range_equalize',1);
% vs_setup=filldefault(vs_setup,'shuff_quantiles',[0.01 0.05 0.5 0.95 0.99]);
% if isempty(vs_setup.figh)
%     figh=figure;
%     set(gcf,'NumberTitle','off');
%     set(gcf,'Name','consensus analysis');
%     set(gcf,'Position',[100 100 1300 800]);
% else
%     figh=figure(vs_setup.figh);
% end
% ncols=4; %unexplained variance by dataset, unexplained variance by stimulus, unexplained variance and shuffles, explained variance and shuffles
% rms_plot_max1=max([max(abs(vs.rmsdev_setwise(:))),max(abs(vs.rmsdev_stmwise(:)))]); %max possible rms variance per category (dataset or stim)
% rms_plot_max2=max(abs(vs.rmsdev_overall(:)));
% if vs_setup.nshuffs>0
%     rms_plot_max2=max([rms_plot_max2,max(abs(vs.rmsdev_overall_shuff(:)))]);
% end
% rms_plot_max3=max(vs.rmsavail_overall);
% %
% allow_scale=vs.opts_pcon.allow_scale;
% if_normscale=vs.opts_pcon.if_normscale;
% % for allow_scale=0:1
% %     ia=allow_scale+1;
% if (allow_scale==0)
%     scale_string='no scaling';
% else
%     scale_string='scaling';
%     if if_normscale
%         scale_string=cat(2,scale_string,'+norm');
%     end
% end
% %compare rms devs across datasets
% subplot(vs_setup.nrows,ncols,(vs_setup.row-1)*ncols+1);
% imagesc(vs.rmsdev_setwise,[0 rms_plot_max1]);
% set(gca,'TickLabelInterpreter','none');
% xlabel('dataset');
% if ~isempty(vs_setup.dataset_labels)
%     set(gca,'XTick',1:vs_setup.nsets);
%     set(gca,'XTickLabel',vs_setup.dataset_labels);
% end
% ylabel('dim');
% set(gca,'YTick',1:vs_setup.dim_list_in_max);
% title(cat(2,'rms var unex, by set, ',scale_string));
% colorbar;
% %compare rms devs across stimuli
% subplot(vs_setup.nrows,ncols,(vs_setup.row-1)*ncols+2);
% imagesc(vs.rmsdev_stmwise,[0 rms_plot_max1]);
% set(gca,'TickLabelInterpreter','none');
% xlabel('stim');
% if ~isempty(vs_setup.stimulus_labels)
%     set(gca,'XTick',1:vs_setup.nstims);
%     set(gca,'XTickLabel',vs_setup.stimulus_labels);
% end
% ylabel('dim');
% set(gca,'YTick',1:vs_setup.dim_list_in_max);
% title(cat(2,'rms var unex, by stim, ',scale_string));
% colorbar;
% for iue=1:2 %unexplained or explained
% %overall rms as function of dimension, and compare with shuffles
%     if (iue==1)
%         var_string='unexplained'; 
%         iue_sign=1; %logic to add unexplained variance
%         iue_mult=0; %and not include available variance
%         ylim=rms_plot_max2;
%     else
%         var_string='explained';
%         iue_sign=-1; %logic to subtract unexplained variance
%         iue_mult=1; %and add explained variance
%         ylim=rms_plot_max3;
%     end
%     subplot(vs_setup.nrows,ncols,(vs_setup.row-1)*ncols+2+iue);
%     hl=cell(0);
%     hp=plot(vs_setup.dim_list_in,iue_mult*vs.rmsavail_overall(vs_setup.dim_list_in,1)+iue_sign*vs.rmsdev_overall(vs_setup.dim_list_in,1),'k');
%     hl=[hl,hp];
%     ht='consensus, data';
%     hold on;
%     if (iue==2)
%         hp=plot(vs_setup.dim_list_in,vs.rmsavail_overall(vs_setup.dim_list_in,1),'b');
%         hl=[hl,hp];
%         ht=strvcat(ht,'avail');
%     end
%     if vs_setup.nshuffs>0
%         nquantiles=length(vs_setup.shuff_quantiles);
%         for iq=1:nquantiles
%             switch sign(vs_setup.shuff_quantiles(iq)-0.5)
%                 case -1
%                     linetype=':';
%                 case 0
%                     linetype='';
%                 case 1
%                     linetype='--';
%             end
%             hp_last=plot(vs_setup.dim_list_in,iue_mult*vs.rmsavail_overall(vs_setup.dim_list_in,1)+...
%                 iue_sign*quantile(vs.rmsdev_overall_shuff(vs_setup.dim_list_in,1,1,:,1),vs_setup.shuff_quantiles(iq),4),cat(2,'r',linetype));
%             hp_all=plot(vs_setup.dim_list_in,iue_mult*vs.rmsavail_overall(vs_setup.dim_list_in,1)+...
%                 iue_sign*quantile(vs.rmsdev_overall_shuff(vs_setup.dim_list_in,1,1,:,2),vs_setup.shuff_quantiles(iq),4),cat(2,'m',linetype));
%             if iq==round((1+nquantiles)/2)
%                 hl=[hl,hp_last,hp_all];
%                 ht=strvcat(ht,'cons, last','cons, all');
%              end
%         end
%     end
%     set(gca,'XTick',vs_setup.dim_list_in);
%     set(gca,'XTickLabel',vs_setup.dim_list_in);
%     set(gca,'XLim',[0 vs_setup.dim_list_in_max]);
%     xlabel('dim in');
%     set(gca,'YLim',[0 ylim]);
%     ylabel('rms dev');
%     title(cat(2,'rms ',var_string,' overall, ',scale_string));
%     legend(hl,ht,'Location','Best','FontSize',7);
% end %iue
% % end
% if (vs_setup.row==vs_setup.nrows)
%     axes('Position',[0.01,0.04,0.01,0.01]); %for text
%     text(0,0,'consensus analysis','Interpreter','none','FontSize',8);
%     axis off;
%     %show input and output dimensions if they differ
%     if ~all(vs_setup.dim_list_in==vs_setup.dim_list_out)
%         dim_in_string=cat(2,'dim  in:',sprintf(' %1.0f',vs_setup.dim_list_in));
%         dim_out_string=cat(2,'dim out:',sprintf(' %1.0f',vs_setup.dim_list_out));
%     end
% end
% if (vs_setup.nshuffs>0)
%     voff=(vs_setup.nrows-vs_setup.row)/vs_setup.nrows;
%     axes('Position',[0.5,0.04+voff,0.01,0.01]); %for text
%     text(0,0,cat(2,sprintf('quantiles from %5.0f shuffles: ',vs_setup.nshuffs),sprintf('%6.4f ',vs_setup.shuff_quantiles)),...
%         'FontSize',8);
%     axis off;
%     %show input and output dimensions if they differ
%     if ~all(vs_setup.dim_list_in==vs_setup.dim_list_out)
%         dim_in_string=cat(2,'dim  in:',sprintf(' %1.0f',vs_setup.dim_list_in));
%         dim_out_string=cat(2,'dim out:',sprintf(' %1.0f',vs_setup.dim_list_out));
%         axes('Position',[0.75,0.04+voff,0.01,0.01]); %for text
%         text(0,0,dim_in_string,'FontSize',8);
%         axis off;
%         axes('Position',[0.75,0.02+voff,0.01,0.01]); %for text
%         text(0,0,dim_out_string,'FontSize',8);
%         axis off;
%     end
% end
% if vs_setup.range_equalize & vs_setup.row==vs_setup.nrows
%     %set equal scales on colorbars
%     for col=1:2
%         cl=0;
%         for row=1:vs_setup.nrows
%             cl=max(cl,max(get(subplot(vs_setup.nrows,ncols,col+(row-1)*ncols),'Clim')));
%         end
%         for row=1:vs_setup.nrows
%             set(subplot(vs_setup.nrows,ncols,col+(row-1)*ncols),'Clim',[0 cl]);
%         end
%     end
%     %set equal scales on line plots
%     for col=3:4
%         cl=0;
%         for row=1:vs_setup.nrows
%             cl=max(cl,max(get(subplot(vs_setup.nrows,ncols,col+(row-1)*ncols),'YLim')));
%         end
%         for row=1:vs_setup.nrows
%             set(subplot(vs_setup.nrows,ncols,col+(row-1)*ncols),'YLim',[0 cl]);
%         end
%     end
% end

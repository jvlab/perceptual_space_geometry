function [opts_used,tables_selected,fighs,res]=btcsel_like_analtable(table_like,opts)
% [opts_used,tables_selected,fighs,res]=btcsel_like_analtable(table_like,opts): analyze results of log likelihoods of btcsel (btc, with selected stimuli)
%
% table_like: likelihood table, typically created by psg_like_maketable.
% opts: options -- all can be omitted and if omitted, will be requested
%  thr_type_choice; threshold type (min, max, avg)
%  frac_keep_choices: 1-> first value of frac keep, 2-> second value etc,
%  others retained from psg_like_analtable
%
% opts_used: options used
% tables_selected: cell array, tables_selected{illr,k} is table_like, after selecting subjects, threshold type, threshold values,
%   and log likelihood type (1: sym, 2: umi, 3: adt)
% fighs: figure handles
% res: analysis results (not used at present)
%
% some categorical variables are converted to tokens, as defined in the
% 'tokens' structure and kept as UserData in table_like.
%
% llr quantities for umi are corrected, i.e., have log(h) subtracted
%
%   See also:  PSG_UMI_TRIPLIKE_DEMO, PPSG_TENTLIKE_DEMO, PSG_COLORS_LIKE, PSG_LIKE_MAKETABLE, PSG_LIKE_ANALTABLE, PSG_SELECT_CHOICEDATA.
%
res=[];
%
if (nargin<1)
    table_like=table();
end
if (nargin<2)
    opts=struct;
end
opts=filldefault(opts,'fn_table_def','psg_like_maketable_btc_24Jun23.mat');
%
colors_like=psg_colors_like;
paradigm_colors=colors_like.paradigm_colors;
subj_symbs_res=colors_like.subj_symbs_res; %reserved symbols
subj_fills_res=colors_like.subj_fills_res; %reserved fills
subj_symbs_unres=colors_like.subj_symbs_unres; %unreserved symbols
subj_fills_unres=colors_like.subj_fills_unres; %unreserved fills
apriori_symb=colors_like.apriori_symb; %symbol for a priori
%
opts=filldefault(opts,'paradigm_colors',paradigm_colors);
opts=filldefault(opts,'subj_symbs_res',subj_symbs_res);
opts=filldefault(opts,'subj_symbs_unres',subj_symbs_unres);
opts=filldefault(opts,'subj_fills_res',subj_fills_res);
opts=filldefault(opts,'subj_fills_unres',subj_fills_unres);
opts=filldefault(opts,'apriori_symb',apriori_symb);
opts=filldefault(opts,'if_surrogates',1); %at most one of [if_surrogates, if_sub_flip_all, if_sub_flip_any, if_sub_flip_one} can be 1
opts=filldefault(opts,'if_sub_flip_all',0);
opts=filldefault(opts,'if_sub_flip_any',0);
opts=filldefault(opts,'if_sub_flip_one',0); 
opts=filldefault(opts,'if_plot_conform',0); %1 to plot conforming surrogate rather than data
opts=filldefault(opts,'if_apriori',1);
opts=filldefault(opts,'if_stdevs',1);
opts=filldefault(opts,'if_stdevs_data',opts.if_stdevs);
opts=filldefault(opts,'if_stdevs_surrogates',opts.if_stdevs);
paradigm_colors=opts.paradigm_colors;
subj_symbs_res=opts.subj_symbs_res;
subj_fills_res=opts.subj_fills_res;
subj_symbs_unres=opts.subj_symbs_unres;
subj_fills_unres=opts.subj_fills_unres;
apriori_symb=opts.apriori_symb;
%
if opts.if_surrogates+opts.if_sub_flip_all+opts.if_sub_flip_any + opts.if_sub_flip_one >1
    warning('at most one of [if_surrogates, if_sub_flip_all, if_sub_flip_any, if_sub_flip_one} can be 1');
    return
end
%plot formatting
opts=filldefault(opts,'box_halfwidth',0.02); %half-width of boxes for s.d. of surrogates
opts=filldefault(opts,'plotrange_a',[0 1.25]);
opts=filldefault(opts,'plotrange_llr',[-2 .1]);
opts=filldefault(opts,'plotrange_llr_sub',[-1 1]);
opts=filldefault(opts,'if_plot_ah_llr',1); %1 to also plot log likelihood of Dirichlet fit
opts=filldefault(opts,'if_plot3d_h',0); %1 to plot h as third axis 
opts=filldefault(opts,'view3d',[-62 13]);
opts=filldefault(opts,'plotrange_h',[0 .2]);
if isempty(table_like)
    fn_table=getinp('likelihood table file name','s',[],opts.fn_table_def);
    load(fn_table);
    disp(sprintf('loaded table with %6.0f rows and %3.0f columns.',size(table_like)));
    disp(table_like.Properties.UserData.tokens);
    disp(table_like.Properties.UserData.notes);
    opts_used.fn_table=fn_table;
else
    opts_used.fn_table=[];
end
%
%modify variables for forward compatibility with conform, which takes
%variable names directly from the output of psg_[umi_trip|tent]like_demo
%
table_like.Properties.VariableNames=strrep(table_like.Properties.VariableNames,'llr_data','orig_data');
table_like.Properties.VariableNames=strrep(table_like.Properties.VariableNames,'llr_flip_','flip_');
%
tokens=table_like.Properties.UserData.tokens;
%
table_selected=table_like;

%select paradigm type
paradigm_types_avail=unique(table_like{:,'paradigm_type'});
opts=filldefault(opts,'paradigm_type_choice',[]);
if isempty(opts.paradigm_type_choice)
    for ip=1:length(paradigm_types_avail)
        disp(sprintf('%2.0f->%s',ip,paradigm_types_avail{ip}));
    end
    opts.paradigm_type_choice=getinp('choice','d',[1 length(paradigm_types_avail)]);
end
paradigm_type_choice=opts.paradigm_type_choice;
paradigm_type=paradigm_types_avail{paradigm_type_choice};
table_selected=table_selected(strmatch(paradigm_type,table_selected.paradigm_type,'exact'),:);
%
%select subject(s)
subj_ids_avail=unique(table_like{:,'subj_id'});
opts=filldefault(opts,'subj_id_choice',[]);
if isempty(opts.subj_id_choice)
    for ip=1:length(subj_ids_avail)
        disp(sprintf('%2.0f->%s',ip,subj_ids_avail{ip}));
    end
    opts.subj_id_choice=getinp('choice numbers','d',[1 length(subj_ids_avail)],[1:length(subj_ids_avail)]);
end
subj_id_choice=opts.subj_id_choice;
subj_id_list=table2cell(table_like(:,'subj_id'));
subj_sel=[];
for ic=1:length(subj_id_choice)
    subj_sel=[subj_sel;strmatch(subj_ids_avail{subj_id_choice(ic)},subj_id_list,'exact')];
end
table_selected=table_selected(subj_sel,:);
%
%select threshold type
thr_types=tokens.thr_type;
opts=filldefault(opts,'thr_type_choice',[]);
if isempty(opts.thr_type_choice)
    for ip=1:length(thr_types)
        disp(sprintf('%2.0f->%s',ip,thr_types{ip}));
    end
    opts.thr_type_choice=getinp('choice','d',[1 length(thr_types)]);
end
thr_type_choice=opts.thr_type_choice;
%
table_selected=table_selected(table_selected.thr_type==thr_type_choice,:);
%
%select fraction to keep
frac_keeps=flipud(unique(table_selected.frac_keep));
opts=filldefault(opts,'frac_keep_choices',[]);
if isempty(opts.frac_keep_choices)
    for ip=1:length(frac_keeps)
        disp(sprintf('%2.0f->keep fraction %8.6f',ip,frac_keeps(ip)));
    end
    opts.frac_keep_choices=getinp('choice(s)','d',[1 length(frac_keeps)]);
end
frac_keep_choices=intersect(opts.frac_keep_choices,[1:length(frac_keeps)]);
frac_keep_list=frac_keeps(frac_keep_choices);
%
nc_plot=length(tokens.llr_type)+opts.if_plot_ah_llr;
nr_plot=length(tokens.ipchoice);
tables_selected=cell(length(tokens.llr_type),length(frac_keep_list));
%
for ifk_ptr=1:length(frac_keep_list)
    frac_keep=frac_keep_list(ifk_ptr);
    table_fk=table_selected(table_selected.frac_keep==frac_keep,:);
    for illr=1:length(tokens.llr_type)
        table_fk_llr=table_fk(table_fk.llr_type==illr,:);
        tables_selected{illr,ifk_ptr}=table_fk_llr;
        disp(sprintf('frac keep = %7.5f, llr type: %s',frac_keep,tokens.llr_type{illr}));
    end
end
for illr=1:length(tokens.llr_type)
    table_llr=table_selected(table_selected.llr_type==illr,:);
    %find unique paradigm names 
    paradigm_names_avail=unique(table_llr{:,'paradigm_name'});
    for ipn=1:length(paradigm_names_avail)
         pn_sel=strmatch(paradigm_names_avail{ipn},table_llr.paradigm_name,'exact');
         table_llr_pn=table_llr(pn_sel,:);
         ip_list=unique(table_llr_pn{:,'ipchoice'});
         for ipp=1:length(ip_list)
            table_llr_pn_ip=table_llr_pn(table_llr_pn.ipchoice==ip_list(ipp),:);
            disp(sprintf(' llr type: %5s    paradigm %20s, fit type %10s',tokens.llr_type{illr},paradigm_names_avail{ipn},tokens.ipchoice{ipp}));
            if_header=0;
            for ifk_ptr=1:length(frac_keep_list)
                frac_keep=frac_keep_list(ifk_ptr);
                table_llr_pn_ip_fk=table_llr_pn_ip(table_llr_pn_ip.frac_keep==frac_keep,:);
                    if size(table_llr_pn_ip_fk,1)>0 %at this point, only difference in metadata should be subject ID
                    table_analyze=table_llr_pn_ip_fk;
                    nsubjs=length(unique(table_analyze{:,'subj_id'}));
                    if size(table_analyze,1)~=nsubjs
                        warning('more than one entry per subject');
                        disp(table_analyze);
                    else
                        if (if_header==0)
                            disp('keep frac  nsubjs   ntriads(sem)  ntrials(sem)     a (sem)           h (sem)     a_priori_llr(sem)');
                            if_header=1;
                        end
                        disp(sprintf('%7.5f  %5.0f   %7.1f(%4.2f) %7.1f(%4.2f)  %6.4f (%6.4f)   %6.4f (%6.4f)   %6.4f (%6.4f)',frac_keep,nsubjs,...
                            mean(table_analyze{:,'ntriads'}),std(table_analyze{:,'ntriads'})/sqrt(nsubjs),...
                            mean(table_analyze{:,'ntrials'}),std(table_analyze{:,'ntrials'})/sqrt(nsubjs),...
                            mean(table_analyze{:,'a'}),std(table_analyze{:,'a'})/sqrt(nsubjs),...
                            mean(table_analyze{:,'h'}),std(table_analyze{:,'h'})/sqrt(nsubjs),...
                            mean(table_analyze{:,'apriori_llr'}),std(table_analyze{:,'apriori_llr'})/sqrt(nsubjs)...
                            ));
                    end
                end
            end %ifk_ptr           
         end %ipp
    end %ipn
end %illr
             
    
%     %find unique paradigm names 
%         paradigm_names_avail=unique(table_fk_llr{:,'paradigm_name'});
%         for ipn=1:length(paradigm_names_avail)
%             pn_sel=strmatch(paradigm_names_avail{ipn},table_fk_llr.paradigm_name,'exact');
%             table_fk_llr_pn=table_fk_llr(pn_sel,:);
%             subjs_avail=unique(table_fk_llr_pn{:,'subj_id'});
%             disp(sprintf('   for %20s, number of subjects: %2.0f',paradigm_names_avail{ipn},length(subjs_avail)));
%             %
%             ip_list=unique(table_fk_llr_pn{:,'ipchoice'});
%             for ipp=1:length(ip_list)
%                 table_analyze=table_fk_llr_pn(table_fk_llr_pn.ipchoice==ip_list(ipp),:);
%                 table_vals=table_analyze(:,{'a','h','apriori_llr','thr_val','ntriads','ntrials','orig_data','orig_data_sd','flip_all','flip_all_sd','flip_any','flip_any_sd'}); 
%                 disp(table_vals)
%             end
%         end %ipn
%     end %llr_type
    

    
%     %
%     %make plots
%     
%     %each plot has a panel for the three kinds of llr_types, and the two kinds of ipchoice
%     %
%     tstring=sprintf('%s: threshold based on %s, frac keep %8.6f',paradigm_type,thr_types{thr_type_choice},frac_keep);
%     fighs(end+1)=figure;
%     set(gcf,'Position',[50 50 1200 850]);
%     set(gcf,'NumberTitle','off');
%     set(gcf,'Name',tstring);
%     subj_symbs=subj_symbs_res; %list of subject symbols, starting with reserved list
%     subj_fills=subj_fills_res;
%     for ipchoice=1:length(tokens.ipchoice)
%         for llr_type=(1-opts.if_plot_ah_llr):length(tokens.llr_type)
%             table_plot=table_fk(intersect(find(table_fk.ipchoice==ipchoice),find(table_fk.llr_type==max(1,llr_type))),:);
%             subj_ids=unique(table_plot.subj_id);
%             paradigm_names=unique(table_plot.paradigm_name);
%             if llr_type==strmatch('umi',tokens.llr_type,'exact')
%                 llr_label_suffix=' - log(h)';
%             else
%                 llr_label_suffix=' ';
%             end
%             sub_name='';
%             if llr_type>0
%                 if opts.if_sub_flip_all
%                     llr_label_suffix=cat(2,llr_label_suffix,' - flip all');
%                     sub_name='flip_all';
%                 end
%                 if opts.if_sub_flip_any
%                     llr_label_suffix=cat(2,llr_label_suffix,' - flip any');
%                     sub_name='flip_any';
%                 end
%                 if opts.if_sub_flip_one
%                     llr_label_suffix=cat(2,llr_label_suffix,' - flip one (conform)');
%                     sub_name='flip_one';
%                 end
%             end
%             if size(table_plot,1)>0
%                 paradigms_shown=[];
%                 subjs_shown=[];
%                 isub=llr_type+opts.if_plot_ah_llr+(ipchoice-1)*nc_plot;
%                 subplot(nr_plot,nc_plot,isub);
%                 nsubj_unres=0; %number of subjects with unreserved symbols
%                 legh=[];
%                 legt=[];
%                 for ipt=1:size(table_plot,1)
%                     data=table_plot(ipt,:);
%                     paradigm_name=data.paradigm_name{1};
%                     paradigm_color=paradigm_colors.(paradigm_name);
%                     %
%                     subj_name=data.subj_id{1};
%                     if ~isempty(strmatch(subj_name,fieldnames(subj_symbs),'exact'))
%                         subj_symb=subj_symbs.(subj_name);
%                         subj_fill=subj_fills.(subj_name);
%                     else
%                         nsubj_unres=nsubj_unres+1;
%                         subj_symb=subj_symbs_unres(1+mod(nsubj_unres-1,length(subj_symbs_unres)));
%                         subj_symbs.(subj_name)=subj_symb;
%                         subj_fill=subj_fills_unres(1+mod(nsubj_unres-1,length(subj_symbs_unres)));
%                         subj_fills.(subj_name)=subj_fill;
%                     end
%                     if (llr_type>0)
%                         %
%                         switch opts.if_plot_conform
%                             case 0
%                                 llr_plot=data.orig_data;
%                                 sd_plot=data.orig_data_sd;
%                                 title_suffix='';
%                             case 1
%                                 llr_plot=data.flip_one;
%                                 sd_plot=data.flip_one_sd;
%                                 title_suffix='(flip one conform)';
%                         end
%                         llr_plot_sub=0;
%                         if ~isempty(sub_name)
%                             llr_plot_sub=data.(sub_name);
%                             sd_plot=sqrt(sd_plot.^2+data.(cat(2,sub_name,'_sd')).^2); %revise standard dev of a difference from Gaussian approx
%                         end
%                         hp=btcsel_like_plot(data.a,llr_plot-llr_plot_sub,cat(2,'k',subj_symb),data.h,d23);
%                         set(hp,'Color',paradigm_color);
%                         if (subj_fill)
%                             set(hp,'MarkerFaceColor',paradigm_color);
%                         end
%                         if any(sd_plot(:)>0) & opts.if_stdevs_data %only plot standard devs if present
%                             hp=btcsel_like_plot(repmat(data.a,1,2),llr_plot-llr_plot_sub+sd_plot*[-1 1],'k',data.h,d23);
%                             set(hp,'Color',paradigm_color);
%                         end
%                         %plot surrogates
%                         if opts.if_surrogates
%                             hs=btcsel_like_plot(repmat(data.a,1,3),[llr_plot data.flip_all,data.flip_any],'k',data.h,d23);
%                             set(hs,'Color',paradigm_color);
%                             if opts.if_stdevs_surrogates
%                                 hs=btcsel_like_plot(data.a+opts.box_halfwidth*[-1 1 1 -1 -1],data.flip_all+data.flip_all_sd*[1 1 -1 -1 1],'k',data.h,d23);
%                                 set(hs,'Color',paradigm_color);
%                                 hs=btcsel_like_plot(data.a+opts.box_halfwidth*[-2 2 2 -2 -2],data.flip_any+data.flip_any_sd*[1 1 -1 -1 1],'k',data.h,d23);
%                                 set(hs,'Color',paradigm_color);
%                             end
%                         end
%                         if opts.if_apriori
%                             %plot a priori
%                             ha=btcsel_like_plot(data.a,data.apriori_llr-llr_plot_sub,cat(2,'k',apriori_symb),data.h,d23);
%                             set(ha,'Color',paradigm_color);
%                         end
%                     else
%                         hp=btcsel_like_plot(data.a,data.ah_llr,cat(2,'k',subj_symb),data.h,d23);
%                         set(hp,'Color',paradigm_color);
%                         if (subj_fill)
%                             set(hp,'MarkerFaceColor',paradigm_color);
%                         end
%                     end
%                     if_addleg=0;
%                     if isempty(strmatch(paradigm_name,paradigms_shown,'exact'))
%                         paradigms_shown=strvcat(paradigms_shown,paradigm_name);
%                         if_addleg=1;
%                     end
%                     if isempty(strmatch(subj_name,subjs_shown,'exact'))
%                         subjs_shown=strvcat(subjs_shown,subj_name);
%                         if_addleg=1;
%                     end
%                     if (if_addleg)
%                         legt=strvcat(legt,cat(2,paradigm_name,' ',subj_name));
%                         legh=[legh;hp];
%                     end
%                     if (llr_type>0)
%                         title(cat(2,tokens.llr_type{llr_type},' ',tokens.ipchoice{ipchoice},' ',title_suffix));
%                     else
%                         title(cat(2,'llr for Dirichlet',' ',tokens.ipchoice{ipchoice}));
%                     end
%                     xlabel('Dirichlet a');
%                     set(gca,'XTick',[opts.plotrange_a(1):.25:opts.plotrange_a(2)]);
%                     set(gca,'XLim',opts.plotrange_a);
%                     llr_lim=opts.plotrange_llr;
%                     if (llr_type>0) & (opts.if_sub_flip_all | opts.if_sub_flip_any)
%                         llr_lim=opts.plotrange_llr_sub;
%                     end
%                     switch opts.if_plot3d_h
%                         case 0
%                             legend(legh,legt,'FontSize',7,'Location','SouthEast','Interpreter','none');
%                             ylabel(cat(2,'llr',' ',llr_label_suffix));
%                             set(gca,'YTick',[llr_lim(1):.5:llr_lim(2)]);
%                             set(gca,'YLim',llr_lim);
%                         case 1
%                             legend(legh,legt,'FontSize',7,'Location','SouthWest','Interpreter','none');
%                             zlabel(cat(2,'llr',' ',llr_label_suffix));
%                             set(gca,'ZTick',[llr_lim(1):.5:llr_lim(2)]);
%                             set(gca,'ZLim',llr_lim);
%                             ylabel('h');
%                             set(gca,'YTick',[opts.plotrange_h(1):.1:opts.plotrange_h(2)]);
%                             set(gca,'YLim',opts.plotrange_h);
%                             box on;
%                             grid on;
%                             set(gca,'View',opts.view3d);
%                     end %if_plot3d_h
%                     if (isub>1)
%                         legend off;
%                     end
%                 end
%             end %anything to plot?
%         end %llr_type
%     end %ipchoice
%     axes('Position',[0.01,0.01,0.01,0.01]); %for text
%     text(0,0,tstring,'Interpreter','none','FontSize',10);
%     axis off;
%ifk_ptr
% 
% function btcselb_like_plot(alist,llr_list,plot_symb,h,d23)
% %utility plotting: 2d or 3d
% switch d23
%     case 2
%         hp=plot(alist,llr_list,plot_symb);
%     case 3
%         hp=plot3(alist,repmat(h,length(alist),1),llr_list,plot_symb);
% end
% hold on;
return

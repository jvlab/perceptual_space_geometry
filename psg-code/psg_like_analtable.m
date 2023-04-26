function [opts_used,fighs,res]=psg_like_analtable(table_like,opts)
% [opts_used,fighs,res]=psg_like_analtable(table_like,opts) analyzes a table of likelihood ratio data
%
% table_like: likelihood table, typically created by psg_like_maketable.
% opts: options -- all can be omitted and if omitted, will be requested
%  paradigm_type_choice: : 1->animals 2->btc
%  thr_type_choice; threshold type (min, max, avg)
%  frac_keep_choices: 1-> first value of frac keep, 2-> second value etc,
%
% opts_used: options used
% fighs: figure handles
% res: analysis results (not used at present)
%
%
% some categorical variables are converted to tokens, as defined in the
% 'tokens' structure and kept as UserData in table_like.
%
% llr quantities for umi are corrected, i.e., have log(h) subtracted
%
% 09Apr23: convert from script to function, option to plot h as third dimension
% 24Apr23: add alternate terms for intermediate animal paradigms; filled in and empty symbols; option to not plot surrogates
% 25Apr23: option to plot a priori; options for a range and h range
%
%   See also:  PSG_UMI_TRIPLIKE_DEMO, PSG_TENTLIKE_DEMO, PSG_UMI_TRIP_LIKE_RUN, PSG_LIKE_MAKETABLE, PSG_COLORS_LIKEc.
%
res=[];
%
if (nargin<1)
    table_like=table();
end
if (nargin<2)
    opts=struct;
end
opts=filldefault(opts,'fn_table_def','psg_like_maketable_06Apr23.mat');
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
opts=filldefault(opts,'if_surrogates',1);
opts=filldefault(opts,'if_apriori',1);
paradigm_colors=opts.paradigm_colors;
subj_symbs_res=opts.subj_symbs_res;
subj_fills_res=opts.subj_fills_res;
subj_symbs_unres=opts.subj_symbs_unres;
subj_fills_unres=opts.subj_fills_unres;
apriori_symb=opts.apriori_symb;
%
%plot formatting
opts=filldefault(opts,'box_halfwidth',0.02); %half-width of boxes for s.d. of surrogates
opts=filldefault(opts,'plotrange_a',[0 1.25]);
opts=filldefault(opts,'plotrange_llr',[-2 .1]);
opts=filldefault(opts,'if_plot_ah_llr',1); %1 to also plot log likelihood of Dirichlet fit
opts=filldefault(opts,'if_plot3d_h',0); %1 to plot h as third axis 
opts=filldefault(opts,'view3d',[-62 13]);
opts=filldefault(opts,'plotrange_h',[0 .2]);
%
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
tokens=table_like.Properties.UserData.tokens;
%
table_selected=table_like;
%
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
%
switch opts.if_plot3d_h
    case 0
        d23=2;
    case 1
        d23=3;
end
%
fighs=[];
opts_used=opts;
%
for ifk_ptr=1:length(frac_keep_list)
    frac_keep=frac_keep_list(ifk_ptr);
    table_fk=table_selected(table_selected.frac_keep==frac_keep,:);
    %
    %make plots
    %
    %each plot has a panel for the three kinds of llr_types, and the two kinds of ipchoice
    %
    tstring=sprintf('%s: threshold based on %s, frac keep %8.6f',paradigm_type,thr_types{thr_type_choice},frac_keep);
    fighs(end+1)=figure;
    set(gcf,'Position',[50 50 1200 850]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',tstring);
    subj_symbs=subj_symbs_res; %list of subject symbols, starting with reserved list
    subj_fills=subj_fills_res;
    for ipchoice=1:length(tokens.ipchoice)
        for llr_type=(1-opts.if_plot_ah_llr):length(tokens.llr_type)
            table_plot=table_fk(intersect(find(table_fk.ipchoice==ipchoice),find(table_fk.llr_type==max(1,llr_type))),:);
            subj_ids=unique(table_plot.subj_id);
            paradigm_names=unique(table_plot.paradigm_name);
            if llr_type==strmatch('umi',tokens.llr_type,'exact')
                llr_label_suffix=' - log(h)';
            else
                llr_label_suffix=' ';
            end
            if size(table_plot,1)>0
                paradigms_shown=[];
                subjs_shown=[];
                isub=llr_type+opts.if_plot_ah_llr+(ipchoice-1)*nc_plot;
                subplot(nr_plot,nc_plot,isub);
                nsubj_unres=0; %number of subjects with unreserved symbols
                legh=[];
                legt=[];
                for ipt=1:size(table_plot,1)
                    data=table_plot(ipt,:);
                    paradigm_name=data.paradigm_name{1};
                    paradigm_color=paradigm_colors.(paradigm_name);
                    %
                    subj_name=data.subj_id{1};
                    if ~isempty(strmatch(subj_name,fieldnames(subj_symbs),'exact'))
                        subj_symb=subj_symbs.(subj_name);
                        subj_fill=subj_fills.(subj_name);
                    else
                        nsubj_unres=nsubj_unres+1;
                        subj_symb=subj_symbs_unres(1+mod(nsubj_unres-1,length(subj_symbs_unres)));
                        subj_symbs.(subj_name)=subj_symb;
                        subj_fill=subj_fills_unres(1+mod(nsubj_unres-1,length(subj_symbs_unres)));
                        subj_fills.(subj_name)=subj_fill;
                    end
                    if (llr_type>0)
                        %
                        hp=psg_like_plot(data.a,data.llr_data,cat(2,'k',subj_symb),data.h,d23);
                        set(hp,'Color',paradigm_color);
                        if (subj_fill)
                            set(hp,'MarkerFaceColor',paradigm_color);
                        end
                        %plot surrogates
                        if opts.if_surrogates
                            hs=psg_like_plot(repmat(data.a,1,3),[data.llr_data data.llr_flip_all,data.llr_flip_any],'k',data.h,d23);
                            set(hs,'Color',paradigm_color);
                            hs=psg_like_plot(data.a+opts.box_halfwidth*[-1 1 1 -1 -1],data.llr_flip_all+data.llr_flip_all_sd*[1 1 -1 -1 1],'k',data.h,d23);
                            set(hs,'Color',paradigm_color);
                            hs=psg_like_plot(data.a+opts.box_halfwidth*[-2 2 2 -2 -2],data.llr_flip_any+data.llr_flip_any_sd*[1 1 -1 -1 1],'k',data.h,d23);
                            set(hs,'Color',paradigm_color);
                        end
                        if opts.if_apriori
                            %plot a priori
                            ha=psg_like_plot(data.a,data.apriori_llr,cat(2,'k',apriori_symb),data.h,d23);
                            set(ha,'Color',paradigm_color);
                        end
                    else
                        hp=psg_like_plot(data.a,data.ah_llr,cat(2,'k',subj_symb),data.h,d23);
                        set(hp,'Color',paradigm_color);
                        if (subj_fill)
                            set(hp,'MarkerFaceColor',paradigm_color);
                        end
                    end
                    if_addleg=0;
                    if isempty(strmatch(paradigm_name,paradigms_shown,'exact'))
                        paradigms_shown=strvcat(paradigms_shown,paradigm_name);
                        if_addleg=1;
                    end
                    if isempty(strmatch(subj_name,subjs_shown,'exact'))
                        subjs_shown=strvcat(subjs_shown,subj_name);
                        if_addleg=1;
                    end
                    if (if_addleg)
                        legt=strvcat(legt,cat(2,paradigm_name,' ',subj_name));
                        legh=[legh;hp];
                    end
                    if (llr_type>0)
                        title(cat(2,tokens.llr_type{llr_type},' ',tokens.ipchoice{ipchoice}));
                    else
                        title(cat(2,'llr for Dirichlet',' ',tokens.ipchoice{ipchoice}));
                    end
                    xlabel('Dirichlet a');
                    set(gca,'XTick',[opts.plotrange_a(1):.25:opts.plotrange_a(2)]);
                    set(gca,'XLim',opts.plotrange_a);
                    switch opts.if_plot3d_h
                        case 0
                            legend(legh,legt,'FontSize',7,'Location','SouthEast','Interpreter','none');
                            ylabel(cat(2,'llr',' ',llr_label_suffix));
                            set(gca,'YTick',[opts.plotrange_llr(1):.5:opts.plotrange_llr(2)]);
                            set(gca,'YLim',opts.plotrange_llr);
                        case 1
                            legend(legh,legt,'FontSize',7,'Location','SouthWest','Interpreter','none');
                            zlabel(cat(2,'llr',' ',llr_label_suffix));
                            set(gca,'ZTick',[opts.plotrange_llr(1):.5:opts.plotrange_llr(2)]);
                            set(gca,'ZLim',opts.plotrange_llr);
                            ylabel('h');
                            set(gca,'YTick',[opts.plotrange_h(1):.1:opts.plotrange_h(2)]);
                            set(gca,'YLim',opts.plotrange_h);
                            box on;
                            grid on;
                            set(gca,'View',opts.view3d);
                    end %if_plot3d_h
                    if (isub>1)
                        legend off;
                    end
                end
            end %anything to plot?
        end %llr_type
    end %ipchoice
    axes('Position',[0.01,0.01,0.01,0.01]); %for text
    text(0,0,tstring,'Interpreter','none','FontSize',10);
    axis off;
end %ifk_ptr

function hp=psg_like_plot(alist,llr_list,plot_symb,h,d23)
%utility plotting: 2d or 3d
switch d23
    case 2
        hp=plot(alist,llr_list,plot_symb);
    case 3
        hp=plot3(alist,repmat(h,length(alist),1),llr_list,plot_symb);
end
hold on;
return

%psg_like_analtable: analyze a table from the automated outputs of psg_[umi_trip|tent]_like_demo.
%
% table is created by psg_like_maketable.
%
% some categorical variables are converted to tokens, as defined in the
% 'tokens' structure and kept as UserData in table_like.
%
% llr quantities for umi are corrected, i.e., have log(h) subtracted
%
%   See also:  PSG_UMI_TRIPLIKE_DEMO, PSG_TENTLIKE_DEMO, PSG_UMI_TRIP_LIKE_RUN, PSG_LIKE_MAKETABLE.
%
if ~exist('fn_table_def') fn_table_def='psg_like_maketable_06Apr23.mat'; end
%
paradigm_colors=struct;
paradigm_colors.texture=     [0.7 0.0 0.0];
paradigm_colors.texture_like=[0.5 0.5 0.0];
paradigm_colors.image_like=  [0.0 0.5 0.0];
paradigm_colors.image=       [0.0 0.0 0.7];
paradigm_colors.word=        [0.4 0.0 0.4];
%
paradigm_colors.bgca3pt=     [1.0 0.0 0.0];
paradigm_colors.bc6pt=       [0.7 0.0 0.7];
paradigm_colors.bcpm3pt=     [0.0 0.0 0.7];
paradigm_colors.bdce3pt=     [0.0 1.0 0.0];
paradigm_colors.tvpm3pt=     [0.6 0.6 0.0];
%
%reserved subject symbols
subj_symbs_res.MC='o';
subj_symbs_res.SAW='s';
subj_symbs_res.BL='+';
subj_symbs_res.ZK='*';
subj_symbs_unres='*dv^<>ph'; %other available symbols
apriori_symb='.';
%
barwidth=0.02;
%
fn_table=getinp('likelihood table file name','s',[],fn_table_def);
load(fn_table);
disp(sprintf('loaded table with %6.0f rows and %3.0f columns.',size(table_like)));
disp(table_like.Properties.UserData.tokens);
disp(table_like.Properties.UserData.notes);
tokens=table_like.Properties.UserData.tokens;
%
table_selected=table_like;
%
%select paradigm type
paradigm_types_avail=unique(table_like{:,'paradigm_type'});
for ip=1:length(paradigm_types_avail)
    disp(sprintf('%2.0f->%s',ip,paradigm_types_avail{ip}));
end
paradigm_type_choice=getinp('choice','d',[1 length(paradigm_types_avail)]);
paradigm_type=paradigm_types_avail{paradigm_type_choice};
table_selected=table_selected(strmatch(paradigm_type,table_selected.paradigm_type,'exact'),:);
%
%select threshold type
thr_types=tokens.thr_type;
for ip=1:length(thr_types)
    disp(sprintf('%2.0f->%s',ip,thr_types{ip}));
end
thr_type_choice=getinp('choice','d',[1 length(thr_types)]);
table_selected=table_selected(table_selected.thr_type==thr_type_choice,:);
%
%select fraction to keep
frac_keeps=flipud(unique(table_selected.frac_keep));
for ip=1:length(frac_keeps)
    disp(sprintf('%2.0f->keep fraction %8.6f',ip,frac_keeps(ip)));
end
frac_keep_choices=getinp('choice(s)','d',[1 length(frac_keeps)]);
frac_keep_list=frac_keeps(frac_keep_choices);
for ifk_ptr=1:length(frac_keep_list)
    frac_keep=frac_keep_list(ifk_ptr);
    table_fk=table_selected(table_selected.frac_keep==frac_keep,:);
    %
    %make plots
    %
    %each plot has a panel for the three kinds of llr_types, and the two kinds of ipchoice
    %
    nc_plot=length(tokens.llr_type);
    nr_plot=length(tokens.ipchoice);
    tstring=sprintf('%s: thr %s, frac keep %8.6f',paradigm_type,thr_types{thr_type_choice},frac_keep);
    figure;
    set(gcf,'Position',[50 50 1200 850]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',tstring);
    subj_symbs=subj_symbs_res; %list of subject symbols, starting with reserved list
    for ipchoice=1:length(tokens.ipchoice)
        for llr_type=1:length(tokens.llr_type);
            table_plot=table_fk(intersect(find(table_fk.ipchoice==ipchoice),find(table_fk.llr_type==llr_type)),:);
            subj_ids=unique(table_plot.subj_id);
            paradigm_names=unique(table_plot.paradigm_name);
            legh=[];
            legt=[];
            if strcmp(tokens.llr_type{llr_type},'umi')
                ylabel_suffix=' - log(h)';
            else
                ylabel_suffix=' ';
            end
            if size(table_plot,1)>0
                paradigms_shown=[];
                subjs_shown=[];
                subplot(nr_plot,nc_plot,llr_type+(ipchoice-1)*nc_plot);
                nsubj_unres=0; %number of subjects with unreserved symbols
                for ipt=1:size(table_plot,1)
                    data=table_plot(ipt,:);
                    paradigm_name=data.paradigm_name{1};
                    paradigm_color=paradigm_colors.(paradigm_name);
                    %
                    subj_name=data.subj_id{1};
                    if ~isempty(strmatch(subj_name,fieldnames(subj_symbs),'exact'))
                        subj_symb=subj_symbs.(subj_name);
                    else
                        nsubj_unres=nsubj_unres+1;
                        subj_symb=subj_symbs_unres(1+mod(nsubj_unres-1,length(subj_symbs_unres)));
                        subj_symbs.(subj_name)=subj_symb;
                    end
                    hp=plot(data.a,data.llr_data,cat(2,'k',subj_symb));
                    set(hp,'Color',paradigm_color);
                    hold on;
                    %plot surrogates
                    hs=plot(repmat(data.a,1,3),[data.llr_data data.llr_flip_all,data.llr_flip_any],'k');
                    set(hs,'Color',paradigm_color);
                    hs=plot(data.a+barwidth*[-1 1 1 -1 -1],data.llr_flip_all+data.llr_flip_all_sd*[1 1 -1 -1 1],'k');
                    set(hs,'Color',paradigm_color);
                    hs=plot(data.a+barwidth*[-2 2 2 -2 -2],data.llr_flip_any+data.llr_flip_any_sd*[1 1 -1 -1 1],'k');
                    set(hs,'Color',paradigm_color);
                    %plot a priori
                    ha=plot(data.a,data.apriori_llr,cat(2,'k',apriori_symb));
                    set(ha,'Color',paradigm_color);
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
                end
                legend(legh,legt,'FontSize',7,'Location','SouthEast');
                title(cat(2,tokens.llr_type{llr_type},' ',tokens.ipchoice{ipchoice}));
                xlabel('Dirichlet a');
                ylabel(cat(2,'llr',' ',ylabel_suffix));
                set(gca,'XLim',[0 1]);
                set(gca,'YLim',[-2 0]);
            end %anything to plot?
        end %llr_type
    end %ipchoice
    axes('Position',[0.01,0.01,0.01,0.01]); %for text
    text(0,0,tstring,'Interpreter','none','FontSize',10);
    axis off;
end %ifk_ptr
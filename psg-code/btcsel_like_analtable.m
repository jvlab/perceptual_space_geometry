function [opts_used,tables_selected,table_stats]=btcsel_like_analtable(table_like,opts)
% [opts_used,tables_selected,table_stats]=btcsel_like_analtable(table_like,opts): analyze results of log likelihoods of btcsel (btc, with selected stimuli)
%  (but will also run on unselected and animal databases)
% table_like: likelihood table, typically created by psg_like_maketable.
% opts: options -- all can be omitted and if omitted, will be requested
%  thr_type_choice; threshold type (min, max, avg)
%  frac_keep_choices: 1-> first value of frac keep, 2-> second value etc,
%  others retained from psg_like_analtable
%
% opts_used: options used
% tables_selected: cell array, tables_selected{illr,k} is table_like, after selecting subjects, threshold type, threshold values,
%   and log likelihood type (1: sym, 2: umi, 3: adt). User data (containing tokens) and variable names are the same as in table_like.
% tables_stats: tables of summary statistics averaged across subjects.
%   For quantities (such as a and h) in which there is no sd available within subject, then the sem is computed in the standard way across subjects
%   For quantitites (such as orig_data llr) in which there is an sd available within subject, then sem is computed based on both the within-subject sd and the across-subject sd
%    as detailed in ../math/MeanVarCombineNotes.pdf.
%  User data (containing tokens) are the same as in table_like and variable names, when possibe, are the same. 
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
%
%create tables with selected keep fraction and llr type
%
tables_selected=cell(length(tokens.llr_type),length(frac_keep_list));
%
for ifk_ptr=1:length(frac_keep_list)
    frac_keep=frac_keep_list(ifk_ptr);
    table_fk=table_selected(table_selected.frac_keep==frac_keep,:);
    for illr=1:length(tokens.llr_type)
        table_fk_llr=table_fk(table_fk.llr_type==illr,:);
        tables_selected{illr,ifk_ptr}=table_fk_llr;
        tables_selected{illr,ifk_ptr}.Properties.UserData.tokens=tokens;
        disp(sprintf('frac keep = %7.5f, llr type: %s',frac_keep,tokens.llr_type{illr}));
    end
end
%
summary_variables={'llr_type','thr_type','ipchoice','paradigm_name','frac_keep','nsubjs',...
    'ntriads','ntriads_sem',...
    'ntrials','ntrials_sem',...
    'a','a_sem',...
    'h','h_sem',...
    'apriori_llr','apriori_llr_sem',...
    'orig_data','orig_data_sem',...
    'flip_all','flip_all_sem',...
    'flip_any','flip_any_sem'};
%
%create tables with summary statistics
%
table_stats=array2table(zeros(0,length(summary_variables)));
table_stats.Properties.VariableNames=summary_variables;
table_stats.Properties.UserData.tokens=tokens;
%
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
            table_meta=[illr,opts.thr_type_choice,ip_list(ipp)];
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
                            disp('keep frac  nsubjs   ntriads(sem)  ntrials(sem)     a (sem)           h (sem)     apriori_llr(sem)  orig_data_llr(sem)  flip_all_llr(sem)  flip_any_llr(sem)');
                            if_header=1;
                        end
                        %
                        %quantities like orig_data have standard devs estimated for each subject; they are combined ,
                        %as per .../math/MeanVarCombineNotes.pdf.  This variance has two terms: 
                        %within subject: the mean of the variances from each subject (found by squaring their sd's)
                        %across subject: the variance of the subjects' means, found by var with an "N" normalization
                        %Note that for one subject, the within-subject term is the only contributor and the across-subject term is zero, as it should be
                        %After the variance is computed, the sem is computed by sqrt(var)/nsubjs
                        %
                        table_row_vals=[frac_keep,nsubjs,...
                            mean(table_analyze{:,'ntriads'}),std(table_analyze{:,'ntriads'})/sqrt(nsubjs),...
                            mean(table_analyze{:,'ntrials'}),std(table_analyze{:,'ntrials'})/sqrt(nsubjs),...
                            mean(table_analyze{:,'a'}),std(table_analyze{:,'a'})/sqrt(nsubjs),...
                            mean(table_analyze{:,'h'}),std(table_analyze{:,'h'})/sqrt(nsubjs),...
                            mean(table_analyze{:,'apriori_llr'}),std(table_analyze{:,'apriori_llr'})/sqrt(nsubjs),...
                            mean(table_analyze{:,'orig_data'}),sqrt(mean(table_analyze{:,'orig_data_sd'}.^2)+var(table_analyze{:,'orig_data'},0))/sqrt(nsubjs),...
                            mean(table_analyze{:,'flip_all'}),sqrt(mean(table_analyze{:,'flip_all_sd'}.^2)+var(table_analyze{:,'flip_all'},0))/sqrt(nsubjs),...
                            mean(table_analyze{:,'flip_any'}),sqrt(mean(table_analyze{:,'flip_any_sd'}.^2)+var(table_analyze{:,'flip_any'},0))/sqrt(nsubjs)];
                        disp(sprintf('%7.5f  %5.0f   %7.1f(%4.2f) %7.1f(%4.2f)  %6.4f (%6.4f)   %6.4f (%6.4f)   %6.4f (%6.4f)   %6.4f (%6.4f)   %6.4f (%6.4f)   %6.4f (%6.4f)',table_row_vals));
                        %
                        table_strings=array2table({paradigm_names_avail{ipn}});
                        lastrow=[array2table(table_meta) table_strings array2table(table_row_vals)];
                        lastrow.Properties.VariableNames=summary_variables;
                        table_stats=[table_stats;lastrow];
                        %table_stats(end,'paradigm_name')=paradigm_names_avail{ipn};
                    end
                end
            end %ifk_ptr           
         end %ipp
    end %ipn
end %illr
return

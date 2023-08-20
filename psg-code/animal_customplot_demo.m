%animal_customplot_demo: demonstrate custom plotting of btcsel data from tables
% animal stimulus sets
%
%   See also: PSG_LIKE_ANALTABLE, BTCSEL_CUSTOMPLOT_DEMO2.
%
dirichlet_range=[0.0 2.0];
%
if ~exist('opts_plot_def')
    opts_plot_def=struct;
end
%
opts_plot_def=filldefault(opts_plot_def,'paradigm_type_choice',[1:5]);
opts_plot_def=filldefault(opts_plot_def,'thr_type_choice',1); %min
opts_plot_def=filldefault(opts_plot_def,'frac_keep_choices',[1 5]); %fraction to keep
opts_plot_def=filldefault(opts_plot_def,'box_halfwidth',0.02*diff(dirichlet_range)/1.25); %rescale box half-width to match abscissa
opts_plot_def=filldefault(opts_plot_def,'abscissa_alt',0);
opts_plot_def=filldefault(opts_plot_def,'abscissa_para_space',1.8);
%
%read face psg table and get subject count
paradigm_type_all='animal';
table_all=getfield(load('psg_like_maketable_animals_18Jun23.mat'),'table_like');
opts_plot_def=filldefault(opts_plot_def,'subj_id_choice',[1:length(unique(table_all.subj_id))]); %all subjects
%
%suffix_afixed='_a05';
%rows_afixed=find(contains(table_all.paradigm_name,suffix_afixed));
rows_afixed=[]; %no computations in this dataset with a fixed
%
pgroups=cell(0);
pgroups{1}.title='all';
pgroups{1}.list_base='';
pgroups{1}.list_suff={'image','intermediate_object','intermediate_texture','texture','word'};
pgroups{1}.names=pgroups{1}.list_suff; %no renaming
pgroups{1}.abscissa_para_order=[2 3 4 5 1];
%
for if_afixed=-1:0 %[-1:a not fixed, but plot as function of paradigm; 1: a fixed, plot as function of paradigm
    switch if_afixed
        case {-1,0}
            table_sel=table_all(setdiff([1:size(table_all,1)],rows_afixed),:);
            suffix='';
        case 1
            table_sel=table_all(rows_afixed,:);
            suffix=suffix_afixed;
    end
    disp(' ');
    for ig=1:length(pgroups)
        disp(' ');
        disp(sprintf('plotting group %s for if_afixed=%2.0f',pgroups{ig}.title,if_afixed));
        table_plot=table;
        opts_plot=opts_plot_def;
        table_plot_unchanged=table; %for debugging
    	for ip=1:length(pgroups{ig}.list_suff)
            para=cat(2,pgroups{ig}.list_base,pgroups{ig}.list_suff{ip},suffix);
            disp(sprintf('paradigm %s renamed to %s',para,pgroups{ig}.names{ip}));
            table_add=table_sel(strmatch(para,table_sel.paradigm_name,'exact'),:);
            table_plot_unchanged=[table_plot_unchanged;table_add]; %add selected data without changing paradigm name
            table_add.paradigm_name(:)={pgroups{ig}.names{ip}};
            table_add.paradigm_type(:)={paradigm_type_all};
            table_plot=[table_plot;table_add];
        end
        disp(sprintf('unique paradigms included:'))
        disp(unique(table_plot_unchanged.paradigm_name));
        if if_afixed~=0
            opts_plot.abscissa_alt=1;
            opts_plot.abscissa_para_order=pgroups{ig}.abscissa_para_order;
        end
        [ou,fh,res]=psg_like_analtable(table_plot,opts_plot);
        %customize individual plots
        ch=get(gcf,'Children');
        titles=get(ch,'Title');
        tstrings=cell(0);for k=1:length(titles),tstrings{k}=titles{k}.String;end %get title strings
%        umi_fixed=strmatch('umi fixed h ',tstrings,'exact');
%        axes(ch(umi_fixed));
        for ich=1:length(titles)
            if strcmp(get(ch(ich),'Type'),'axes')
                axes(ch(ich))
                if opts_plot.abscissa_alt==0
                    set(gca,'XLim',dirichlet_range);
                    set(gca,'XTick',[min(dirichlet_range):0.2:max(dirichlet_range)]);
                end
                set(gca,'YLim',[-1.0 0.25]);
                set(gca,'YTick',[-1.0:0.25:0.25]);
            end
        end
        %
        fig_name=cat(2,pgroups{ig}.title,' ',suffix);
        set(gcf,'Name',fig_name);
        axes('Position',[0.01,0.05,0.01,0.01]); %for text
        text(0,0,fig_name,'Interpreter','none','FontSize',10);
        axis off;
    end %ig
end %if_afixed

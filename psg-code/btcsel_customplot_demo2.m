%btcsel_customplot_demo2: demonstrate custom plotting of btcsel data from tables
% from bc6 stimulus set
%
%   See also: PSG_LIKE_ANALTABLE, BTCSEL_CUSTOMPLOT_DEMO, FACES_CUSTOMPLOT_DEMO, PSG_TYPENAMES2COLORS,
%    BTCSEL_CUSTOMPLOT_DEMO3.
%
dirichlet_range=[0.0 1.0];
%
if ~exist('opts_plot_def')
    opts_plot_def=struct;
end
%override symbol choice: all filled, differ by shape
subj_fills_res=struct;
subj_fills_res.BL=1;
subj_fills_res.MC=1;
subj_fills_res.SAW=1;
subj_fills_res.ZK=1;
subj_symbs_res=struct;
subj_symbs_res.BL='^';
subj_symbs_res.MC='s';
subj_symbs_res.SAW='v';
subj_symbs_res.ZK='o';
%
opts_plot_def=filldefault(opts_plot_def,'paradigm_type_choice',1);
opts_plot_def=filldefault(opts_plot_def,'thr_type_choice',1); %min
opts_plot_def=filldefault(opts_plot_def,'frac_keep_choices',1); %fraction to keep
opts_plot_def=filldefault(opts_plot_def,'box_halfwidth',0.02*diff(dirichlet_range)/1.25); %rescale box half-width to match abscissa
opts_plot_def=filldefault(opts_plot_def,'abscissa_alt',0);
opts_plot_def=filldefault(opts_plot_def,'subj_fills_res',subj_fills_res);
opts_plot_def=filldefault(opts_plot_def,'subj_symbs_res',subj_symbs_res);
%
%read face psg table and get subject count
paradigm_type_all='btcsel';
table_all=getfield(load('psg_like_maketable_btcsel_14Aug23.mat'),'table_like');
opts_plot_def=filldefault(opts_plot_def,'subj_id_choice',[1:length(unique(table_all.subj_id))]); %all subjects
%
suffix_afixed='_a05';
rows_afixed=find(contains(table_all.paradigm_name,suffix_afixed));
%
colors_def.b=psg_typenames2colors({'bp1000'});
colors_def.c=psg_typenames2colors({'cp1000'});
%
pgroups=cell(0);
pgroups{1}.list_base='bc6pt';
pgroups{1}.title='bc6';
pgroups{1}.list_suff={'_bXcXrand','_bXrand','_cXrand','_bpXcpXrand','_bmXcmXrand'};
pgroups{1}.names={'all','b','c','pos','neg'};
pgroups{1}.colors=[0.00 0.00 0.00;round(colors_def.b);colors_def.c;0.00 0.80 0.00;0.80 0.00 0.00]; %brighter b, default color for c; green for pos, red for neg
pgroups{1}.abscissa_para_order=[1 2 3 5 4];
%
% pgroups{2}.list_base='mpi_en2_fc-';
% pgroups{2}.title='all and by age, M or F';
% pgroups{2}.list_suff={'_','_y_fX_m_f','_y_fX_o_f','_m_fX_o_f','_y_mX_m_m','_y_mX_o_m','_m_mX_o_m'};
% pgroups{2}.names={'all','F_y_m','F_y_o','F_m_o','M_y_m','M_y_o','M_m_o'};
% pgroups{2}.colors=[0.00 0.00 0.00;...
%     1.00 0.50 0.50;0.70 0.35 0.35;0.50 0.30 0.30;...
%     0.50 0.50 1.00;0.35 0.35 0.70; 0.30 0.30 0.50]; %black; reds for F, blues for M
%
for if_afixed=-1:1 %[-1:a not fixed, but plot as function of paradigm; 1: a fixed, plot as function of paradigm
    switch if_afixed
        case {-1,0}
            table_sel=table_all(setdiff([1:size(table_all,1)],rows_afixed),:);
            suffix='';
        case 1
            table_sel=table_all(rows_afixed,:);
            suffix=suffix_afixed;
    end
    disp(' ');
    disp(sprintf('paradigms included in all groups:'));
    disp(unique(table_sel.paradigm_name));
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
            opts_plot.paradigm_colors.(pgroups{ig}.names{ip})=pgroups{ig}.colors(ip,:);
        end
        disp(sprintf('unique paradigms included:'))
        disp(unique(table_plot_unchanged.paradigm_name));
        if if_afixed~=0
            opts_plot.abscissa_alt=1;
            opts_plot.abscissa_para_space=0.6;
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

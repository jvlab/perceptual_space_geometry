%btcsel_customplot_demo: demonstrate custom plotting of btcsel data from tables
%
%   See also: PSG_LIKE_ANALTABLE, FACES_CUSTOMPLOT_DEMO.
%
%merge tables to plot bc6 combined and b and c individual-ray data (7 points each) from bc,
table_all=getfield(load('psg_like_maketable_btc_18Jun23.mat'),'table_like');
table_bc=table_all(strmatch('bc6pt',table_all.paradigm_name),:);
%
table_bc.paradigm_type(:)={'btcsel'}; %so that it all plots in same paradigm
%
table_sel=getfield(load('psg_like_maketable_btcsel_24Jun23.mat'),'table_like');
table_sel_bp=table_sel(strmatch('bc6pt_bp_rand',table_sel.paradigm_name),:);
table_sel_bm=table_sel(strmatch('bc6pt_bm_rand',table_sel.paradigm_name),:);
table_sel_cp=table_sel(strmatch('bc6pt_cp_rand',table_sel.paradigm_name),:);
table_sel_cm=table_sel(strmatch('bc6pt_cm_rand',table_sel.paradigm_name),:);
%
opts.paradigm_colors=struct;
opts.paradigm_colors.bc6pt=[0 0 0];
opts.paradigm_colors.bc6pt_bp_rand=[1.0 0.2 0.2];
opts.paradigm_colors.bc6pt_bm_rand=[0.6 0.0 0.0];
opts.paradigm_colors.bc6pt_cp_rand=[0.2 1.0 0.2];
opts.paradigm_colors.bc6pt_cm_rand=[0.0 0.6 0.0];
table_plot=[table_bc;table_sel_bp;table_sel_bm;table_sel_cp;table_sel_cm];
%
[ou,fh,res]=psg_like_analtable(table_plot,opts);

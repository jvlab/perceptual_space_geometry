%script to run psg_umi_triplike_demo and psg_tentlike on brightness dataset, with selection of specific stimuli
%
%all umi_trip, then all tent
%
%   See also:  PSG_UMI_TRIPLIKE_DEMO, PSG_TENTLIKE_DEMO, PSG_SELECT_CHOICEDATA.
%
%%%%%%%%%%%%%
%GA, umi_trip
%%%%%%%%%%%%%
%
%no selection
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_umi_triplike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='c|s'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'));
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
figs_to_ps_bestfit;
close all
%
%select by c=1
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_umi_triplike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='c01'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'));
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
figs_to_ps_bestfit;
close all
%
%select by c=2
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_umi_triplike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='c02'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'));
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
figs_to_ps_bestfit;
close all
%
%select by s=1368
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_umi_triplike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='s01|s03|s06|s08'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'));
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
figs_to_ps_bestfit;
close all
%
%select by s=1357
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_umi_triplike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='s01|s03|s05|s07'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'));
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
figs_to_ps_bestfit;
close all
%
%select by s=2468
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_umi_triplike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='s02|s04|s06|s08'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'));
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
figs_to_ps_bestfit;
close all
%
%%%%%%%%
%fixed a
%%%%%%%%
%
%fixed a, no selection
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fixa=2;
auto.a_fixval=0.1;
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_umi_triplike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='c|s'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'),'_a01');
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
figs_to_ps_bestfit;
close all
%
%fixed a, select by c=1
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fixa=2;
auto.a_fixval=0.1;
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_umi_triplike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='c01'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'),'_a01');
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
figs_to_ps_bestfit;
close all
%
%fixed a, select by c=2
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fixa=2;
auto.a_fixval=0.1;
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_umi_triplike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='c02'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'),'_a01');
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
figs_to_ps_bestfit;
close all
%
%fixed a, select by s=1368
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fixa=2;
auto.a_fixval=0.1;
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_umi_triplike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='s01|s03|s06|s08'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'),'_a01');
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
figs_to_ps_bestfit;
close all
%
%fixed a, select by s=1357
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fixa=2;
auto.a_fixval=0.1;
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_umi_triplike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='s01|s03|s05|s07'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'),'_a01');
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
figs_to_ps_bestfit;
close all
%
%fixed a, select by s=2468
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fixa=2;
auto.a_fixval=0.1;
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_umi_triplike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='s02|s04|s06|s08'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'),'_a01');
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
figs_to_ps_bestfit;
close all
%
%%%%%%%%%%%%%
%GA, tent
%%%%%%%%%%%%%
%
%no selection
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_tentlike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='c|s'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'));
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
figs_to_ps_bestfit;
close all
%
%select by c=1
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_tentlike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='c01'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'));
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
figs_to_ps_bestfit;
close all
%
%select by c=2
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_tentlike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='c02'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'));
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
figs_to_ps_bestfit;
close all
%
%select by s=1368
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_tentlike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='s01|s03|s06|s08'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'));
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
figs_to_ps_bestfit;
close all
%
%select by s=1357
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_tentlike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='s01|s03|s05|s07'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'));
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
figs_to_ps_bestfit;
close all
%
%select by s=2468
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_tentlike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='s02|s04|s06|s08'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'));
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
figs_to_ps_bestfit;
close all
%
%%%%%%%%
%fixed a
%%%%%%%%
%
%fixed a, no selection
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fixa=2;
auto.a_fixval=0.1;
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_tentlike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='c|s'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'),'_a01');
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
figs_to_ps_bestfit;
close all
%
%fixed a, select by c=1
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fixa=2;
auto.a_fixval=0.1;
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_tentlike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='c01'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'),'_a01');
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
figs_to_ps_bestfit;
close all
%
%fixed a, select by c=2
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fixa=2;
auto.a_fixval=0.1;
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_tentlike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='c02'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'),'_a01');
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
figs_to_ps_bestfit;
close all
%
%fixed a, select by s=1368
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fixa=2;
auto.a_fixval=0.1;
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_tentlike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='s01|s03|s06|s08'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'),'_a01');
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
figs_to_ps_bestfit;
close all
%
%fixed a, select by s=1357
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fixa=2;
auto.a_fixval=0.1;
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_tentlike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='s01|s03|s05|s07'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'),'_a01');
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
figs_to_ps_bestfit;
close all
%
%fixed a, select by s=2468
%
clear
ds_base='bright_c02s08_choices_GA_oddoneout';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fixa=2;
auto.a_fixval=0.1;
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=-1; %do not reorder
auto.db_file='bright_tentlike_db_14Sep23.mat';
auto.sel_apply=1;
auto.sel_string='s02|s04|s06|s08'; %selection
auto.sel_desc=cat(2,'sel_',strrep(auto.sel_string,'|','X'),'_a01');
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
figs_to_ps_bestfit;
close all

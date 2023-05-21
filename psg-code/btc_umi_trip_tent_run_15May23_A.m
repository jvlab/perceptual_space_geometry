%script to run psg_umi_triplike_demo and psg_tentlike
%can create many copies of this, each with a different value of ds_base
%note that the workspace is cleared twice, once for umi_triplike, once for tent
%
%  04Apr23: makes use of autosaving into database in psg_[umi_trip|tent]like_demo
%  24Apr23: incorpores all datasets
%  05May23: incorporates conform surrogate
%  15May23: runs combined across subjects
%
%   See also:  PSG_UMI_TRIPLIKE_DEMO, PSG_TENTLIKE_DEMO.
%
%%%%%
%All subjects to date -- 15 May23 (3 or 4)
%%%%%
clear
ds_base='bgca3pt_choices_A03_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_umi_triplike_db_15May23.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
%
clear
ds_base='bgca3pt_choices_A03_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_tentlike_db_15May23.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
%
figs_to_ps_bestfit;
close all
%%%%%
clear
ds_base='bc6pt_choices_A03_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_umi_triplike_db_15May23.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
%
clear
ds_base='bc6pt_choices_A03_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_tentlike_db_15May23.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
%
figs_to_ps_bestfit;
close all
%%%%%
clear
ds_base='bcpm3pt_choices_A04_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_umi_triplike_db_15May23.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
%
clear
ds_base='bcpm3pt_choices_A04_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_tentlike_db_15May23.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
%
figs_to_ps_bestfit;
close all
%%%%%
clear
ds_base='bdce3pt_choices_A03_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_umi_triplike_db_15May23.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
%
clear
ds_base='bdce3pt_choices_A03_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_tentlike_db_15May23.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
%
figs_to_ps_bestfit;
close all
%%%%%
clear
ds_base='tvpm3pt_choices_A04_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_umi_triplike_db_15May23.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
%
clear
ds_base='tvpm3pt_choices_A04_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_tentlike_db_15May23.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
%
figs_to_ps_bestfit;
close all

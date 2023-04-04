%script to run psg_umi_triplike_demo and psg_tentlike
%can create many copies of this, each with a different value of ds_base
%note that the workspace is cleared and then saved twice, once for umi_triplike, once for tent
%
%this file is for btc data
%
%  04Apr23: makes use of autosaving into database in psg_[umi_trip|tent]like_demo
%
%   See also:  PSG_UMI_TRIPLIKE_DEMO, PSG_TENTLIKE_DEMO.
%
clear
ds_base='bgca3pt_choices_MC_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_umi_triplike_db.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
%save(cat(2,ds_base,'_umi_trip')); data are saved by psg_umi_trip|tent]_like.m
%
clear
ds_base='bgca3pt_choices_MC_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_tentlike_db.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
%save(cat(2,ds_base,'_tent')); ; data are saved by psg_umi_trip|tent]_like.m
%
figs_to_ps_bestfit;
close all
%
clear
ds_base='bc6pt_choices_MC_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_umi_triplike_db.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
%save(cat(2,ds_base,'_umi_trip')); data are saved by psg_umi_trip|tent]_like.m
%
clear
ds_base='bc6pt_choices_MC_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_tentlike_db.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
%save(cat(2,ds_base,'_tent')); ; data are saved by psg_umi_trip|tent]_like.m
%
figs_to_ps_bestfit;
close all
%
clear
ds_base='bcpm3pt_choices_MC_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_umi_triplike_db.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
%save(cat(2,ds_base,'_umi_trip')); data are saved by psg_umi_trip|tent]_like.m
%
clear
ds_base='bcpm3pt_choices_MC_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_tentlike_db.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
%save(cat(2,ds_base,'_tent')); ; data are saved by psg_umi_trip|tent]_like.m
%
figs_to_ps_bestfit;
close all
%
clear
ds_base='bdce3pt_choices_MC_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_umi_triplike_db.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
%save(cat(2,ds_base,'_umi_trip')); data are saved by psg_umi_trip|tent]_like.m
%
clear
ds_base='bdce3pt_choices_MC_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_tentlike_db.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
%save(cat(2,ds_base,'_tent')); ; data are saved by psg_umi_trip|tent]_like.m
%
figs_to_ps_bestfit;
close all
%
clear
ds_base='tvpm3pt_choices_MC_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_umi_triplike_db.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
%save(cat(2,ds_base,'_umi_trip')); data are saved by psg_umi_trip|tent]_like.m
%
clear
ds_base='tvpm3pt_choices_MC_sess01_10';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=1; %reorder for btc
auto.db_file='btc_tentlike_db.mat';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
%save(cat(2,ds_base,'_tent')); ; data are saved by psg_umi_trip|tent]_like.m
%
figs_to_ps_bestfit;
close all
%


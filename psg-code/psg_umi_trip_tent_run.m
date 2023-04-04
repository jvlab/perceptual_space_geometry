%script to run psg_umi_triplike_demo and psg_tentlike
%can create many copies of this, each with a different value of ds_base
%note that the workspace is cleared and then saved twice, once for umi_triplike, once for tent
%
%  04Apr23: makes use of autosaving into database in psg_[umi_trip|tent]like_demo
%
%   See also:  PSG_UMI_TRIPLIKE_DEMO, PSG_TENTLIKE_DEMO.
%
clear
ds_base='SJ_word_choices';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=0; %do not reorder (change this for btc psg expts)
auto.db_file='animals_umi_triplike_db';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_umi_triplike_demo;
%save(cat(2,ds_base,'_umi_trip')); data are saved by psg_umi_trip|tent]_like.m
%
clear
ds_base='SJ_word_choices';
data_fullname=cat(2,'./psg_data/',ds_base);
if_auto=1;
auto=struct; %set non-default options
auto.if_fast=-1; %fast and do not run private
auto.if_reorder=0; %do not reorder (change this for btc psg expts)
auto.db_file='animals_tentlike_db';
h_fixlist=[0 .001 .01 .1]; %only a few forced values of h
psg_tentlike_demo;
%save(cat(2,ds_base,'_tent')); ; data are saved by psg_umi_trip|tent]_like.m
%
figs_to_ps_bestfit;
close all
%


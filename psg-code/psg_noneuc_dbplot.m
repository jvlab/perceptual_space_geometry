%psg_noneuc_dbplot: plot non-euclidean log likelihood fits across subjects, paradigms, etc.
%
% Reads one or more tables created by psg_noneuc_summ, analyzes, and plots
% After running, can also save the concatenated data tables
% 
%  See also: PSG_NONEUC_SUMM, PSG_LLFITS_SUMM, PSG_RAYSTATS_DBPLOT.
%
if ~exist('ui_filter') ui_filter='noneuc_*.mat'; end
if ~exist('ui_filter_gen') ui_filter_gen='noneuc_*_ddmmmyy.mat'; end
%to enable plotting of multiple values against a parameter, within a model dimension
if_multival=1; 
multival_param_source='noneuc_metadata'; %field in Properties.UserData
psg_raystats_dbplot;

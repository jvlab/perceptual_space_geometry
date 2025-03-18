%psg_llfits_dbplot: plot log likelihood fits across subjects, paradigms, etc.
%
% Reads one or more tables created by psg_llfits_summ, analyzes, and plots
% After running, can also save the concatenated data tables
% 
%  See also: PSG_LLFITS_SUMM, PSG_RAYSTATS_DBPLOT.
%
dbtype='llfits';
if_norays=1; %rays are irrelevant
psg_raystats_dbplot;

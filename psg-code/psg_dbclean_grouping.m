%psg_dbclean_grouping: remove spurious entries from databases for
%psg_[raystats|llfits]_dbplot, created by either psg_raystats_summ, psg_llfits_summ, etc.
%
%Spurious entries are those for which standard grouping has sessions 1-10;
%These are remmoved as they are to be supplanted by standard grouping with sessions 1-20
%
%  See also: PSG_RAYSTATS_DBPLOT, PSG_RAYSTATS_SUMM, PSG_LLFITS_SUMM.
%
if ~exist('db_old') db_old='psg_data/llfits_cumulative_13Jan25.mat'; end
if ~exist('db_new') db_new='psg_data/llfits_cumulative_11Mar25.mat'; end
%
disp(sprintf(' loading %s',db_old))
load(db_old);
whos
%
%process t_meta_all
sess_range_meta=cell2mat(t_meta_all{:,'sess_range'});
rows_gp_meta=strmatch('gp',t_meta_all{:,'expt_name'},'exact');
disp('metadata rows with experiment name = gp  and session ranges are:')
disp([rows_gp_meta(:) sess_range_meta(rows_gp_meta,:)]);
rows_delete_meta=intersect(rows_gp_meta,find(sess_range_meta(:,2)==10));
rows_keep_meta=setdiff([1:size(t_meta_all,1)],rows_delete_meta);
t_meta_all_old=t_meta_all;
t_meta_all=t_meta_all_old(rows_keep_meta,:);
disp(sprintf(' %4.0f rows to be deleted from t_meta_all',length(rows_delete_meta)));
%
%process t_all
sess_range=cell2mat(t_all{:,'sess_range'});
rows_gp=strmatch('gp',t_all{:,'expt_name'},'exact');
%disp('metadata rows with experiment name = gp  and session ranges are:')
%disp([rows_gp(:) sess_range(rows_gp,:)]);
rows_delete=intersect(rows_gp,find(sess_range(:,2)==10));
rows_keep=setdiff([1:size(t_all,1)],rows_delete);
t_all_old=t_all;
t_all=t_all_old(rows_keep,:);
disp(sprintf(' %4.0f rows to be deleted from t_all',length(rows_delete)));
if getinp(sprintf('1 if ok to write new t_meta, t_meta_all into %s',db_new),'d',[0 1])
    save(db_new,'t_meta_all','t_all');
end




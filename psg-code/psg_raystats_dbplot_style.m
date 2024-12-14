%psg_raystats_dbplot_style
% utility for psg_raystats_dbplot to set plot style and colors
%
% uses tick_posits,tick_labels,xlims,vars_avail,vars_sel,ivar_sel,if_angle
%
%  See also: PSG_RAYSTATS_DBPLOT.
%

%symbols for each subject, taken when possible from psg_colors_like
if ~exist('subj_model_ID_marker')
    subj_model_ID_marker=struct;
    subj_model_ID_marker.avg='*';
    subj_model_ID_marker.mc='s';
    subj_model_ID_marker.saw='d';
    subj_model_ID_marker.zk='^';
    subj_model_ID_marker.cme='x';
    subj_model_ID_marker.bl='v';
    subj_model_ID_marker.nf='p';
    subj_model_ID_marker.sn='h';
end
if ~exist('expt_grp_colors')
    expt_grp_colors=struct;
    expt_grp_colors.brightness='c';
    expt_grp_colors.constrained_grouping=[0 0.5 0];
    expt_grp_colors.dis_similarity='m';
    expt_grp_colors.similarity='b';
    expt_grp_colors.threshold='k';
    expt_grp_colors.unconstrained_grouping=[0 1 0];
    expt_grp_colors.working_memory='r';
end
if ~exist('expt_uid_linestyle')
    expt_uid_linestyle=struct;
    expt_uid_linestyle.bc6pt=':';
    expt_uid_linestyle.bc55qpt=':';
    expt_uid_linestyle.bcpm3pt=':';
    expt_uid_linestyle.bdce3pt='--';
    expt_uid_linestyle.bgca3pt='-';
    expt_uid_linestyle.dgea3pt='-.';
    expt_uid_linestyle.tvpm3pt=':';
end
%
%check that only one subject, group, and uid is on each plot
%
if size(subj_model_ID,1)~=1
    warning('expecting a single subject ID, found')
    disp(subj_model_ID);
    subj_model_ID=subj_model_ID(1,:);
end
if size(expt_grp,1)~=1
    warning('expecting a single expt grp, found')
    disp(expt_grp);
    expt_grp=expt_grp(1,:);
end
if size(expt_uid,1)~=1
    warning('expecting a single expt grp, found')
    disp(expt_uid);
    expt_uid=expt_uid(1,:);
end
%
%assign symbols, colors, and line types
%
if isfield(subj_model_ID_marker,subj_model_ID)
    set(hp,'Marker',subj_model_ID_marker.(subj_model_ID));                                   
end
if isfield(expt_grp_colors,expt_grp)
    set(hp,'Color',expt_grp_colors.(expt_grp));
end
if isfield(expt_uid_linestyle,expt_uid)
    set(hp,'LineStyle',expt_uid_linestyle.(expt_uid));
end
%capture symbol, color, and line style but remove line if only one point
style=struct;
style.marker=get(hp,'Marker');
style.color=get(hp,'Color');
style.line=get(hp,'LineStyle');
if length(get(hp,'XData'))<=1
    set(hp,'LineStyle','none');
end

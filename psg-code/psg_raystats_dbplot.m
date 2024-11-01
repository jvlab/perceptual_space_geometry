%psg_raystats_dbplot: plot ray statstics across subjects, paradigms, etc.
%
% Reads one or more tables created by psg_raystats_summ, analyzes, and plots
%
%  See also: PSG_RAYSTATS_SUMM.
%
%
if ~exist('ui_filter') ui_filter='raystats_*.mat'; end
%read one or ore data table and keep originals
%optionally apply relabelling of model id ('avg-bl' -> 'bl','avg-zk'->'zk')
%
%table definitions for t_meta_all
meta_variable_names={'psy_model','subj_model_ID','expt_grp','expt_name','expt_param','expt_uid','coords_source','coords_file','choices_source','choices_file','sess_range'};
%
% fields that do not change within a dataset
% psy_model: {'psy'|'mdl'}
% subj_model_ID: 'mc', ..., for a subject; 'qfm' for a generic quadratic model, qfm_mc for a customized model
% expt_grp: {'threshold,'similarity','dis_similarity','working_memory','constrained_grouping','unconstrained_grouping','brightness'};
% expt_name: {'','','dis','wm','gp','gm',br'}; 
% expt_param: NaN for anything but wm, 1000 for wm1000
% expt_uid: 'bc6pt','bc55qpt','bcpm3pt','bgca3pt',etc.
% coords_source: if psychophysical datas: coordinate file name with path, if qform model, a label that has file name and params
% coords_file: coordinate file name
% choices_source: choice file name with path
% choices_file: choice filename
% sess_range: range of session numbers [lo, hi]
%
%table definitions for t_data and t_all
data_variable_names={'dim','jit_crit','if_bid','neg_pos_bid','ray_num','ray2_num','ray_label','ray2_label','var_name','value_data','value_eblo','value_ebhi','value_sem'};
stats_needed={'data','clo','chi','sem'}; %correspondence between value_[data|eblo|ebhi|sem] and fields of angles_stats, mults_stats
% dim: dimension of model
% jit_crit: rms jitter per coordinate for requested pval (depends on dimension but not on the geometry param)
% if_bid: 0 for a measurement from a (unidirectional) ray, 1 for a bidirectional ray
% neg_pos_bid: m for neg, p for pos, z for bidirectional
neg_pos_bid_text={'neg[-]','pos[+]','bid[-+]'};
%ray_num: ray number (for mult), 
%ray2_no: second ray number (for angle), 0 for mult
%ray_label: ray label (endpoint, eg., g0.40)
%ray2_label: ray lebel for second ray (for angle), '' for mult
%var_name: variable name
%value: value of the variable
%value_eblo, value_ebhi, value_sem: lower and upper confidence limit and
%   standard error of measurement, NaN statistics not calculates
%
underscore='_';
dash='-';
%
criterion_names={'subj_model_ID','expt_grp','expt_uid'}; %ways to group or plot
%
%read and concatenate one or more database files, checking for duplicates
%
if_ok=0;
while (if_ok==0)
    ntabs=0;
    t_meta_all_orig=cell(0);
    t_all_orig=cell(0);
    files_orig=cell(0);
    if_more=1;
    while (if_more==1)
        files_req=cell(0);
        n_add_sign=getinp('>0 to add one or more tables by name, <=0 to add one or more from dialog box, 0 if done','d',[-100 100]);
        if (n_add_sign>0) %from console
            n_toadd=abs(n_add_sign);
            files_req=cell(1,n_toadd);
            for iadd=1:abs(n_add_sign)
                files_req{iadd}=getinp(sprintf('table file name, typically %s',ui_filter),'s',[]);
            end
        else %dialog box with logic for specified or unspecified number of files
            [filenames_short,pathname,filter_index]=uigetfile(ui_filter,'select table files','Multiselect','on');            %dialog box here
            if ~iscell(filenames_short) filenames_short={filenames_short}; end
            if (n_add_sign<0) %specified number of files
                if length(filenames_short)~=abs(n_add_sign)
                    n_toadd=0; %will force the dialog box to be ignored
                    disp(sprintf('no files added; %3.0f requested, %3.0f provided',abs(n_add_sign),length(filenames_short)));
                else
                    n_toadd=abs(n_add_sign);
                end
            elseif filter_index==0
                n_toadd=0;
            else
                n_toadd=length(filenames_short);
            end
            for ireq=1:n_toadd
                files_req{ireq}=cat(2,pathname,filesep,filenames_short{ireq});
            end
        end
        for ireq=1:n_toadd
            files_req{ireq}=strrep(files_req{ireq},'\','/');
            if ~exist(files_req{ireq},'file')
                disp(sprintf('%s not found',files_req{ireq}));
            elseif ~isempty(strmatch(files_req{ireq},files_orig,'exact'))
                disp(sprintf('%s already added.',files_req{ireq}));
            else
                s=load(files_req{ireq});
                if ~isfield(s,'t_meta_all') | ~isfield(s,'t_all')
                    disp(sprintf('%s is missing either t_meta_all or t_all',files_req{ireq}));
                else
                    ntabs=ntabs+1;
                    files_orig{ntabs}=files_req{ireq};
                    t_meta_all_orig{ntabs}=s.t_meta_all;
                    t_all_orig{ntabs}=s.t_all;
                    if (ntabs==1)
                        t_meta_all=s.t_meta_all;
                        t_all=s.t_all;
                    else
                        t_meta_all=[t_meta_all;s.t_meta_all];
                        t_all=[t_all;s.t_all];
                    end
                    disp(sprintf('%s added (%4.0f metadata rows, %6.0f data rows), total so far: %4.0f metadata rows, %6.0f data rows.',...
                        files_req{ireq},size(t_meta_all_orig{ntabs},1),size(t_all_orig{ntabs},1),size(t_meta_all,1),size(t_all,1)));
                end %file found?
            end %fields found?
        end %ireq
        disp(sprintf('So far, %3.0f table files concatenated',ntabs));
        if ntabs>0
            if_more=getinp('1 for more files','d',[0 1]);
        end
    end %if_more
    t_meta_all.Properties.UserData.files_orig=files_orig;
    t_all.Properties.UserData.files_orig=files_orig;
    if_ok=getinp('1 if ok','d',[0 1]);   
end
%
%criterion_names={'subj_model_ID','expt_grp','expt_uid'}; %ways to group or plot
%
unique_crits=struct;
for icrit=1:length(criterion_names)
    crit=criterion_names{icrit};
    meta_avail=t_meta_all.(crit);
    data_avail=t_meta_all.(crit);
    unique_crits.(crit).values=unique([meta_avail;data_avail]);
    disp(sprintf(' for %s',crit));   
    for k=1:length(unique_crits.(crit).values) %point to all rows in metadata and data that have this value of the criterion      
        unique_crits.(crit).meta_pointers{k,1}=strmatch(unique_crits.(crit).values{k},t_meta_all.(crit),'exact');
        unique_crits.(crit).meta_counts(1,k)=length(unique_crits.(crit).meta_pointers{k});
        unique_crits.(crit).data_pointers{k,1}=strmatch(unique_crits.(crit).values{k},t_all.(crit),'exact');
        unique_crits.(crit).data_counts(1,k)=length(unique_crits.(crit).data_pointers{k});
        disp(sprintf('     %25s: %5.0f in metadata, %5.0f in data',...
            unique_crits.(crit).values{k},unique_crits.(crit).meta_counts(1,k),unique_crits.(crit).data_counts(1,k)));
    end   
    disp(sprintf('   total:                      %6.0f             %6.0f',...
        sum(unique_crits.(crit).meta_counts),sum(unique_crits.(crit).data_counts)))
end


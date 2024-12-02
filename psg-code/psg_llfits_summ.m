%psg_llfits_summ: summarize log likelihood fits
% use psg_raystats_dbplot for plotting summaries from the tables generated here
%
% For multiple datasets:
% * reads coord files and extracts log likelihood statistics
% * tabulates, for multiple datasets, the log likelihoods of the fits
%
% In comparison with psg_raystats_summ:
%   * does not read choice files
%   * does not depend on jitters
%   * cannot be used on qform model data (no log likelihood fit)
%
% Can set opts_read.ui_filter to filter the data files, e.g., 'bc6*coords*0.mat'
% 
% t_all: a table with all metadata and data, a single line for each
%   model dimension, axis (or axis pair) and each variable measured
%   Overall analysis parameters are in t_all.Properties.UserData
% t_meta_all: a table with the key metadata from each dataset, one line per
%  dataset -- this is duplicated in the left-most columns of t_all
% Overall analysis parameters (pval, etc.) are in t_[meta_]all.Properties.UserData
% t_meta_set(iset): metadata for a single dataset (one line)
%
% When possible, table column headers are compatible with those of [mtc|ramp]_mgm_maketables,
%  which makes mtc_mgm_ramp_tables.mat, with fields including 
%  psy_model, subj_model_ID, expt_grp, expt_name, expt_uid, [varname], [varname]_eblo, [varname]_ebhi
%
% Sample use of the table: extract data from subject ms and dimension 3:
% t_all(intersect(strmatch('mc',t_all.subj_model_ID,'exact'),find(cell2mat(t_all.dim)==3)),:)
% 
%  See also: PSG_GET_COORDSETS, PSG_READ_COORDDATA, PSG_RAYSTATS_SUMM,
%  PSG_DEFOPTS, BTC_DEFINE, PSG_PARSE_FILENAME,
%  MTC_MGM_MAKETABLES, RAMP_MGM_MAKETABLES, PSG_RAYSTATS_DBPLOT.
%
%
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays and psg_get_coordsets
if ~exist('opts_qpred') opts_qpred=struct(); end %for psg_qformpred
if ~exist('data_fullname') data_fullname=[]; end
if ~exist('setup_fullname') setup_fullname=[]; end
opts_read.input_type=1; %qform model not allowed
%
tag_coords='_coords_';
tag_choices='_choices_';
%
aux_fields={'bestModelLL','biasEstimate','debiasedRelativeLL','rawLLs','metadata'}; %required fields
aux_fields_scalar={'bestModelLL'};
aux_fields_bydim={'bestModelLL','biasEstimate','debiasedRelativeLL','rawLLs'};
%
%table definitions for t_meta_all
meta_variable_names={'psy_model','subj_model_ID','expt_grp','expt_name','expt_param','expt_uid','coords_source','coords_file','choices_source','choices_file','sess_range'};
expt_grps=struct;
expt_grps.dis='dis_similarity';
expt_grps.wm='working_memory';
expt_grps.gp='constrained_grouping';
expt_grps.gm='unconstrained_grouping';
expt_grps.br='brightness';
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
%table definitions for t_data and t_all (many for compatibility with psg_raystats_summ)
data_variable_names={'dim','jit_crit','if_bid','neg_pos_bid','ray_num','ray2_num','ray_label','ray2_label','var_name','value_data','value_eblo','value_ebhi','value_sem'};
stats_needed={'data','clo','chi','sem'}; %correspondence between value_[data|eblo|ebhi|sem] and fields of angles_stats, mults_stats
% dim: dimension of model
% jit_crit: NaN
% if_bid: 1
% neg_pos_bid: z
%ray_num: 0 
%ray2_num: 0
%ray_label: ''
%ray2_label: '' for mult
%var_name: variable name
%value: value of the variable
%value_eblo, value_ebhi, value_sem: lower and upper confidence limit and standard error of measurement, NaN if statistics not calculated
%
underscore='_';
dash='-';
zstring='0.00'; %a string to remove from labels, along with leading char
%
opts_read=filldefault(opts_read,'if_log',1);
[sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred);
nsets=length(ds);
%
no_aux=cell(0);
no_aux_set=zeros(0);
%
t_meta_set=cell(nsets,1);
if_t_all=0; %set to zero to initialize table 
for iset=1:nsets
    nrays=rayss{iset}.nrays;
    disp(' ');
    disp(sprintf('set %1.0f: %s',iset,sets{iset}.label))
    nstims=size(ds{iset}{1},1);
    ndims=length(ds{iset});
    dim_list=sets{iset}.dim_list;
    model_dim_max=max(dim_list);
    %
    coords_source=opts_read_used{iset}.data_fullname;
    if isempty(coords_source)
        coords_source='';
    end
    disp(sprintf('processing file: %s',strrep(coords_source,'\','/')));
    if_havechoice=0;
    choices_source='';       
    %
    %set up metadata: parse file full name file name, subject ID, paradigm, etc
    %
    psy_model='unknown';
    subj_model_ID='unknown';
    expt_uid='unknown';
    expt_name='unknown';
    expt_grp='unknown';
    expt_param=NaN;
    sess_range=NaN(1,2);
    coords_file='';
    choices_file='';
    switch sets{iset}.type
        case 'data' % psychophysical data files
            psy_model='psy';
            %
            coords_namestart=max([0,max(find(coords_source=='/')),max(find(coords_source=='\'))]);
            coords_file=coords_source((1+coords_namestart):end);
            choices_namestart=max([0,max(find(choices_source=='/')),max(find(coords_source=='\'))]);
            choices_file=choices_source((1+choices_namestart):end);
            %
            %coords_file is like bc6pt_coords_CME-wm1000_sess01_10.mat or bc6pt_coords_BL_sess01_10.mat
            underscores=find(coords_file==underscore);
            sess_range=zeros(0,2);
            if length(underscores)<4
                warning(sprintf('file name %s cannot be parsed.',coords_file));
            else
                parsed=psg_parse_filename(coords_file);
                expt_uid=parsed.paradigm_name;
                subj_expt=parsed.subj_id; %break into subj_id, epxt_name
                subj_expt_dash=find(subj_expt==dash);
                if ~isempty(subj_expt_dash)
                    subj_model_ID=lower(subj_expt(1:subj_expt_dash(1)-1));
                    expt_name_full=subj_expt([subj_expt_dash(1)+1]:end);
                    if isempty(expt_name_full)
                        expt_name_full='';
                    end
                    expt_name_num=min(regexp(expt_name_full,'[0-9]'));
                    if ~isempty(expt_name_num)
                        expt_name=expt_name_full(1:expt_name_num-1);
                        expt_param=str2num(expt_name_full(expt_name_num:end));
                    else
                        expt_name=expt_name_full;
                        expt_param=NaN;
                    end
                else
                    subj_model_ID=lower(subj_expt);
                    expt_name='';
                    expt_grp='similarity'; %would be threshold if it is a model dataset
                    expt_param=NaN;
                end
                if ~isempty(expt_name)
                    if isfield(expt_grps,expt_name);
                        expt_grp=expt_grps.(expt_name);
                    else
                        expt_grp='unknown';
                    end
                end
                sess_range(1)=str2num(strrep(coords_file(underscores(3)+1:underscores(4)-1),'sess',''));
                sess_range(2)=str2num(coords_file(underscores(4)+1:max(regexp(coords_file,'[0-9]'))));
            end %end parsing
         otherwise
            warning(sprintf('set type for set %2.0f (%s) not recognized or not allowed',iset,sets{iset}.type));
            disp(sets{iset})
    end
    metadata_cell={psy_model,subj_model_ID,expt_grp,expt_name,expt_param,expt_uid,coords_source,coords_file,choices_source,choices_file,sess_range};
    t_meta_set{iset}=array2table(metadata_cell);
    t_meta_set{iset}.Properties.VariableNames=meta_variable_names;
    disp(t_meta_set{iset});
    if (iset==1)
        t_meta_all=t_meta_set{iset};
    else
        t_meta_all=[t_meta_all;t_meta_set{iset}];
    end
    %
    %read auxiliary data in this dataset
    %
    if_aux=0;
    aux_metadata=[];
    if exist(coords_source,'file')
        aux=load(coords_source);
        disp('auxiliary data loaded');
        if_aux=1;
        for ifn=1:length(aux_fields)
            if ~isfield(aux,aux_fields{ifn})
                if_aux=0;
                disp(sprintf('field %s missing',aux_fields{ifn}))
            end
        end
        if if_aux
            if isempty(aux_metadata)
                aux_metadata=aux.metadata;
            else
                if ~strcmp(aux_metadata,aux.metadata)
                    if_aux=0;
                    disp('metadata mismatch');
                end
            end
        end
    else
        disp(sprintf('cannot reload coordinate file %s',strrep(coords_source,'\','/')));
    end
    if (if_aux==0)
        no_aux{end+1}=coords_source;
        no_aux_set(end+1)=iset;
    else
        %process the auxiliary data
    end
%aux_fields={'bestModelLL','biasEstimate','debiasedRelativeLL','rawLLs','metadata'}; %required fields
%aux_fields_scalar={'bestModelLL'};
%aux_fields_bydim={={'bestModelLL','biasEstimate','debiasedRelativeLL','rawLLs'};

    % %
    % %for each dimension model, find best-fitting signed and unsigned rays, including the origin
    % %
    % for idimptr=1:length(dim_list) %compute ray fits and angles
    %     idim=dim_list(idimptr);
    %     for if_bid=0:1 %uni- and bi-directional
    %         opts_stats_use=setfield(opts_stats,'if_bid',if_bid);
    %         if if_havechoice==0
    %             opts_stats_use=setfield(opts_stats_use,'nsurrs',0);
    %         end
    %         jit_use=jit_rms_list(idim);
    %         if isnan(jit_use) %will only happen if jitters are all NaN for all dimensions
    %             jit_use=0;
    %         end
    %         [angles{iset,1+if_bid}{idim},mults{iset,1+if_bid}{idim},angles_stats{iset,1+if_bid}{idim},mults_stats{iset,1+if_bid}{idim},rayfit{iset,1+if_bid}{idim}]=...
    %             psg_raystats(ds{iset}{idim},sas{iset},rayss{iset},jit_rms_list(idim),opts_stats_use);
    %         %
    %     end %if_bid
    % end %idim_ptr
    % %
    % %nicely formatted output of mults and angles to console and add to table
    % %
    % disp(' ');
    % disp('gains along each ray, and for bidirectional fit to each axis')
    % disp(cat(2,' dim    ',mult_labels{iset}{:}));
    % for idimptr=1:length(dim_list) %display ray fits and angles
    %     idim=dim_list(idimptr);
    %     stat_names=fieldnames(mults_stats{iset,1}{idim});
    %     nstats=length(stat_names); %will be zero if nsurrs=0 or choice data file not found
    %     v=cell(nstats+1,2);
    %     stats_have=cell(1,nstats+1); %this is to map fields from mults_stats into clo, chi, sem
    %     for iv=0:nstats
    %         if (iv==0)
    %             text_string=sprintf('%3.0f   ',idim);
    %             for ib=1:2
    %                 v{iv+1,ib}=mults{iset,ib}{idim};
    %             end
    %             stats_have{1}='data';
    %         else
    %             stat_name=stat_names{iv};
    %             stats_have{iv+1}=stat_name;
    %             text_string=sprintf('%6s',stat_name);
    %             for ib=1:2
    %                 v{iv+1,ib}=mults_stats{iset,ib}{idim}.(stat_name);
    %             end
    %         end
    %         for iray=1:nrays
    %             text_string=cat(2,text_string,sprintf(' %15.4f',v{iv+1,1}.dist_gain(iray,1)),'     '); %gain on negative ray
    %             text_string=cat(2,text_string,sprintf(' %15.4f',v{iv+1,1}.dist_gain(iray,2)),'     '); %gain on positive ray
    %             text_string=cat(2,text_string,sprintf(' %15.4f',v{iv+1,2}.dist_gain(iray,1)),'     '); %bidirectional gain
    %         end
    %         disp(text_string);
    %     end %iv 0=data, >=1: stats
    %     disp(' ');
    %     %add to table
    %     mult_vals=NaN(4,3,nrays); %d1: value, clo, chi, sem; d2: unipolar neg, unipolar pos, bidir; d3: each ray
    %     for iv=0:nstats
    %         ivptr=strmatch(stats_needed{iv+1},stats_have,'exact'); %look for clo, chi, and sem in output of psg_raystats
    %         if length(ivptr)==1
    %             for iray=1:nrays
    %                 mult_vals(iv+1,1,iray)=v{ivptr,1}.dist_gain(iray,1);
    %                 mult_vals(iv+1,2,iray)=v{ivptr,1}.dist_gain(iray,2);
    %                 mult_vals(iv+1,3,iray)=v{ivptr,2}.dist_gain(iray,1);
    %             end
    %         end
    %     end %iv
    %     vname='dist_gain';
    %     for inpb=1:3 %1: unipolar negative, 2: unipolar positive, 3: bidirectional
    %         for iray=1:nrays
    %             if_bid=ismember(inpb,[3]);
    %             neg_pos_bid=neg_pos_bid_text{inpb};
    %             values=mult_vals(:,inpb,iray)';
    %             data_cell=[{idim,jit_rms_list(idim),if_bid,neg_pos_bid,iray,0,ray_labels{iset}{iray},'',vname} num2cell(values)];
    %             t_data=array2table(data_cell);
    %             t_data.Properties.VariableNames=data_variable_names;
    %             if (if_t_all==0)
    %                 t_all=[t_meta_set{iset},t_data];
    %                 if_t_all=1;
    %             else
    %                 t_all=[t_all;[t_meta_set{iset},t_data]];
    %             end
    %         end %iray
    %     end %inpb 
    % end %idim_ptr
    % disp(' ');
      %     end %iray
    % end %idim_ptr
end %iset
%
disp(sprintf(' auxiliary data not found for %2.0f datasets',length(no_aux)));
for k=1:length(no_aux)
    disp(sprintf('missing set %2.0f (set %2.0f of inputs): %s',k,no_aux_set(k),no_aux{k}));
end
%
disp('adding overall settings to UserData of tables')
settings=struct;
settings.ll_metadata=aux_metadata;
settings_fields=fieldnames(settings);
disp(settings);
for k=1:length(settings_fields)
    if isstruct(settings.(settings_fields{k}))
        disp(settings_fields{k});
        disp(settings.(settings_fields{k}));
    end
end
t_meta_all.Properties.UserData=settings;
t_all.Properties.UserData=settings;
%
disp('suggest saving t_meta_all and t_all in a file such as raystats_*.mat')


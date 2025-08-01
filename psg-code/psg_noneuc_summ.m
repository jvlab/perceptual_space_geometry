%psg_noneuc_summ: summarize non-Euclidean curvature analysis
% use psg_llfits_dbplot->psg_raystats_dbplot for plotting summaries from the tables generated here
%
% For multiple datasets:
% * reads coord files and extracts log likelihood statistics
% * tabulates, for multiple datasets, curvature data
%
% Derived from psg_llfits_summ, and, like psg_llfits_summ, is driven by reading coordinate files.
% This is used for metadata and to determine the name of a csv file with curvature data.
% In contrast to psg_llfits_sum: 
%    * no dimension-independent quantities (aux_fields_nodim={})
%    * dimension list is extracted from csv file into dim_list_curve (typically 2:5), not dim_list (typically 1:7) in coordinate file
%
% List of curvature values in csv file (lambda if hyperbolic, 2*mu if
% spherical) checked for conformance to curve_lambda_mu_list, and is stored in t_meta_all.Properties.UserData
%  If there are missing values or extra values, processing will continue
%  but a warning will be issued
%
% Can set opts_read.ui_filter to filter the data files, e.g., 'bc6*coords*0.mat'
% 
% t_all: a table with all metadata and data, a single line for each
%   model dimension, axis (or axis pair) and each variable measured
%   Overall analysis parameters are in t_all.Properties.UserData
% t_meta_all: a table with the key metadata from each dataset, one line per
%  dataset -- this is duplicated in the left-most columns of t_all
% t_meta_set(iset): metadata for a single dataset (one line)
%
% When possible, table column headers are compatible with those of [mtc|ramp]_mgm_maketables,
%  which makes mtc_mgm_ramp_tables.mat, with fields including 
%  psy_model, subj_model_ID, expt_grp, expt_name, expt_uid, [varname], [varname]_eblo, [varname]_ebhi
%
% Sample use of the table: extract data from subject ms and dimension 3:
% t_all(intersect(strmatch('mc',t_all.subj_model_ID,'exact'),find(cell2mat(t_all.dim)==3)),:)
% 
% 17Jun25: fix for wm1000 -> wm, add directories for grouping and working memory
%
%  See also: PSG_LLFITS_SUMM, PSG_DEFOPTS, BTC_DEFINE, PSG_PARSE_FILENAME,
%  MTC_MGM_MAKETABLES, RAMP_MGM_MAKETABLES, PSG_NONEUC_DBPLOT, PSG_LLFITS_DBPLOT, PSG_RAYSTATS_DBPLOT.
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
tag_curve='_curve_';
%
curv_constant_names={'Domain','Sigma','Subject','Task'}; %these should be constant within a curvature file
%
aux_fields={'rawLLs','biasEstimate','debiasedRelativeLL','debiasedLLs'};
aux_fields_nodim={};
aux_fields_bydim=aux_fields;
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
% dim: dimension of model, or 0 if independent of model dimension
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
%quantities specific to curvature data
%to construct file names with non-Euclidean analyses
if ~exist('curve_paths')
    curve_paths={...
%        '.\psg_data\psg_data_NonEuc',... %this path had preliminary analyses with a different file name and the Task field missing
        '.\psg_data_NonEuc\standard\std',...
        '.\psg_data_NonEuc\standard\combos',...
        '.\psg_data_NonEuc\brightness',...
        '.\psg_data\psg_data_NonEuc\standard\std',...
        '.\psg_data\psg_data_NonEuc\standard\combos',...
        '.\psg_data\psg_data_NonEuc\brightness',...
        '.\psg_data\psg_data_NonEuc\grouping',...
        '.\psg_data\psg_data_NonEuc\workingmemory'};
end
if ~exist('curve_infixs')
    curve_infixs={'','_dim2-5'};
end
if ~exist('curve_ext') curve_ext='.csv'; end
curve_expt_uid_suffix='9'; 
% if these values are not all found, or, if other values are found, a
% warning will be issued but datasets will be processed
%
if ~exist('curve_lambda_mu_list')
    curve_lambda_mu_list_raw=[0:.125:1 1.25:.25:2.5 0:-.25:-2 -2.5:-.5:-5]; %list of values of lambda and mu in typical csv file
    curve_lambda_mu_list=unique(curve_lambda_mu_list_raw); %remove duplicates and sort from low to high
end
curve_nvals=length(curve_lambda_mu_list);
curve_list=zeros(1,curve_nvals);
curve_list(curve_lambda_mu_list<0)=  curve_lambda_mu_list(curve_lambda_mu_list<0);
curve_list(curve_lambda_mu_list>0)=2*curve_lambda_mu_list(curve_lambda_mu_list>0);
%conversion of csv column names to database variable names, corresponding, where possible, to llfits database
curve_csv_headers=struct;
curve_csv_headers.rawLLs='Log Likelihood'; %as in llfits database
curve_csv_headers.biasEstimate='Bias Estimate'; %as in llfits database
curve_csv_headers.debiasedRelativeLL='Corrected LL Relative to Best Model'; %as in llfits database
curve_csv_headers.debiasedLLs='Corrected LLs'; %no corresponding quantity in llfits database; equal to rawLLs+biasEstimate
curve_varnames=fieldnames(curve_csv_headers);
curve_nvars=length(curve_varnames);
%
curve_tol=10^-3; %tolerance for parameters matching expected values; this is forgiving b/o possible rounding errors
%
opts_read=filldefault(opts_read,'if_log',1);
[sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred);
nsets=length(ds);
%
no_aux=cell(0);
no_aux_set=zeros(0);
%
warn_aux=cell(0);
warn_aux_set=zeros(0);
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
    %read curvature data for this dataset
    %
    if_aux=0;
    if_aux_warn=0;
    noneuc_metadata=[];
    sep_last=max([0,max(find(coords_source=='/')),max(find(coords_source=='\'))]);
    coords_source_base=coords_source(1+sep_last:end);
    coords_source_base=strrep(coords_source_base,'.mat','');
    curve_source_base=strrep(coords_source_base,tag_coords,tag_curve);
    curve_source_base=strrep(curve_source_base,'wm1000','wm'); %17Jun25
    curve_file_list=cell(0);
    curve_file_nmatches=0;
    curve_file='';
    if_aux=0;
    for curve_infix=1:length(curve_infixs)
        for curve_path=1:length(curve_paths)
            curve_source=cat(2,curve_paths{curve_path},filesep,curve_source_base,curve_infixs{curve_infix},curve_ext);
%            curve_source
            if exist(curve_source,'file')
                curve_file_nmatches=curve_file_nmatches+1;
                curve_file_list{curve_file_nmatches}=curve_source;
            end
        end
    end
    if curve_file_nmatches>1
        for imatch=1:curve_file_nmatches
            disp(sprintf(' for %30s, match %1.0f is %s',coords_source_base,imatch,curve_file_list{imatch}));
        end
        match_choice=getinp('choice','d',[1 curve_file_nmatches]);
        curve_file=curve_file_list{match_choice};
        if_aux=1;
    elseif curve_file_nmatches==1
        curve_file=curve_file_list{1};
        if_aux=1;
    end
    if if_aux==1
        disp(sprintf(' for %30s, loading curvature data from %s',coords_source_base,curve_file));
        t_curve=readtable(curve_file,'VariableNamingRule','preserve');
    else
        disp(sprintf('cannot find curvature data file for %s',strrep(coords_source,'\','/')));
    end
    %
    %check self-consistency
    %
    if if_aux==1
        for ic=1:length(curv_constant_names)
            if length(unique(t_curve{:,curv_constant_names{ic}}))>1
                disp(sprintf('curvature file inconsistency: %10s varies within the file',curv_constant_names{ic}));
                if_aux=0;
            end
        end
    end
    %
    %if self-consistent, check consistency with coord file
    %
    if (if_aux==1)
        %
        Sigma=t_curve{1,'Sigma'};
        if Sigma~=1
            disp(sprintf('Sigma should be 1, found %5.3f',Sigma));
            if_aux=0;
        end
        %
        subj_model_ID=char(t_meta_set{iset}{1,'subj_model_ID'});
        Subject=char(t_curve{1,'Subject'});
        if ~strcmp(lower(Subject),subj_model_ID)
            disp(sprintf('Subject should be %s, found %s',subj_model_ID,Subject));
            if_aux=0;
        end
        %
        expt_uid=char(t_meta_set{iset}{1,'expt_uid'});
        Domain=char(t_curve{1,'Domain'});
        if ~strcmp(Domain,expt_uid) & ~strcmp(Domain,cat(2,expt_uid,curve_expt_uid_suffix))
            disp(sprintf('Domain should be %s, found %s',expt_uid,Domain));
            if_aux=0;
        end
        %
        expt_grp=char(t_meta_set{iset}{1,'expt_grp'});
        Task=char(t_curve{1,'Task'});
        if ~strcmp(Task,expt_grp) & ~strcmp(strrep(Task,'mem','memory'),expt_grp) %working_mem or working_memory is OK, 17Jun25
            disp(sprintf('Task should be %s, found %s',expt_grp,Task));
            if_aux=0;
        end
    end
    %
    % read and check computed values
    %
    if (if_aux==1)
        dim_list_curve=unique(t_curve{:,'Dimension'});
        ndim_list_curve=length(dim_list_curve);
        lambda_mu_vals=cell(ndim_list_curve,1);
        curve_vals=cell(ndim_list_curve,1);
        ndim_rows=zeros(1,ndim_list_curve);
        for idim_ptr=1:ndim_list_curve
            dim_val=dim_list_curve(idim_ptr);
            dim_rows=find(t_curve{:,'Dimension'}==dim_val);
            ndim_rows(idim_ptr)=length(dim_rows);
            lambda_mu_vals{idim_ptr}=t_curve{dim_rows,'Lambda-Mu'}';
            curve_vals{idim_ptr}=t_curve{dim_rows,'Curvature of Space'}';
        end
        if min(ndim_rows)~=max(ndim_rows)
            disp('Unequal number of values of lambda_mu and curvature across dimensions');
            if_aux=0;
        else %check consistency lf lambda_mu and curvature, across and within dimensions
            for idim_ptr=1:ndim_list_curve
                if any(lambda_mu_vals{idim_ptr}~=lambda_mu_vals{1})
                    disp('Unequal lambda_mu values across dimensions');
                    if_aux=0;
                end
                if any(curve_vals{idim_ptr}~=curve_vals{1})
                    if_aux=0;
                    disp('Unequal curvature values across dimensions');
                end
                if any(curve_vals{idim_ptr}(lambda_mu_vals{idim_ptr}<0)~=lambda_mu_vals{idim_ptr}(lambda_mu_vals{idim_ptr}<0))
                    if_aux=0;
                    disp(sprintf('Inconsistent curvature and lambda_mu values for negative curvature in dimension %2.0f',dim_list_curve(idim_ptr)));
                end
                if any(curve_vals{idim_ptr}(lambda_mu_vals{idim_ptr}>0)~=2*lambda_mu_vals{idim_ptr}(lambda_mu_vals{idim_ptr}>0))
                    if_aux=0;
                    disp(sprintf('Inconsistent curvature and lambda_mu values for positive curvature in dimension %2.0f',dim_list_curve(idim_ptr)));
                end
           end
        end %consistency check for lambda and mu
        %
        %establish pointers from curve_lambda_mu_list into data
        %       
        if (if_aux==1)
            curve_ifused=zeros(1,ndim_rows(1));
            curve_ptrs=cell(1,curve_nvals);
            curve_found=zeros(1,curve_nvals);
            for ival=1:curve_nvals
                matches=find(abs(curve_vals{1}-curve_list(ival))<=curve_tol);
                curve_ptrs{ival}=matches;
                curve_ifused(matches)=curve_ifused(matches)+1;
                if length(matches)>0
                    curve_found(ival)=curve_found(ival)+1;
                end
            end
            if any(curve_ifused>1)
                disp('Multiple matches to a standard curvature value, standard curvature values may be duplicated');
                if_aux=0;
            end
            if any(curve_ifused==0)
                disp('Some available curvature values in csv file not used');
                if_aux_warn=1;
            end
            if any(curve_found==0)
                disp('Some expected curvature values not found in csv file');
                if_aux_warn=1;
            end
        end
    end
    if (if_aux==1)
        %exctract log likelihoods, and check that if the curvature value is
        %the same, that the log likelihoods agree within curve_tol
        curve_varvals=NaN(ndim_list_curve,length(curve_list),curve_nvars); % d1: model dimension, d2: curvature value (sorted neg to pos), d3: raw ll, bias estimate, debiased rel llcorrected ll
        for idim_ptr=1:ndim_list_curve
            dim_val=dim_list_curve(idim_ptr);
            dim_rows=find(t_curve{:,'Dimension'}==dim_val);
            for ivar=1:curve_nvars
                varname=curve_varnames{ivar};
                csv_vals_all=t_curve{dim_rows,curve_csv_headers.(varname)}';
                for icl=1:length(curve_list)
                    csv_ptrs=curve_ptrs{icl};
                    if ~isempty(csv_ptrs)
                        csv_vals=csv_vals_all(curve_ptrs{icl});
                        curve_varvals(idim_ptr,icl,ivar)=mean(csv_vals);
                        if abs(min(csv_vals)-max(csv_vals))>curve_tol
                            disp(sprintf('Mismatch in csv file for curvature value %7.4f on dimension %2.0f of %s',...
                                curve_list(icl),dim_val,varname));
                            if_aux_warn=1;
                        end %mismatch
                    end
                end %icl
            end %ivar
        end %idim_ptr      
    end
    if (if_aux_warn==1)
        warn_aux{end+1}=coords_source;
        warn_aux_set(end+1)=iset;
    end
    if (if_aux==0)
        no_aux{end+1}=coords_source;
        no_aux_set(end+1)=iset;
    else %process the auxiliary data
        % data_cell={idim,jit_crit,if_bid,neg_pos_bid,ray_num,ray2num,ray_label,ray2_label,var_name,value_data,value_eblo,value_ebhi,value_sem};
        %
        %fields that are dimension-independent (none)
        % for ifn=1:length(aux_fields_nodim)
        %     fn=aux_fields_nodim{ifn};
        %     data_cell=[{0,NaN,1,'z',0,0,'','',fn} num2cell(aux.(fn)) NaN NaN NaN];
        %     t_data=array2table(data_cell);
        %     t_data.Properties.VariableNames=data_variable_names;
        %     if (if_t_all==0)
        %         t_all=[t_meta_set{iset},t_data];
        %         if_t_all=1;
        %     else
        %         t_all=[t_all;[t_meta_set{iset},t_data]];
        %     end
        % end %dimension-independent fields
        %
        %fields that depend on dimension
        for ifn=1:length(aux_fields_bydim)
            fn=aux_fields_bydim{ifn};
            for idimptr=1:ndim_list_curve %extract log likelihoods for nonEuclidean fits
                idim=dim_list_curve(idimptr);
                data_cell=[{idim,NaN,1,'z',0,0,'','',fn} mat2cell(curve_varvals(idimptr,:,ifn),1,length(curve_list)) NaN NaN NaN];
                t_data=array2table(data_cell);
                t_data.Properties.VariableNames=data_variable_names;
                if (if_t_all==0)
                    t_all=[t_meta_set{iset},t_data];
                    if_t_all=1;
                else
                    t_all=[t_all;[t_meta_set{iset},t_data]];
                end
            end %dimptr
        end %dimension-independent fields
    end %if_aux
end %iset
%
disp(sprintf(' auxiliary data not found or invalid for %2.0f datasets',length(no_aux)));
for k=1:length(no_aux)
    disp(sprintf('missing or invalid set %2.0f (set %2.0f of inputs): %s',k,no_aux_set(k),no_aux{k}));
end
disp(sprintf(' auxiliary data warning for %2.0f datasets',length(warn_aux)));
for k=1:length(warn_aux)
    disp(sprintf('warning for set %2.0f (set %2.0f of inputs): %s',k,warn_aux_set(k),warn_aux{k}));
end

%
disp('adding overall settings to UserData of tables')
settings=struct;
settings.noneuc_metadata.curve_lambda_mu_list=curve_lambda_mu_list;
settings.noneuc_metadata.curve_list=curve_list;
%
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
disp('suggest saving t_meta_all and t_all in a file such as noneuc_*_ddmmmyy.mat')


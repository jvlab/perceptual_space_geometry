%psg_raystats_summ: summarize ray angle statistics (see psg_visualize_demo for plotting)
%
% For multiple datasets:
% * reads choice files to find critical jitter
% * propagates the confidence limits based on critical jitter to ray
%   geometry
% * tabulates, for multiple datasets
%     angles between positive and negative rays on each axis (using psg_rayfit with bid=1)
%     angles between best-fitting line for each pair of axes (using psg_rayfit with bid=0)
%     multipliers (gains) on each axis
%
% This is designed for datasets that have positive and negative extents on each axis.
% Computation of ray statistics for datasets in one quadrant (bcpp55, bcpm55, bcmp55, bcmm55), or circular (bc24)
%    will generate many warnings. 
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
% t_all(intersect(strmatch('mc',cell2mat(t_all.subj_model_ID),'exact'),find(cell2mat(t_all.dim)==3)),:)
% 
%  See also: PSG_GET_COORDSETS, PSG_READ_COORDDATA, PSG_FINDRAYS,
%  PSG_DEFOPTS, BTC_DEFINE, PSG_PARSE_FILENAME,
%  PSG_RAYFIT, PSG_RAYANGLES, PSG_RAYMULTS, PSG_RAYSTATS, PSG_LLJIT_CRIT, MTC_MGM_MAKETABLES, RAMP_MGM_MAKETABLES.
%
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays and psg_get_coordsets
if ~exist('opts_fit') opts_fit=struct(); end %for psg_rayfit
if ~exist('opts_stats') opts_stats=struct(); end %for psg_raystats
if ~exist('opts_lljit') opts_lljit=struct(); end %for psg_lljit_crit
if ~exist('opts_qpred') opts_qpred=struct(); end %for psg_qformpred
if ~exist('data_fullname') data_fullname=[]; end
if ~exist('setup_fullname') setup_fullname=[]; end
%
njit_types=2; %two types of critical jitter, see psj_lljit_crit
%
tag_coords='_coords_';
tag_choices='_choices_';
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
zstring='0.00'; %a string to remove from labels, along with leading char
%
opts_lljit=filldefault(opts_lljit,'ndraws',100);
opts_lljit.ndraws=getinp('number of draws for critical jitter calculation','d',[1 10^4],opts_lljit.ndraws);
%
opts_stats=filldefault(opts_stats,'nsurrs',100);
opts_stats.nsurrs=getinp('number of surrogates for final confidence limits (0 to omit)','d',[0 10^4],opts_stats.nsurrs);
%
if ~exist('pval') pval=0.05; end
pval=getinp('p-value for confidence limits','f',[0 1],pval); 
%
jit_crit_choice=getinp('choice of critical jitter or 0 to combine','d',[0 njit_types],0);
%
opts_read=filldefault(opts_read,'if_log',1);
[sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred);
nsets=length(ds);
ray_labels=cell(nsets,1);
ray_pair_labels=cell(nsets,1);
mult_labels=cell(nsets,1);
%
angles=cell(nsets,2);
mults=cell(nsets,2);
angles_stats=cell(nsets,2);
mults_stats=cell(nsets,2);
rayfit=cell(nsets,2);
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
    disp(sprintf('processing file: %s',strrep(coords_source,'/','\')));
    %if nsurrs>0, read choice file so that surrogates for error bars can be created.
    %
    if_havechoice=1;
    if contains(coords_source,'_coords_') & opts_stats.nsurrs>0
        choices_source=strrep(coords_source,tag_coords,tag_choices);
        disp(sprintf('   choice file:  %s',strrep(choices_source,'/','\')));
        if exist(choices_source,'file')
            c=load(choices_source);
            disp(sprintf('    nstims: %3.0f dims: %3.0f, cols in responses: %3.0f',nstims,ndims,size(c.responses,2)));
        else
            disp(sprintf('choice file not found: %s',strrep(choices_source,'/','\')));
            if_havechoice=0;
        end
    else
        if_havechoice=0;
        choices_source='';       
    end
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
        case 'qform' %quadratic form model
            psy_model='mdl';
            coords_source=sets{iset}.label_long;
            expt_grp='threshold'; %these are always threshold models
            expt_name='';
            % parse label_long to obtain expt_uid and subject ID (i.e., source of model)
            %   label_long: './psg_data/bgca3pt9.mat c:aug q:../stim/btc_allraysfixedb_avg_100surrs_madj.mat m:qform'
            label_long=sets{iset}.label_long;
            qform_file_ptr=strfind(label_long,'btc_allraysfixedb_');
            qform_mat_ptr=strfind(label_long,'.mat');
            if_qform_ok=1;
            if (length(qform_file_ptr)~=1) | (length(qform_mat_ptr)~=2)
                if_qform_ok=0;
            else
                %extract expt_uid from setup file name
                qform_setup_fullname=label_long(1:qform_mat_ptr(1)-1);
                qform_setup_file_start=max([0,max(find(qform_setup_fullname=='/')),max(find(qform_setup_fullname=='\'))]);
                expt_uid=qform_setup_fullname(1+qform_setup_file_start:end);
                expt_uid_pt=strfind(expt_uid,'pt');
                if ~isempty(expt_uid_pt)
                    expt_uid=expt_uid(1:expt_uid_pt+1); %delete anything after pt
                end
                %extract subject ID that is source of model from qform file name
                qform_desc=label_long(qform_file_ptr:end); %something like 'btc_allraysfixedb_avg_100surrs_madj.mat m:qform';           
                qform_underscores=find(qform_desc==underscore);
                if length(qform_underscores)<3
                    if_qform_ok=0;
                else
                    subj_model_ID=cat(2,'qform',dash,qform_desc(qform_underscores(2)+1:qform_underscores(3)-1));
                end
            end

            if (if_qform_ok==0)
                warning(sprintf('quadratic form model label %s cannot be parsed',label_long));
            end
        otherwise
            warning(sprintf('set type for set %2.0f (%s) not recognized',iset,sets{iset}.type));
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
    %compute critical jitters
    %
    jit_rms_list=NaN(model_dim_max,1);
    jit_crits=NaN(model_dim_max,njit_types); %type 1 and type 2 jits
    if if_havechoice
        for idimptr=1:length(dim_list) %compute ray fits and angles
            idim=dim_list(idimptr);
            jit_crits(idim,:)=psg_lljit_crit(pval,ds{iset}{idim},sas{iset}.typenames,c.responses,c.stim_list,opts_lljit);
        end       
        disp('critical jitters by type')
        for itype=1:njit_types
            jits_nan=find(isnan(jit_crits(:,itype)));
            if ~isempty(jits_nan) %if any of the jits are NaN, replace by maximum of non-nan jits
                jits_nonan=setdiff(dim_list,jits_nan);
                if ~isempty(jits_nonan)
                    jit_max=max(jit_crits(jits_nonan,itype));
                    jit_crits(jits_nan,itype)=jit_max;
                end
                disp(cat(2,' critical jitters were NaN for model dimensions:',sprintf(' %2.0f ',jits_nan))); ...
            end
            disp(sprintf(' %8.6f ',jit_crits(:,itype)))
        end
        if jit_crit_choice>0
            jit_rms_list=jit_crits(:,jit_crit_choice);
        else
            jit_rms_list=sqrt(sum(jit_crits.^2,2));
        end
        disp('rms jitter used')
        disp(sprintf(' %8.6f ',jit_rms_list));
    else
        disp('critical jitters not calculated, and confidence limits will be skipped');
    end
    %
    disp(sprintf('stimulus coordinates group along %2.0f rays',nrays));
    %
    ray_counts=full(sparse(rayss{iset}.whichray(rayss{iset}.whichray>0),ones(sum(rayss{iset}.whichray>0),1),1,nrays,1));
    for iray=1:nrays
        disp(sprintf('ray %2.0f: %2.0f points; endpoint: %s',iray,ray_counts(iray),sprintf('%5.2f',rayss{iset}.endpt(iray,:))));
    end
    %find stimulus label at end of ray, in positive direction when posssible
    ray_labels{iset}=cell(1,nrays);
    for iray=1:nrays
        mults_ray=rayss{iset}.mult(rayss{iset}.whichray==iray);
        maxend=intersect(find(abs(rayss{iset}.mult)==max(abs(mults_ray))),find(rayss{iset}.whichray==iray));
        maxend=maxend(find(rayss{iset}.mult(maxend)==max(rayss{iset}.mult(maxend)))); %choose positive direction if possible
        ray_labels{iset}{iray}=strrep(strrep(sas{iset}.spec_labels{maxend},' ',''),'=',''); %strip = and space
        zstart=strfind(ray_labels{iset}{iray},zstring); %0.00 to remove?
        if ~isempty(zstart)
            ray_labels{iset}{iray}=ray_labels{iset}{iray}([1:zstart-2 zstart+length(zstring):end]);
        end
        disp(sprintf('ray %2.0f label: %s',iray,ray_labels{iset}{iray})); 
        mult_labels{iset}{iray}=...
            sprintf('%12s%5s    ',ray_labels{iset}{iray},'[-]',ray_labels{iset}{iray},'[+]',ray_labels{iset}{iray},'[-+]');
    end
    ray_pair_labels{iset}=cell(1,nrays*(nrays-1)/2);
    ilab=0;
    for iray=1:nrays-1
        for jray=iray+1:nrays
            ilab=ilab+1;
            ray_pair_labels{iset}{ilab}=cat(2,ray_labels{iset}{iray},':',ray_labels{iset}{jray});
        end
    end
    %
    %for each dimension model, find best-fitting signed and unsigned rays, including the origin
    %
    for idimptr=1:length(dim_list) %compute ray fits and angles
        idim=dim_list(idimptr);
        for if_bid=0:1 %uni- and bi-directional
            opts_stats_use=setfield(opts_stats,'if_bid',if_bid);
            if if_havechoice==0
                opts_stats_use=setfield(opts_stats_use,'nsurrs',0);
            end
            jit_use=jit_rms_list(idim);
            if isnan(jit_use) %will only happen if jitters are all NaN for all dimensions
                jit_use=0;
            end
            [angles{iset,1+if_bid}{idim},mults{iset,1+if_bid}{idim},angles_stats{iset,1+if_bid}{idim},mults_stats{iset,1+if_bid}{idim},rayfit{iset,1+if_bid}{idim}]=...
                psg_raystats(ds{iset}{idim},sas{iset},rayss{iset},jit_rms_list(idim),opts_stats_use);
            %
        end %if_bid
    end %idim_ptr
    %
    %nicely formatted output of mults and angles to console and add to table
    %
    disp(' ');
    disp('gains along each ray, and for bidirectional fit to each axis')
    disp(cat(2,' dim    ',mult_labels{iset}{:}));
    for idimptr=1:length(dim_list) %display ray fits and angles
        idim=dim_list(idimptr);
        stat_names=fieldnames(mults_stats{iset,1}{idim});
        nstats=length(stat_names); %will be zero if nsurrs=0 or choice data file not found
        v=cell(nstats+1,2);
        stats_have=cell(1,nstats+1); %this is to map fields from mults_stats into clo, chi, sem
        for iv=0:nstats
            if (iv==0)
                text_string=sprintf('%3.0f   ',idim);
                for ib=1:2
                    v{iv+1,ib}=mults{iset,ib}{idim};
                end
                stats_have{1}='data';
            else
                stat_name=stat_names{iv};
                stats_have{iv+1}=stat_name;
                text_string=sprintf('%6s',stat_name);
                for ib=1:2
                    v{iv+1,ib}=mults_stats{iset,ib}{idim}.(stat_name);
                end
            end
            for iray=1:nrays
                text_string=cat(2,text_string,sprintf(' %15.4f',v{iv+1,1}.dist_gain(iray,1)),'     '); %gain on negative ray
                text_string=cat(2,text_string,sprintf(' %15.4f',v{iv+1,1}.dist_gain(iray,2)),'     '); %gain on positive ray
                text_string=cat(2,text_string,sprintf(' %15.4f',v{iv+1,2}.dist_gain(iray,1)),'     '); %bidirectional gain
            end
            disp(text_string);
        end %iv 0=data, >=1: stats
        disp(' ');
        %add to table
        mult_vals=NaN(4,3,nrays); %d1: value, clo, chi, sem; d2: unipolar neg, unipolar pos, bidir; d3: each ray
        for iv=0:nstats
            ivptr=strmatch(stats_needed{iv+1},stats_have,'exact'); %look for clo, chi, and sem in output of psg_raystats
            if length(ivptr)==1
                for iray=1:nrays
                    mult_vals(iv+1,1,iray)=v{ivptr,1}.dist_gain(iray,1);
                    mult_vals(iv+1,2,iray)=v{ivptr,1}.dist_gain(iray,2);
                    mult_vals(iv+1,3,iray)=v{ivptr,2}.dist_gain(iray,1);
                end
            end
        end %iv
        vname='dist_gain';
        for inpb=1:3 %1: unipolar negative, 2: unipolar positive, 3: bidirectional
            for iray=1:nrays
                if_bid=ismember(inpb,[3]);
                neg_pos_bid=neg_pos_bid_text{inpb};
                values=mult_vals(:,inpb,iray)';
                data_cell=[{idim,jit_rms_list(idim),if_bid,neg_pos_bid,iray,0,ray_labels{iset}{iray},'',vname} num2cell(values)];
                t_data=array2table(data_cell);
                t_data.Properties.VariableNames=data_variable_names;
                if (if_t_all==0)
                    t_all=[t_meta_set{iset},t_data];
                    if_t_all=1;
                else
                    t_all=[t_all;[t_meta_set{iset},t_data]];
                end
            end %iray
        end %inpb 
    end %idim_ptr
    disp(' ');
    disp('cosines of angles between pos and neg rays, and between bidirectional fits to axes')
    disp(cat(2,' dim    ',sprintf('%20s ',ray_labels{iset}{:}),sprintf('%30s ',ray_pair_labels{iset}{:}))); %header
    for idimptr=1:length(dim_list) %display ray fits and angles
        idim=dim_list(idimptr);
        stat_names=fieldnames(angles_stats{iset,1}{idim});
        nstats=length(stat_names);
        v=cell(nstats+1,2);
        for iv=0:nstats
            if (iv==0)
                text_string=sprintf('%3.0f   ',idim);
                for ib=1:2
                    v{iv+1,ib}=angles{iset,ib}{idim};
                end
                stats_have{1}='data';
            else
                stat_name=stat_names{iv};
                stats_have{iv+1}=stat_name;
                text_string=sprintf('%6s',stat_name);
                for ib=1:2
                    v{iv+1,ib}=angles_stats{iset,ib}{idim}.(stat_name);
                end
            end
            for iray=1:nrays
                text_string=cat(2,text_string,sprintf(' %15.4f',v{iv+1,1}.cosangs(iray,iray,1,2)),'     '); %cosine of angle between pos and neg direction on each axis
            end
            for iray=1:nrays-1
                for jray=iray+1:nrays
                    text_string=cat(2,text_string,sprintf(' %25.4f',v{iv+1,2}.cosangs(iray,jray)),'     '); %cosine of angle between bidirectional fits of two axes
                end
            end
            disp(text_string);
        end %iv 0=data, >=1: stats
        disp(' ');
        %add to table
        nraypairs=nrays*(nrays-1)/2; %angles between pos and neg, and then all pairs
        ang_self=NaN(4,nrays); %d1: value, clo, chi, sem; d2: angle between neg and pos ray at origin
        ang_pair=NaN(4,nraypairs);  %d1: value, clo, chi, sem; d2: each ray pair
        for iv=0:nstats
            ivptr=strmatch(stats_needed{iv+1},stats_have,'exact'); %look for clo, chi, and sem in output of psg_raystats
            if length(ivptr)==1
                for iray=1:nrays
                    ang_self(iv+1,iray)=v{ivptr,1}.cosangs(iray,iray,1,2);
                end
                ijray=0;
                for iray=1:nrays-1
                    for jray=iray+1:nrays
                        ijray=ijray+1;
                        ang_pair(iv+1,ijray)=v{ivptr,2}.cosangs(iray,jray);
                    end %jray
                end %iray
            end %found a match
        end %iv
        vname='cosang_self';
        if_bid=0;
        neg_pos_bid='';
        for iray=1:nrays
            if_bid=0;
            values=ang_self(:,iray)';
            data_cell=[{idim,jit_rms_list(idim),if_bid,neg_pos_bid,iray,0,ray_labels{iset}{iray},'',vname} num2cell(values)];
            t_data=array2table(data_cell);
            t_data.Properties.VariableNames=data_variable_names;
            t_all=[t_all;[t_meta_set{iset},t_data]];
        end %iray
        vname='cosang_pair';
        if_bid=1;
        neg_pos_bid=neg_pos_bid_text{3};
        ijray=0;
        for iray=1:nrays-1
            for jray=iray+1:nrays
                ijray=ijray+1;
                values=ang_pair(:,ijray)';
                data_cell=[{idim,jit_rms_list(idim),if_bid,neg_pos_bid,iray,jray,ray_labels{iset}{iray},ray_labels{iset}{jray},vname} num2cell(values)];
                t_data=array2table(data_cell);
                t_data.Properties.VariableNames=data_variable_names;
                t_all=[t_all;[t_meta_set{iset},t_data]];
            end %jray
        end %iray
    end %idim_ptr
end %iset
%
disp('adding overall settings to UserData of tables')
settings=struct;
settings.opts_stats=opts_stats;
settings.opts_lljit=opts_lljit;
settings.pval=pval;
settings.jit_crit_choice=jit_crit_choice;
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


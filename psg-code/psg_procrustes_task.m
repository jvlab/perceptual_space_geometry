%psg_procrustes_task: special-purpose script to compare perceptual spaces
%via Procrustes across subjects and tasks
%
%  See also: PSG_RAYSTATS_DBPLOT, PSG_GET_COORDSETS, PROCRUSTES.
%
if ~exist('setup_path') setup_path='../psg/psg_data/'; end
if ~exist('data_path') psg_path='../psg/psg_data/'; end
if ~exist('qform_path') qform_path='../stim/'; end
%
if ~exist('task_list') %tasks to analyze, using psg_raystats_dbplot
    task_list={'threshold','similarity','brightness','working_memory','unconstrained_grouping'}; %a subset of expt_grp
end
if ~exist('subj_list') %corresponds to subj_model_ID in psg_raystats_dbplot databases
    subj_list={'bl','mc','nf','sn','zk'}; %could also add cme, saw, but more incomplet
end
if ~exist('stim_list') %corresponds to expt_uid in psg_raystats_dbplot databases
    stim_list={'bgca3pt','dgea3pt'}; %could also add bc63pt, bcpm3pt, etc
end
if ~exist('opts_read')
    opts_read=struct;
end
opts_read.if_auto=1;
opts_read=filldefault(opts_read,'if_log',0);
%
subjs_have_personal={'df','dt','jd','mc'}; %subjects with personalizd models
subjs_have_custavg={'bl','nf','saw','sn','zk'}; %subjects with customized averge models
qform_pref='btc_allraysfixedb_';
qform_suff='_100surrs_madj.mat';
% 
%from raystats_dbplot
%
% values_short.brightness='bright';
% values_short.constrained_grouping='paired';
% values_short.unconstrained_grouping='group';
% values_short.similarity='standard';
% values_short.dis_similarity='dissim';
% values_short.working_memory='workmem';
% values_short.threshold='thresh';
%
%thi section can be used generically to read in all available data
expt_grp={'threshold','similarity','dis_similarity','working_memory','constrained_grouping','unconstrained_grouping','brightness'};
expt_name={'','','dis','wm','gp','gm','br'};
wm_dur_string='1000';
setup_suffix='9';
%
ntasks=length(task_list);
nsubjs=length(subj_list);
nstims=length(stim_list);
task_infix=cell(1,ntasks);
for itask=1:ntasks
    imatch=strmatch(task_list{itask},expt_grp,'exact');
    task_infix{itask}=expt_name{imatch};
    if ~isempty(expt_name{imatch})
        task_infix{itask}=cat(2,'-',task_infix{itask});
    end
    if strmatch(expt_name{imatch},'wm','exact')
        task_infix{itask}=cat(2,task_infix{itask},wm_dur_string);
    end
end
%read all datassets
ds=cell(ntasks,nsubjs,nstims);
sets=cell(ntasks,nsubjs,nstims);
sas=cell(ntasks,nsubjs,nstims);
opts_read_used=cell(ntasks,nsubjs,nstims);
dims_avail=NaN(ntasks,nsubjs,nstims);
for itask=1:ntasks
    switch task_list{itask}
        case 'threshold'
            input_type=2;
        otherwise
            input_type=1;
    end
    for isubj=1:nsubjs
        for istim=1:nstims
            setup_file=cat(2,stim_list{istim},setup_suffix,'.mat');
            if input_type==1
                data_path=psg_path;
                data_file=cat(2,stim_list{istim},'_coords_',upper(subj_list{isubj}),task_infix{itask},'_sess01_10.mat');
            end
            if input_type==2
                data_path=qform_path;
                if strmatch(subj_list{isubj},subjs_have_personal,'exact') 
                    data_file=cat(2,qform_pref,subj_list{isubj},qform_suff);
                elseif strmatch(subj_list{isubj},subjs_have_custavg,'exact')
                    data_file=cat(2,qform_pref,'avg-',subj_list{isubj},qform_suff);
                else
                    data_file=cat(2,qform_pref,'avg,',qform_suff);
                end
            end
            disp(sprintf('task %30s       subj %4s       stim %5s:   %20s, %s',task_list{itask},subj_list{isubj},stim_list{istim},setup_file,data_file))
            if_missing=0;
            setup_full=cat(2,setup_path,filesep,setup_file);
            setup_full=strrep(strrep(setup_full,'\','/'),'//','/');
            data_full=cat(2,data_path,filesep,data_file);
            data_full=strrep(strrep(data_full,'\','/'),'//','/');
            if ~exist(setup_full,'file')
                disp(sprintf('setup file %s not found.',setup_full));
                if_missing=1;
            end
            if ~exist(data_full,'file')
                disp(sprintf('data file %s not found.',data_full));
                if_missing=1;
            end
            %now read the coords, and keep track of max dimension (0 for no data)
            if (if_missing==0)
                opts_read_use=opts_read;
                opts_read_use.input_type=input_type;
                opts_read_use.data_fullnames{1}=data_full;
                opts_read_use.setup_fullnames{1}=setup_full;
                [set,d,sa,rayss,oru]=psg_get_coordsets(opts_read_use,struct(),struct(),1);
                sets{itask,isubj,istim}=set{1};
                ds{itask,isubj,istim}=d{1};
                sas{itask,isubj,istim}=sa{1};
                dims_avail(itask,isubj,istim)=length(d{1});
                opts_read_used{itask,isubj,istim}=oru;
            end
        end %stim
    end %subj
end %task
nsets_avail=sum(~isnan(dims_avail(:)));
max_dims_avail=min(dims_avail(:),[],'omitnan');
disp(' ');
disp(sprintf('%4.0f datasets found, out of %4.0f (%4.0f tasks x %4.0f subjects x %4.0f stimulus sets)',nsets_avail,...
    ntasks*nsubjs*nstims,ntasks,nsubjs,nstims));
disp(sprintf('maximum dimension available across all datasets: %3.0f',max_dims_avail));
%dss and full list of subjects, paradigms, and coord sets, along with stim labels, are useful for general exprting
%may want to check that all dimensions are present

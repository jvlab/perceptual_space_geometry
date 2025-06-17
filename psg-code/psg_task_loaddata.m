function [sets,ds,sas,opts_read_used,paths_used,dlists_used]=psg_task_loaddata(dlists,paths,opts)
% [sets,ds,sas,opts_read_used,paths_used,dlists_used]=psg_task_loaddata(dlists,paths,opts)
% is a utiltiy to load perceptual space coordinates for the task study
%
% dlists is a structure with entries task_list,subj_list,stim_list
% paths is a structure with data paths
% opts is a structure of options for reading files
%
%[sets,ds,sas]{itask,isubj,istim} contains the corresponding structures returned by psg_get_coordsets
%
%  See also: PSG_READ_COORDDATA, PSG_PROCRUSTES_TASK, FILLDEFAULT.
%
dlists=filldefault(dlists,'task_list',{'threshold','similarity','brightness','working_memory','unconstrained_grouping'}); %a subset of expt_grp
dlists=filldefault(dlists,'subj_list',{'bl','mc','nf','sn','zk'}); %could also add cme, saw, but these subjs are more incomplete
dlists=filldefault(dlists,'stim_list',{'bgca3pt','dgea3pt'}); %could also add bc63pt, bcpm3pt, etc
dlists_used=dlists;
%
paths=filldefault(paths,'setup_path','../psg/psg_data/');
paths=filldefault(paths,'psg_path','../psg/psg_data/');
paths=filldefault(paths,'qform_path','../stim/');
paths_used=paths;
%
opts_read=struct;
opts_read.if_auto=1;
opts=filldefault(opts,'opts_read',opts_read);
opts_read=opts.opts_read;
opts_read=filldefault(opts_read,'if_log',0);
%
task_list=dlists.task_list;
subj_list=dlists.subj_list;
stim_list=dlists.stim_list;
%
setup_path=paths.setup_path;
psg_path=paths.psg_path;
qform_path=paths.qform_path;
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
%read all datasets
ds=cell(ntasks,nsubjs,nstims);
sets=cell(ntasks,nsubjs,nstims);
sas=cell(ntasks,nsubjs,nstims);
opts_read_used=cell(ntasks,nsubjs,nstims);
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
                    data_file=cat(2,qform_pref,'avg',qform_suff);
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
                opts_read_used{itask,isubj,istim}=oru;
            end
        end %stim
    end %subj
end %task
return

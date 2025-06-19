function [sets,ds,sas,opts_read_used,paths_used,dlists_used,choices,expt_labels]=psg_task_loaddata(dlists,paths,opts)
% [sets,ds,sas,opts_read_used,paths_used,dlists_used,choices,expt_labels]=psg_task_loaddata(dlists,paths,opts)
% is a utiltiy to load perceptual space coordinates for the task study
%
% dlists is a structure with entries task_list,subj_list,stimset_list
% paths is a structure with data paths
% opts is a structure of options for reading files
%   opts.if_choices=1 to attempt to read choice data, defaults to 0
%
% [sets,ds,sas]{itask,isubj,istimset} contains the corresponding structures returned by psg_get_coordsets 
% choices{itask,isubj,istimset} contains the corresponding structures returned by psg_read_choicedata
% expt_labels conains expt_grp, expt_name, other details
%
%  See also: PSG_READ_COORDDATA, PSG_PROCRUSTES_TASK, FILLDEFAULT, PSG_READ_CHOICEDATA, PSG_DIMSTAT_TASK.
%
dlists=filldefault(dlists,'task_list',{'threshold','similarity','brightness','working_memory','unconstrained_grouping'}); %a subset of expt_grp
dlists=filldefault(dlists,'subj_list',{'bl','mc','nf','sn','zk'}); %could also add cme, saw, but these subjs are more incomplete
dlists=filldefault(dlists,'stimset_list',{'bgca3pt','dgea3pt'}); %could also add bc63pt, bcpm3pt, etc
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
opts=filldefault(opts,'if_choices',0);
%
task_list=dlists.task_list;
subj_list=dlists.subj_list;
stimset_list=dlists.stimset_list;
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
expt_labels.expt_grp=expt_grp;
expt_labels.expt_name=expt_name;
expt_labels.wm_dur_string=wm_dur_string;
expt_labels.setup_suffix=setup_suffix;
%
ntasks=length(task_list);
nsubjs=length(subj_list);
nstimsets=length(stimset_list);
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
ds=cell(ntasks,nsubjs,nstimsets);
sets=cell(ntasks,nsubjs,nstimsets);
sas=cell(ntasks,nsubjs,nstimsets);
choices=cell(ntasks,nsubjs,nstimsets);
opts_read_used=cell(ntasks,nsubjs,nstimsets);
for itask=1:ntasks
    switch task_list{itask}
        case 'threshold'
            input_type=2;
        otherwise
            input_type=1;
    end
    for isubj=1:nsubjs
        for istimset=1:nstimsets
            setup_file=cat(2,stimset_list{istimset},setup_suffix,'.mat');
            if input_type==1
                data_path=psg_path;
                data_file=cat(2,stimset_list{istimset},'_coords_',upper(subj_list{isubj}),task_infix{itask},'_sess01_10.mat');
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
            disp(sprintf('task %30s       subj %4s       stim %5s:   %20s, %s',task_list{itask},subj_list{isubj},stimset_list{istimset},setup_file,data_file))
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
                sets{itask,isubj,istimset}=set{1};
                ds{itask,isubj,istimset}=d{1};
                sas{itask,isubj,istimset}=sa{1};
                opts_read_used{itask,isubj,istimset}=oru;
                if opts.if_choices & input_type==1
                    choices_full=strrep(data_full,'coords','choices');
                    if ~exist(choices_full,'file')
                       disp(sprintf('choice file %s not found.',choices_full));
                    else
                        choices{itask,isubj,istimset}=psg_read_choicedata(choices_full,setup_full,opts_read_use);
                    end
                end %if_choices
            end %data file present
        end %stim
    end %subj
end %task
return

%psg_procrustes_task: special-purpose script to compare perceptual spaces
%via Procrustes across subjects and tasks
%
%  See also: PSG_GET_COORDSETS, PROCRUSTES, PSG_TASK_LOADDATA.
%
if ~exist('dlists') dlists=struct; end
dlists=filldefault(dlists,'task_list',{'threshold','similarity','brightness','working_memory','unconstrained_grouping'}); %a subset of expt_grp
dlists=filldefault(dlists,'subj_list',{'bl','mc','nf','sn','zk'}); %could also add cme, saw, but these subjs are more incomplete
dlists=filldefault(dlists,'stim_list',{'bgca3pt','dgea3pt'}); %could also add bc63pt, bcpm3pt, etc
%
if ~exist('paths') paths=struct; end
paths=filldefault(paths,'setup_path','../psg/psg_data/');
paths=filldefault(paths,'psg_path','../psg/psg_data/');
paths=filldefault(paths,'qform_path','../stim/');
%
if ~exist('opts_read') opts_read=struct; end
opts_read.if_auto=1;
opts_read=filldefault(opts_read,'if_log',0);
if ~exist('opts') opts=struct; end
opts=filldefault(opts,'opts_read',opts_read);
%
subjs_have_personal={'df','dt','jd','mc'}; %subjects with personalizd models
subjs_have_custavg={'bl','nf','saw','sn','zk'}; %subjects with customized averge models
qform_pref='btc_allraysfixedb_';
qform_suff='_100surrs_madj.mat';
%
[sets,ds,sas,opts_read_used,paths_used,dlists_used]=psg_task_loaddata(dlists,paths,opts);
ntasks=length(dlists.task_list);
nsubjs=length(dlists.subj_list);
nstims=length(dlists.stim_list);
%
dims_avail=NaN(ntasks,nsubjs,nstims);
for itask=1:ntasks
    for isubj=1:nsubjs
        for istim=1:nstims
            if ~isempty(ds{itask,isubj,istim})
                dims_avail(itask,isubj,istim)=length(ds{itask,isubj,istim});
            end
        end %stim
    end %subj
end %task
%
nsets_avail=sum(~isnan(dims_avail(:)));
max_dims_avail=min(dims_avail(:),[],'omitnan');
disp(' ');
disp(sprintf('%4.0f datasets found, out of %4.0f (%4.0f tasks x %4.0f subjects x %4.0f stimulus sets)',nsets_avail,...
    ntasks*nsubjs*nstims,ntasks,nsubjs,nstims));
disp(sprintf('maximum dimension available across all datasets: %3.0f',max_dims_avail));
%

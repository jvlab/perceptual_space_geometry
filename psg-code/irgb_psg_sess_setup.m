%irgb_psg_sess_setup: sets up perceptual space geometry session files
% for a generic stimulus set, designated irgb for convenience but stimuli need not be made
% with irgb_spec_make.  
% This only makes the session files (csv), and saves a mat-file with the
% information needed to create the stimulus files.
%
%  Derived from irgb_psg_setup.
%
% stimulus file names: like mater06_filt_bw_whiten023_216.png,
% mater06 is the paradigm ID, filt_bw_whiten is the manip, 023 designates
% the image id number (1 to nstims), 216 is the instance
% 
% spec_params, if specified, determines the stimulus set, otherwise defaults are taken from irgb_spec_make
%
% See also:  IRGB_PSG_SETUP, PSG_DEFOPTS, PSG_COND_CREATE, PSG_COND_WRITE,
% PSG_SESSCONFIG_MAKE, IRGB_MODIFY, IRGB_PSG_IMGS_SETUP, IRGB_MODIFY_DICT.
%
nrgb=3;
%
%allowed stimulus manipulations
manip_list=fieldnames(irgb_modify_dict());
%
%psg defaults
%
if ~exist('opts_psg') opts_psg=struct; end
opts_psg=psg_defopts(opts_psg);
if ~exist('spec_params') spec_params=struct; end
if ~exist('opts_stim') opts_stim=struct; end
if ~exist('image_pixels') image_pixels=144; end %subsamples for stimuli
%
%defaults for irgb
if ~exist('cond_file_prefix') cond_file_prefix='irgb';end
opts_psg.cond_nstims=37; %default for generic experiment
opts_psg.cond_nstims=getinp('number of stimuli','d',[13 49],opts_psg.cond_nstims);
nstims=opts_psg.cond_nstims;
%
if_frozen_psg=getinp('1 for frozen random numbers, 0 for new random numbers each time for session configuration, <0 for a specific seed','d',[-10000 1]);
%
if (if_frozen_psg~=0)
    rng('default');
    if (if_frozen_psg<0)
        rand(1,abs(if_frozen_psg));
    end
else
    rng('shuffle');
end
%
ifok=0;
while (ifok==0)
    disp(sprintf('current psg spoke_setup: %3.0f augmented stimuli, %3.0f to be replaced, %3.0f comparison stimuli per trial, overlap %3.0f; %3.0f sessions',...
        opts_psg.cond_nstims, opts_psg.cond_nstims_toreplace,opts_psg.cond_ncompares,opts_psg.cond_novlp,opts_psg.cond_nsess));
    if (nstims==opts_psg.cond_nstims-opts_psg.cond_nstims_toreplace)
        if getinp('1 if ok','d',[0 1])
            ifok=1;
        end
    else
        disp(sprintf('current setup invalid: number of stimuli in spokes: %3.0f, in psg: %3.0f',...
            nstims,opts_psg.cond_nstims));
    end
    if (ifok==0)
        opts_psg.cond_nstims=getinp('nstims in augmented setup (prior to replacement)','d',[1 1000],opts_psg.cond_nstims);
        opts_psg.cond_nstims_toreplace=getinp('nstims to replace (for session generation)','d',[0 1000],opts_psg.cond_nstims-nstims);
        opts_psg.cond_ncompares=getinp('ncompares','d',[1 1000],opts_psg.cond_ncompares);
        opts_psg.cond_novlp=getinp('novlp','d',[1 1000],opts_psg.cond_novlp);
        opts_psg.cond_nsess=getinp('number of sessions','d',[1 100],opts_psg.cond_nsess);
    end
end
for k=1:length(opts_psg.refseq_labels)
    disp(sprintf('%1.0f->method for choosing stimuli in overlap: %s',k,opts_psg.refseq_labels{k}));
end
opts_psg.refseq=getinp('choice','d',[1 length(opts_psg.refseq_labels)],opts_psg.refseq);
%
%modified code from faces_mpi_psg_setup to deal with stimulus replacement
%
[sessions_withaug,sessions_withaug_sorted,opts_used,psg_desc]=psg_sessconfig_make(opts_psg);
opts_psg.cond_desc=psg_desc;
s.if_frozen_psg=if_frozen_psg;
%
%accumulate and display statistics of the augmented configuration
%
disp(sprintf('Analyzing the session configuration %s prior to replacement',psg_desc));
session_stats_withaug=psg_session_stats(sessions_withaug,setfield(opts_psg,'if_log',1));
s.sessions_withaug=sessions_withaug;
s.session_stats_withaug=session_stats_withaug;
%
%replace the extra stimuli and do stats
[sessions,opts_psg,warnings]=psg_sessconfig_replace(sessions_withaug,opts_psg);
if ~isempty(warnings)
    disp(warnings);
end
disp(sprintf('Analyzing the session configuration %s after replacement',psg_desc));
%
disp(sprintf('Analyzing spoke_setup with %s',opts_psg.cond_desc));
ntrials=size(sessions,1);
%
%accumulate and display statistics of the configuration
%
session_stats=psg_session_stats(sessions,setfield(opts_psg,'if_log',1));
%
%determine options for stimulus re-use and create stimulus file names with example numbers
%for each session
%
disp('options for stimulus example re-use')
for k=1:length(opts_psg.example_infix_labels)
    disp(sprintf('%1.0f->%s',k,opts_psg.example_infix_labels{k}));
end
opts_psg.example_infix_mode=getinp('mode','d',[1 length(opts_psg.example_infix_labels)],opts_psg.example_infix_mode);
s.paradigm_name=getinp('paradigm name (indicates group of original image files), e.g., mater01','s',[]);
%
for k=1:length(manip_list)
    disp(sprintf('%2.0f->%s',k,manip_list{k}));
end
manip=getinp('image manipulation choice','d',[1 length(manip_list)]);
image_pixels=getinp('number of pixels in final images','d',[1 1024],image_pixels);
%
s.typenames=cell(nstims,1);
for istim=1:nstims
    s.typenames{istim}=cat(2,manip_list{manip},zpad(istim,3));
end
filename_prefix=cat(2,s.paradigm_name,'_'); 
%filename prefix needed since file names from typenames are too generic , e.g. 1cov1_meandir2_meanmult1)
[session_cells,perms_used,examps_used]=psg_cond_create(sessions,s.typenames,setfield(opts_psg,'prefix',filename_prefix));
%
%add fields to s (in contrast to psg_spokes_setup, many fields including nstims already set above)
%
s.opts_psg=opts_psg;
s.session_stats=session_stats;
s.sessions=sessions;
s.session_cells=session_cells;
s.perms_used=perms_used;
s.examps_used=examps_used;
s.if_frozen_psg=if_frozen_psg;
s.image_pixels=image_pixels;
s.image_manipulation_name=manip_list{manip};
s.image_manipulation_params=struct(); %for future
%
disp('key variables')
disp(s)
%
%write out the condition file, stimulus files, etc.
%
ifok=0;
while (ifok==0)
    pathname=getinp('relative path for condition file and stimulus file','s',[],'./');
    %convert to system-specific separator, append to end, and remove duplicates
    pathname=strrep(cat(2,strrep(pathname,'/',filesep),filesep),cat(2,filesep,filesep),filesep);
    if ~exist(cat(2,'./',pathname),'dir')
        disp(sprintf(' path %s does not exist.',pathname))
        if_create=getinp('1 to create','d',[0 1]);
        if (if_create==1)
            [status,result]=dos(sprintf('mkdir %s',pathname));
            if (status==0)
                ifok=1;
            else
                disp(result)
            end
        end
    else
        disp(sprintf('path %s exists',pathname));
        ifok=1;
    end
    if (ifok==1)
        ifok=getinp('1 if ok','d',[0 1]);
    end
end
%
s.creation_time=datestr(now);
filename_base=getinp('file name base (and path), _sess[#].csv will be appended for cond files','s',[],s.paradigm_name);
filename_mat=cat(2,pathname,filename_base,'.mat');
save(filename_mat,'s');
disp(sprintf('key variables saved in %s',filename_mat));
%
for isess=1:opts_psg.cond_nsess
    filename=cat(2,pathname,filename_base,'_sess',zpad(isess,opts_psg.sess_zpad));
    psg_cond_write(filename,session_cells{isess},setfield(opts_psg,'if_log',1));
end
disp(sprintf('max number of stimulus examples required is %5.0f',1+max(s.examps_used(:)))); %s.examps_used starts at 0

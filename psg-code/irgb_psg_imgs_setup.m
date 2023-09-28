%irgb_psg_imgs_setup: sets up perceptual space geometry image files
% for a generic stimulus set, to accompany session files made by irgb_psg_sess_setup.
%
% stimulus file names: like mater06_filt_bw_whiten023_216.png,
% mater06 is the paradigm ID, filt_bw_whiten is the manip, 023 designates
% the image id number (1 to nstims), 216 is the instance
% 
% See also:  PSG_DEFOPTS, IRGB_PSG_SESS_SETUP, PSG_SESSCONFIG_MAKE, IRGB_MODIFY.
%
%%%need to calculate local image stats
%%%which stimuli can be flipped?  different phase randomizations? -- in irgb_modify
%%%save the file list and the stimuli into an augmented setup file
%
if ~exist('opts_psg') opts_psg=struct; end
if ~exist('setup_mat_file_def') setup_mat_file_def='./mater/mater06.mat'; end
if ~exist('orig_img_file_list_def') orig_img_file_list_def='irgb_psg_imgs_file_list_test.mat'; end
if ~exist('orig_img_file_path_def') orig_img_file_path_def='./GieselZaidiImages/'; end
opts_psg=psg_defopts(opts_psg);
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
setup_mat_file=getinp('setup mat file (and path)','s',[],setup_mat_file_def);
s=getfield(load(setup_mat_file),'s');
disp(s);
nexamps=1+max(s.examps_used(:));
nstims=size(s.typenames,1);
nsess=length(s.session_cells);
opts_psg.cond_nstims=nstims;
opts_psg.cond_nsess=nsess;
%
disp(sprintf('number of  stimuli: %5.0f',nstims));
disp(sprintf('number of sessions: %5.0f',nsess));
disp(sprintf('max number of stimulus examples required is %5.0f',nexamps)); %s.examps_used starts at 0
%
orig_img_file=cell(nstims,1); %the file names of images selected for stimuli
imgs_orig=cell(nstims,1); %the images
if_ok=0;
while (if_ok==0)
    orig_img_file_list=getinp('original image file list','s',[],orig_img_file_list_def);
    orig_img_full_list=getfield(load(orig_img_file_list),'file_list');
    orig_img_count=length(orig_img_full_list);
    disp(sprintf('list contains %4.0f files',orig_img_count));
    if (orig_img_count<nstims)
        disp('not enough stimuli.');
    elseif orig_img_count>nstims
        disp('Need to select a subset')
        orig_img_select=getinp(sprintf('list of %4.0f elements',nstims),'d',[1 orig_img_count],[1:nstims]);
        orig_img_select=unique(orig_img_select);
        if length(orig_img_select)==nstims
            if_ok=1;
        else
            disp('Wrong length');
        end
    else
        orig_img_select=[1:nstims];
        if_ok=1;
    end
    if (if_ok==1)
        orig_img_file_path=getinp('original image file path','s',[],orig_img_file_path_def);
        disp('reading original image files ')
        orig_img_file_path_def=orig_img_file_path;
        for istim=1:nstims
            orig_img_file{istim}=cat(2,orig_img_file_path,orig_img_full_list{orig_img_select(istim)});
            if exist(orig_img_file{istim},'file')
                imgs_orig{istim}=imread(orig_img_file{istim});
                disp(sprintf('%40s read, original image for stimulus %2.0f',orig_img_file{istim},istim));
            else
                disp(sprintf('%40s not found',orig_img_file{istim}));
                if_ok=0;
            end
        end
        if (if_ok==1)
            if_ok=getinp('1 if ok','d',[0 1],if_ok);
        end
    end
end
%
%default path for output images is same as path for setup
img_file_start=max(0,max(union(find(setup_mat_file=='/'),find(setup_mat_file=='\'))));
stim_img_file_path_def=setup_mat_file(1:img_file_start); %get path
stim_img_file_path=getinp('stimulus image file path','s',[],stim_img_file_path_def);
%
for istim=1:nstims
    for iexamp=1:nexamps
        file_name_base=cat(2,stim_img_file_path,s.typenames{istim},'*.',opts_psg.stim_filetype);
        %file name: want to change this
        outfile=strrep(file_name_base,'*',cat(2,'_',zpad(iexamp-1,opts_psg.example_infix_zpad)));
    end
    disp(sprintf('%3.0f files (%25s, %s to %s) written, for stimulus %2.0f of %2.0f, from %s.',nexamps,file_name_base,...
        zpad(0,opts_psg.example_infix_zpad),zpad(nexamps-1,opts_psg.example_infix_zpad),istim,nstims,orig_img_file{istim}));
end
%information needed to specify stimuli
s_aug=struct;
s_aug.orig_img_file_list=orig_img_file_list; %list of image files to select from, may be more than nstims
s_aug.orig_img_select=orig_img_select; %pointers to which image files are used
s_aug.orig_img_file=orig_img_file; %list of image files used
s_aug.orig_img_file_path=orig_img_file_path; %path to image files
disp('s_aug');
disp(s_aug);
%above needs to be saved somewhere

% %
% %allowed stimulus manipulations
% manip_list={'orig_bw','bw_whiten','bw_randph','filt_bw','filt_bw_whiten','filt_bw_randph'};
% %

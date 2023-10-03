%irgb_psg_imgs_setup: sets up perceptual space geometry image files
% for a generic stimulus set, to accompany session files made by
% irgb_psg_sess_setup, and calculates image statistics
%
% stimulus file names: like mater06_filt_bw_whiten023_005.png,
% mater06 is the paradigm ID, filt_bw_whiten is the manip, 023 designates
% the image id number (1 to nstims), 005 is the instance
% 
% See also:  PSG_DEFOPTS, IRGB_PSG_SESS_SETUP, PSG_SESSCONFIG_MAKE, IRGB_MODIFY, IRGB_MODIFY_DICT,
%   IRGB_BTCSTATS.
%
if ~exist('setup_mat_file_def') setup_mat_file_def='./mater/mater06.mat'; end
if ~exist('orig_img_file_list_def') orig_img_file_list_def='irgb_ClothChoices37_file_list.mat'; end % was irgb_psg_imgs_file_list_test
if ~exist('orig_img_file_path_def') orig_img_file_path_def='./GieselZaidiImages/'; end
%
%for image stat calculation
btc_dict=btc_define;
btc_n=length(btc_dict.codel); %number of coords, typically 10
if ~exist('downsamples_def') downsamples_def=[1 2 4 8 16 32]; end
if ~exist('opts_detrend') opts_detrend=[]; end
if ~exist('opts_mlis') opts_mlis=[]; end
%
setup_mat_file=getinp('setup mat file (and path)','s',[],setup_mat_file_def);
s=getfield(load(setup_mat_file),'s');
opts_psg=psg_defopts(s.opts_psg); %take from opts_psg as used by irgb_psg_sess_setup
%
if ~exist('opts_modify') opts_modify=struct;end
opts_modify=filldefault(opts_modify,'nrandpase',16);
opts_modify=filldefault(opts_modify,'range',[0 65536]); %for png
%translate manipulation string into the field in irgb_modify
%
manip_dict=irgb_modify_dict;
manip_list=fieldnames(manip_dict);
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
if_write=getinp('1 to write the image files (0 will create but not write)','d',[0 1]);
downsamples=getinp('downsamplings for calculating image statistics (0 to omit)','d',[0 64],downsamples_def);
if all(downsamples==0)
    downsamples=[];
end
%
disp(s);
nexamps=1+max(s.examps_used(:));
nstims=size(s.typenames,1);
nsess=length(s.session_cells);
%
disp(sprintf('number of  stimuli: %5.0f',nstims));
disp(sprintf('number of sessions: %5.0f',nsess));
disp(sprintf('max number of stimulus examples required is %5.0f (to be numbered %5.0f to %5.0f)',nexamps,...
    opts_psg.example_numoffset+[0 nexamps-1])); %s.examps_used starts at 0
%
if ~isfield(manip_dict,s.image_manipulation_name)
    disp(sprintf('image_manipulation %s not recognized',s.image_manipulation_name));
else
    modify_name=manip_dict.(s.image_manipulation_name).modify_name;
    disp(sprintf('image_manipulation: %s, -> %s',s.image_manipulation_name,modify_name));
end
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
                disp(sprintf('%45s read, original image for stimulus %2.0f',orig_img_file{istim},istim));
            else
                disp(sprintf('%45s not found',orig_img_file{istim}));
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
if if_write
    cwstring='created and written';
else
    cwstring='created';
end
%
s_aug=struct; %additional information needed to specify stimuli
s_aug.stack_sel=cell(1,nstims); %a random selection in the stack of manipulated images for each example of this stimulus
s_aug.jit_sel=cell(1,nstims); %jitters on the two coordinates for manipulated images 
s_aug.btcstats_downsamples=downsamples;
s_aug.btcstats=cell(1,nstims);
for istim=1:nstims
    label=strrep(strrep(orig_img_file{istim},orig_img_file_path,''),opts_psg.stim_filetype,'');
    file_name_base=cat(2,stim_img_file_path,s.paradigm_name,'-',s.typenames{istim},'*.',opts_psg.stim_filetype);
    %
    img_modified_stack=getfield(irgb_modify(imgs_orig{istim},setfield(opts_modify,'label',label)),modify_name);
    nstack=size(img_modified_stack,3);   
    njit=1+size(img_modified_stack,2)-s.image_pixels;
    if njit<=0
        disp(sprintf('original image for stimulus  %2.0f (%25s) is too small (edge: %4.0f, needs to be at least %4.0f).  Skipped.',...
            istim,orig_img_file{istim},size(img_modified_stack,1),s.image_pixels));
    else
        %calculate image statistics
        if ~isempty(downsamples)
            s_aug.btcstats{istim}=zeros(length(downsamples),btc_n,nstack);
            for istack=1:nstack
                [s_aug.btcstats{istim}(:,:,istack),opts_mlis_used,opts_detrend_used]=...
                    irgb_btcstats(img_modified_stack(:,:,istack),downsamples,opts_mlis,opts_detrend);
            end
            if (istack==1) & ~isfield(s_aug,'btcstats_opts_mlis')
                s_aug.btcstats_opts_mlis=opts_mlis_used;
                s_aug.btcstats_opts_detrend=opts_detrend;
            end
        end
        s_aug.stack_sel{istim}=max(1,ceil(nstack*rand(nexamps,1)));
        s_aug.jit_sel{istim}=max(1,ceil(njit*rand(nexamps,2)));
        jits=s_aug.jit_sel{istim};
        for iexamp=1:nexamps
            %output file name, underscore separates example number from rest
            outfile=strrep(file_name_base,'*',cat(2,'_',zpad(iexamp-1+opts_psg.example_numoffset,opts_psg.example_infix_zpad)));
            %
            %choose a random location in the stack, and a random offset
            %
            stim_examp=img_modified_stack(jits(iexamp,1)+[0:s.image_pixels-1],jits(iexamp,2)+[0:s.image_pixels-1],s_aug.stack_sel{istim}(iexamp));
            if (if_write)
                imwrite(stim_examp,outfile);
            end
        end
        disp(sprintf('%3.0f files (%25s, %s to %s) %s, for stimulus %2.0f of %2.0f, from %s.',nexamps,strrep(file_name_base,'*','_*'),...
            zpad(opts_psg.example_numoffset,opts_psg.example_infix_zpad),...
            zpad(nexamps-1+opts_psg.example_numoffset,opts_psg.example_infix_zpad),...
            cwstring,istim,nstims,orig_img_file{istim}));
    end
end
s_aug.orig_img_file_list=orig_img_file_list; %list of image files to select from, may be more than nstims
s_aug.orig_img_select=orig_img_select; %pointers to which image files are used
s_aug.orig_img_file=orig_img_file; %list of image files used
s_aug.orig_img_file_path=orig_img_file_path; %path to image files
s_aug.orig_img_file_short=strrep(orig_img_file,orig_img_file_path,'');
disp('s_aug');
disp(s_aug);
s_aug_fields=fieldnames(s_aug);
for k=1:length(s_aug_fields)
    s.(s_aug_fields{k})=s_aug.(s_aug_fields{k});
end
%save original setup with image file lists added
setup_aug_mat_file=cat(2,strrep(setup_mat_file,'.mat',''),'_',strrep(orig_img_file_list,'.mat',''),'.mat');
setup_aug_mat_file=getinp('file and path for augmented setup','s',[],setup_aug_mat_file);
save(setup_aug_mat_file,'s');
disp(sprintf('%s written.',setup_aug_mat_file));
% %
% %allowed stimulus manipulations
% manip_list={'orig_bw','bw_whiten','bw_randph','filt_bw','filt_bw_whiten','filt_bw_randph'};
% %

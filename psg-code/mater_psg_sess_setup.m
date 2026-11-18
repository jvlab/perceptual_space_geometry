%mater_psg_sess_setup: sets up perceptual space geometry session files for material experiments
% 
% Derived from irgb_psg_sess_setup, in turn from irgb_psg_setup.
%
% Key differences with irgb_psg_sess_setup:
%
%  option to exclude immediate repeats of stimulus examples
%  stimulus images are merely renamed from the fabric images (for pilot) or converted from tif to bmp (for definitive).
%  default path for stimulus images are in ...EY7977/Materials/FabricImages/pilot_circular/ for pilot
%  typical names for original stimulus images are [ddd-ex].png, where ddd is an integer, not necessarily consecutive, and ex indicate the example number (typically 1-4)
%  typical names for image files listed in *csv's are circular-13-x_003.png
%    circular is the typename_prefix, defaults to s.paradigm_name
%    13 corresponds is original stimulus image designator, i.e., ddd above, not zero-padded
%    -x is typename_suffix
%    003 corresponds to the stimulus example, i.e., ex above, zero-padded to 3 digits
%  image size defaults to 334 x 334 pixels
%  number of stimuli per experiment defaults to 25 for pilot, 37 for definitive
%  definitive: defaults to splitting each session in halves
%  no augmented stimulus sets, no stimulus sets with replacement
%  option to select stimuli to use from available stimuli
%  s.typename_prefix: for pilot, derived from paradigm name, which defaults to terminal folder of path for original images (path_orig_img)
%   for definitive, they defualt to 'fabric'
%  consistent use of istim_ptr to index stimuli (1-based), and stim_no=stims_select(istim_ptr) for the label in the stimulus file name
%
% image files are written to the same folder as the session conditions (.csv) and setup (.mat) files
% image files are always converted to png
%
% 01Sep26: add option for circular vignetting
% 03Nov26: options for definitive vs pilot
%
% See also:  IRGB_PSG_SETUP, PSG_DEFOPTS, PSG_COND_CREATE, PSG_COND_WRITE,
% PSG_SESSCONFIG_MAKE, IRGB_PSG_SESS_SETUP, PSG_COND_SESS_SPLIT
% MATER_RATINGS_STIMSELECT_RS.
%
nrgb=3;
dash='-';
if ~exist('opts_psg') opts_psg=struct; end
opts_psg=psg_defopts(opts_psg);
%
if_pilot_def=getinp('1 for pilot (orig images are png, 25 stims), 2 for definitive (orig images are tif, 37 stims)','d',[1 2]);
%subsets of stimuli selected via mater_ratings_stimselect_rs, mater_ratings_stimselect_03Nov26.txt:
subsets_builtin=struct();
subsets_builtin.A=[1,2,15,26,27,36,37,38,40,44,49,50,51,55,56,59,66,70,73,78,80,82,83,84,85,88,91,92,96,99,102,103,104,109,123,129,130];
subsets_builtin.B=[1,2,4,15,26,27,30,36,37,38,44,49,50,51,56,59,65,70,73,75,78,80,82,83,84,85,87,88,91,92,102,103,104,109,120,129,130];
%
switch if_pilot_def
    case 1
        pd_string='pilot';
        img_type='png';
        nstims_def=25;
        if ~exist('path_orig_image_def') path_orig_image_def='../Materials/FabricImages/pilot_circular'; end
        nsplit=1;
        stimset=[];
    case 2
        pd_string='definitive';
        img_type='tif';
        nstims_def=37;
        if ~exist('path_orig_image_def') path_orig_image_def='../Materials/FabricImages/fabrics_16Oct26_resized'; end
        if_splitsess=1;
        nsplit=2;
        disp(subsets_builtin);
        stimset=getinp('stimulus set (typically A or B)','s',[]);
end
%subsets of stimuli selected via mater_ratings_stimselect_rs, mater_ratings_stimselect_03Nov26.txt:
subsets_builtin=struct();
subsets_builtin.A=[1,2,15,26,27,36,37,38,40,44,49,50,51,55,56,59,66,70,73,78,80,82,83,84,85,88,91,92,96,99,102,103,104,109,123,129,130];
subsets_builtin.B=[1,2,4,15,26,27,30,36,37,38,44,49,50,51,56,59,65,70,73,75,78,80,82,83,84,85,87,88,91,92,102,103,104,109,120,129,130];
%
nsplit=getinp('number of parts to split each session into','d',[1 12],nsplit);
if ~exist ('image_type_out') img_type_out='png'; end
%
%params specific to materials images
if ~exist('image_pixels') image_pixels=334; end %texture image resolution
if ~exist('if_convert_bw') if_convert_bw=1; end %converts color to b/w if all r=g=b
if ~exist('typename_suffix') typename_suffix='-x'; end %to distinguish from raw file names
if ~exist('exclude_immed_repeats_def') exclude_immed_repeats_def=1; end
if ~exist('if_circvig') if_circvig=0; end
%
status=1;
while status>0
    path_orig_image=getinp('path for original fabric images','s',[0],strrep(path_orig_image_def,'\','/'));
    %
    path_orig_fix=strrep(strrep(path_orig_image,'/',filesep),'\',filesep);
    while path_orig_fix(end)==filesep
        path_orig_fix=path_orig_fix(1:end-1);
    end
    %
    dir_cmd=cat(2,'dir ',path_orig_fix,filesep,'*.',img_type,' /b');
    [status,dir_list]=dos(dir_cmd);
    if ~isempty(status)
        disp(sprintf('directory not found'));
    end
end
%
%parse the directory
%
file_list=cell(0);
while ~isempty(dir_list)
    [token,dir_list]=strtok(dir_list);
    if ~isempty(token)
        file_list{end+1}=token;
    end
end
nfiles_all=length(file_list);
nfiles_ok=0;
pref_suff=zeros(0,2);
for ifile=1:nfiles_all
    if contains(file_list{ifile},cat(2,'.',img_type))
        file_base=strrep(file_list{ifile},cat(2,'.',img_type),'');
        dash_pos=min(find(file_base==dash));
        if ~isempty(dash_pos)
            pref=str2num(file_base(1:dash_pos-1));
            suff=str2num(file_base(dash_pos+1:end));
            if ~isempty(pref) & ~isempty(suff)
                nfiles_ok=nfiles_ok+1;
                pref_suff(nfiles_ok,:)=[pref suff];
            end
        end
    end
end
stims_avail=unique(pref_suff(:,1));
nstims_avail=length(stims_avail);
disp(sprintf(' %3.0f *.%s files identified, %3.0f not identified; %3.0f unique stimulus labels, up to %3.0f examples per stimulus',...
    nfiles_ok,img_type,nfiles_all-nfiles_ok,nstims_avail,max(unique(pref_suff(:,2)))));
%
nstims=getinp('number of stimuli to use','d',[13 min(49,nstims_avail)],nstims_def);
%
stims_avail=stims_avail(:)';
stims_select=[];
while length(stims_select)~=nstims
    stims_select=getinp('stimulus selection (- to exclude, 0 for builtin)','d',[-max(stims_avail) max(stims_avail)],stims_avail(1:nstims));
    subset_builtin=[];
    if any(stims_select<0)
        stims_select=union(stims_select(stims_select>0),setdiff(stims_avail,-stims_select(stims_select<0)));
    elseif any(stims_select==0) & if_pilot_def==2
        disp(subsets_builtin);
        stimset=getinp('field name to use','s',[],stimset);
        stims_select=subsets_builtin.(stimset);
    end
    stims_select=unique(stims_select);
    if length(stims_select)~=nstims
        disp('choice:')
        disp(stims_select);
        disp(sprintf('must choose %3.0f out of %3.0f stimuli',nstims,nstims_avail));
    end
end
%
%vignetting
%
if_circvig=getinp('1 for circular vignetting','d',[0 1],if_circvig);
if if_circvig 
    vigval=getinp('vignetting value','d',[0 255]);
    rmax=(image_pixels-1)/2;
    xy=[-rmax:rmax];
    rsq=xy.^2;
    rsq=repmat(rsq,image_pixels,1)+repmat(rsq',1,image_pixels);
    vigmask=double(rsq>(rmax+0.5)^2);
    vigstring=sprintf(', vignetted, gray val=%3.0f',vigval);
else
    vigstring=[];
end
%
%load the original files
%
image_dbase=struct;
image_dbase.images=cell(nstims,1);
image_dbase.file_base=cell(nstims,1);
image_dbase.im_type=cell(nstims,1); %0: color, 1: b/w, -1: original is color but converted to b/w
image_dbase.examp_list=cell(nstims,1); %list of example numbers
for istim_ptr=1:nstims
    stim_no=stims_select(istim_ptr);
    nexamps_expected=sum(pref_suff(:,1)==stim_no);
    examp_list=pref_suff(pref_suff(:,1)==stim_no,2);
    image_dbase.images{istim_ptr}=cell(0);
    image_dbase.im_type{istim_ptr}=zeros(0);
    nexamps_avail=0;
    image_dbase.examp_list{istim_ptr}=examp_list;
    for iexamp=1:nexamps_expected
        file_base=cat(2,sprintf('%1.0f',stim_no),dash,sprintf('%1.0f',examp_list(iexamp)),'.',img_type);
        if iexamp==1
            file_base_first=file_base;
        end
        if iexamp==nexamps_expected
            file_base_last=file_base;
        end
        file_name_full=cat(2,path_orig_fix,filesep,file_base);
        if exist(file_name_full,'file')
            im_data=imread(cat(2,path_orig_fix,filesep,file_base));
            if size(im_data,1)==image_pixels & size(im_data,2)==image_pixels
                nexamps_avail=nexamps_avail+1;
                if size(im_data,3)==1
                    im_type=1; %b/w
                else
                    if all(all(im_data(:,:,1)==im_data(:,:,2))) & all(all(im_data(:,:,1)==im_data(:,:,3))) & if_convert_bw
                        im_type=-1; %color but convert to b/w
                        im_data=im_data(:,:,1);
                    else
                        im_type=0;
                    end
                end
                if if_circvig
                    for k=1:size(im_data,3)
                        im_slice=im_data(:,:,k);
                        im_slice(vigmask==1)=vigval;
                        im_data(:,:,k)=im_slice;
                    end
                end
                image_dbase.images{istim_ptr}{nexamps_avail,1}=im_data;
                image_dbase.file_base{istim_ptr}{nexamps_avail,1}=file_base;
                image_dbase.im_type{istim_ptr}(nexamps_avail,1)=im_type;
            else
                warning(sprintf('s has the wrong size, expecting %3.0f x %3.0f, found %3.0f x %3.0f',image_pixels,image_pixels,size(im_data,1),size(im_data,2)));
            end
        else
            warning(sprintf('%s does not exist',file_base));
        end
    end
    disp(sprintf(' loaded images for chosen stimulus index %2.0f (stimulus file no: %3.0f), %3.0f examples, files %15s to %15s, [%2.0f b/w, %2.0f color->b/w, %2.0f color]',...
        istim_ptr,stim_no,nexamps_avail,file_base_first,file_base_last,...
        sum(image_dbase.im_type{istim_ptr}== 1),...
        sum(image_dbase.im_type{istim_ptr}==-1),...
        sum(image_dbase.im_type{istim_ptr}== 0)));
    image_dbase.stim_no(istim_ptr,1)=stim_no;
    image_dbase.nexamps_avail(istim_ptr,1)=nexamps_avail;
end
%
%customize for materials experiment
%
switch if_pilot_def
    case 1
        cond_file_prefix_def=path_orig_fix(max(find(cat(2,filesep,path_orig_fix)==filesep)):end);
    case 2
        cond_file_prefix_def=cat(2,'fabric',stimset);
end
cond_file_prefix=getinp('prefix for condition files','s',[],cond_file_prefix_def);
opts_psg.cond_nstims=nstims;
%
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time for session configuration, <0 for a specific seed','d',[-10000 1]);
%
if (if_frozen~=0)
    rng('default');
    if (if_frozen<0)
        rand(1,abs(if_frozen));
    end
else
    rng('shuffle');
end
%
opts_psg.cond_nstims_toreplace=0; %no stimulus replacement
ifok=0;
while (ifok==0)
    disp(sprintf('current psg spoke_setup: %3.0f stimuli, %3.0f comparison stimuli per trial, overlap %3.0f; %3.0f sessions',...
        opts_psg.cond_nstims,opts_psg.cond_ncompares,opts_psg.cond_novlp,opts_psg.cond_nsess));
    if (nstims==opts_psg.cond_nstims-opts_psg.cond_nstims_toreplace)
        if getinp('1 if ok','d',[0 1])
            ifok=1;
        end
    else
        disp(sprintf('current setup invalid: number of stimuli in spokes: %3.0f, in psg: %3.0f',...
            nstims,opts_psg.cond_nstims));
    end
    if (ifok==0)
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
[sessions,sessions_sorted,opts_used,psg_desc]=psg_sessconfig_make(opts_psg);
opts_psg.cond_desc=psg_desc;
s.if_frozen=if_frozen;
%
%accumulate and display statistics
%
disp(sprintf('Analyzing the session configuration %s prior to replacement',psg_desc));
session_stats=psg_session_stats(sessions,setfield(opts_psg,'if_log',1));
if mod(size(sessions,1),nsplit)~=0
    disp(sprintf('cannot split sessions as requested, trials per session (%4.0f) not a multiple of %2.0f; disabling split',size(sessions,1),nsplit));
    nsplit=1;
else
    disp(sprintf('splitting sessions into %2.0f pieces, %4.0f trials each',nsplit,size(sessions,1)/nsplit));
end
%
%determine options for stimulus re-use and create stimulus file names with example numbers for each session
%
disp('options for stimulus example re-use')
for k=1:length(opts_psg.example_infix_labels)
    disp(sprintf('%1.0f->%s',k,opts_psg.example_infix_labels{k}));
end
opts_psg.example_infix_mode=getinp('mode','d',[1 length(opts_psg.example_infix_labels)],opts_psg.example_infix_mode);
opts_psg.exclude_immed_repeats=getinp('1 to exclude immediate repeats','d',[0 1],exclude_immed_repeats_def);
%
switch if_pilot_def
    case 1
        paradigm_name_default=cond_file_prefix;
    case 2
        paradigm_name_default='fabric';
end
s.paradigm_name=getinp('paradigm name','s',[],paradigm_name_default);
s.typename_prefix=getinp('typename prefix','s',[],s.paradigm_name);
%
s.typenames=cell(nstims,1);
for istim_ptr=1:nstims
    s.typenames{istim_ptr}=strrep(cat(2,s.typename_prefix,'-',sprintf('%1.0f',image_dbase.stim_no(istim_ptr)),typename_suffix),'_','-'); %cannot have an underscore in typename, this is used to parse results files
end
filename_prefix=cat(2,s.paradigm_name,'-'); 
%
opts_psg.example_maxnum=image_dbase.nexamps_avail; %number of available stimuli
opts_psg.example_list=image_dbase.examp_list; %list of example numbers
%
if nsplit>1
    sessions_orig=sessions;
    %split along dimension 1, but reshape original d1: trial, 2: stimulus, d3: session
    sessions=reshape(sessions,[size(sessions_orig,1)/nsplit,nsplit,size(sessions_orig,2),size(sessions_orig,3)]);
    sessions=permute(sessions,[1 3 2 4]);
    sessions=reshape(sessions,[size(sessions_orig,1)/nsplit,size(sessions_orig,2),size(sessions_orig,3)*nsplit]);
    opts_psg.cond_nsess=opts_psg.cond_nsess*nsplit;
    opts_psg.cond_desc=cat(2,opts_psg.cond_desc,sprintf(' (each session split into %1.0f)',nsplit));
    session_stats=psg_session_stats(sessions,setfield(opts_psg,'if_log',1)); %recompute session stats
end
[session_cells,perms_used,examps_used,opts_psg_used]=psg_cond_create(sessions,s.typenames,opts_psg);
%
%add fields to s (in contrast to psg_spokes_setup, many fields including nstims already set above)
%
s.opts_psg=opts_psg_used;
s.session_stats=session_stats;
s.sessions=sessions;
s.session_cells=session_cells;
s.perms_used=perms_used;
s.examps_used=examps_used;
s.if_frozen=if_frozen;
% new with mater_psg_sess_setup
s.image_pixels=image_pixels;
s.cond_file_prefix=cond_file_prefix;
s.path_orig_image=path_orig_image;
s.stim_nos=image_dbase.stim_no;
s.nexamps_avail=image_dbase.nexamps_avail;
%
disp('key variables')
disp(s)
%
%write out the condition file, stimulus files, etc.
%
ifok=0;
while (ifok==0)
    pathname=getinp('relative path for cond, setup, and stimulus files','s',[],'./');
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
switch if_pilot_def
    case 1
        filename_base_def=s.paradigm_name;
    case 2
        filename_base_def=cat(2,'fabric',stimset);
end
filename_base=getinp('file name base for cond and setup files, _sess[#].csv will be appended for cond files','s',[],filename_base_def);
filename_mat=cat(2,pathname,filename_base,'.mat');
save(filename_mat,'s');
disp(sprintf('key variables saved in %s',filename_mat));
%
for isess=1:opts_psg_used.cond_nsess
    filename=cat(2,pathname,filename_base,'_sess',zpad(isess,opts_psg.sess_zpad));
    psg_cond_write(filename,session_cells{isess},setfield(opts_psg,'if_log',1));
end
disp(sprintf('max number of stimulus examples used is %5.0f',1+max(s.examps_used(:)))); %s.examps_used starts at 0
%
%write stimulus files from image database
%use params from opts_psg_used to create output file name
%
if getinp('1 if ok to write stimulus files','d',[0 1]);
    for istim_ptr=1:nstims
        stim_no=image_dbase.stim_no(istim_ptr);
        disp(' ');
        disp(sprintf('for stimulus index %3.0f, material sample %3.0f from %s',istim_ptr,stim_no,path_orig_fix));
        for iexamp=1:image_dbase.nexamps_avail(istim_ptr)
            example_number=image_dbase.examp_list{istim_ptr}(iexamp);
            entry=cat(2,opts_psg_used.prefix,s.typenames{istim_ptr},opts_psg_used.example_infix_string,zpad(example_number,opts_psg_used.example_infix_zpad));
            stimfile_base=cat(2,entry,'.',img_type_out); %convert to png
            disp(sprintf(' example %3.0f: writing image data from %20s into %s%s',iexamp,image_dbase.file_base{istim_ptr}{iexamp},stimfile_base,vigstring));
            stimfile=cat(2,pathname,stimfile_base);
            %write to stimfile
            imwrite(image_dbase.images{istim_ptr}{iexamp},stimfile);
        end
    end
end

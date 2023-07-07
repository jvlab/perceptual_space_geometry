%faces_mpi_psg_setup: sets up perceptual space geometry with faces in MPI face database
%
%  Derived from psg_spokes_setup, but customized for mpi faces
%
%  * images downsampled to 144 pixels wide by block averaging
%  * files created for color and grayscale
%  * invokes "augmented" stimuli since these stimuli typically have 24 or 36, not 25 and 37
%
%  * creates several session types at a time, so that same faces are not reused
%     to exclude some id numbers from face library, put this list in
%     opts_faces.id_num_exclude, which defaults to empty
%
% See also:  PSG_SPOKES_SETUP, FACES_MPI_SETUP_CREATE, FACES_MPI_INVENTORY, FACES_MPI_GET_SETUPS,
% PSG_DEFOPTS, PSG_COND_CREATE, PSG_COND_WRITE, PSG_SESSCONFIG_MAKE, PSG_SESSION_STATS,
% IRGB_PSG_SETUP.
%

%psg defaults
%
if ~exist('opts_psg') opts_psg=struct; end
opts_psg=psg_defopts(opts_psg);
%
%defaults for faces
if ~exist('cond_file_prefix') cond_file_prefix='faces_mpi_';end
gyfc_infix={'gy','fc'}; %gray or full color
gyfc_labels={'grayscale','fullcolor'};
if ~exist('opts_faces_inventory') opts_faces_inventory=struct; end
if ~exist('faces_mpi_inventory_filename') faces_mpi_inventory_filename='faces_mpi_inventory.mat'; end
if ~exist('faces_mpi_database_path') faces_mpi_database_path='..\faces\mpi_faces\FacesSetASetB\'; end
%
if ~exist('opts_faces') opts_faces=struct; end
opts_faces=filldefault(opts_faces,'image_width',144);
%
load(faces_mpi_inventory_filename);
disp(sprintf('loaded %s',faces_mpi_inventory_filename));
disp(sprintf('total entries in database from table: %5.0f  from array: %5.0f',...
    size(faces_mpi_table,1),sum(faces_mpi_array(:))));
disp('attributes:')
disp(fieldnames(faces_mpi_attrib_info));
%
ifok=0;
while ifok==0
    [faces_mpi_setups,faces_mpi_info,opts_faces_used]=faces_mpi_setup_create(faces_mpi_array,faces_mpi_attrib_info,faces_mpi_table,opts_faces);
    %
    %display images and create session files
    %
    imgs=cell(1,length(faces_mpi_setups)); %this will hold downsampled images
    s_all=cell(1,length(faces_mpi_setups));
    for isetup_ptr=1:length(faces_mpi_info.setup_choices)
        isetup=faces_mpi_info.setup_choices(isetup_ptr);
        fms=faces_mpi_setups{isetup};
        disp(sprintf('working on orig setup %3.0f (%12s) %50s',isetup,fms.name_brief,fms.name));
        tstring=sprintf(' setup %2.0f->%s',isetup,fms.name_brief);
        imgs{isetup}=cell(fms.nstims,1);
        %
        s_all{isetup}=struct;
        s_all{isetup}.faces_mpi_setup=fms;
        s_all{isetup}.faces_mpi_info=faces_mpi_info;
        s_all{isetup}.opts_faces_used=opts_faces_used;
        %
        s_all{isetup}.nstims=fms.nstims;
        s_all{isetup}.typenames=fms.typenames;
        s_all{isetup}.nchecks=opts_faces.image_width;
        s_all{isetup}.nsubsamp=1;
        s_all{isetup}.specs=fms.specs;
        s_all{isetup}.spec_labels=fms.spec_labels;
        %
        figure;
        set(gcf,'Position',[50 100 1200 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',tstring);
        for istim=1:fms.nstims
            switch fms.display.dim_first
                case 'column' %column first
                    ir=1+mod(istim-1,fms.display.rc(1));
                    ic=1+floor((istim-1)/fms.display.rc(1));
                case 'row' %row first
                    ic=1+mod(istim-1,fms.display.rc(2));
                    ir=1+floor((istim-1)/fms.display.rc(2));
            end
            subplot(fms.display.rc(1),fms.display.rc(2),ic+(ir-1)*fms.display.rc(2));
            %read and downsample by block averaging
            z=imread(cat(2,faces_mpi_database_path,fms.filenames{istim}));
            ndown=floor(size(z,2)/opts_faces.image_width);
            z=z(1:ndown*floor(size(z,1)/ndown),1:ndown*opts_faces.image_width,:);
            z=reshape(z,[ndown,size(z,1)/ndown,ndown,size(z,2)/ndown,size(z,3)]);
            z=uint8(squeeze(mean(mean(z,1),3)));
            image(z);
            imgs{isetup}{istim}=z;
            axis equal;
            axis tight;
            title(cat(2,'s',zpad(istim,2),' ',fms.typenames{istim}),'Interpreter','none');
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
        end
        axes('Position',[0.01,0.04,0.01,0.01]); %for text
        text(0,0,cat(2,tstring,' ',fms.name),'Interpreter','none','FontSize',10);
        axis off;
    end
    ifok=getinp('1 to proceed to make session files','d',[0 1]);
    if (ifok==0)
        if getinp('1 for new random numbers','d',[0 1])
            opts_faces.how_rand='shuffle';
        else
            opts_faces.how_rand='default';
        end
    end
end
%
disp('for session config generation:');
if_frozen_psg=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1]);
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
opts_psg_std=opts_psg; %keep the same starting point
%generate session data -- modified from psg_spokes_setup to allow for
%augmented stimuli, which will be replaced by lower-numbered stimuli
for isetup_ptr=1:length(faces_mpi_info.setup_choices)
    isetup=faces_mpi_info.setup_choices(isetup_ptr);
    fms=faces_mpi_setups{isetup};
    %
    %customizations for mpi faces
    %
    opts_psg=opts_psg_std;
    opts_psg.cond_nstims=fms.nstims+fms.nstims_aug;
    opts_psg.cond_nstims_toreplace=fms.nstims_aug;
    opts_psg.example_infix_mode=4; %same example for each presentation, no infix
    %
    ifok=0;
    while (ifok==0)
        disp(sprintf('current psg session setup: %3.0f stimuli (%3.0f augmented), %3.0f comparison stimuli per trial, overlap %3.0f; %3.0f sessions',...
            opts_psg.cond_nstims,opts_psg.cond_nstims_toreplace,opts_psg.cond_ncompares,opts_psg.cond_novlp,opts_psg.cond_nsess));
        if (fms.nstims+fms.nstims_aug==opts_psg.cond_nstims)
            if getinp('1 if ok','d',[0 1])
                ifok=1;
            end
        else
            disp(sprintf('current setup invalid: number of stimuli in spokes: %3.0f, in psg: %3.0f',...
                nstims,opts_psg.cond_nstims));
        end
        if (ifok==0)
            opts_psg.cond_nstims=getinp('nstims (for session generation)','d',[1 1000],opts_psg.cond_nstims);
            opts_psg.cond_nstims_toreplace=getinp('nstims to replace (for session generation)','d',[1 1000],opts_psg.cond_nstims_toreplace);
            opts_psg.cond_ncompares=getinp('ncompares','d',[1 1000],opts_psg.cond_ncompares);
            opts_psg.cond_novlp=getinp('novlp','d',[1 1000],opts_psg.cond_novlp);
            opts_psg.cond_nsess=getinp('number of sessions','d',[1 100],opts_psg.cond_nsess);
        end
    end
    for k=1:length(opts_psg.refseq_labels)
        disp(sprintf('%1.0f->method for choosing stimuli in overlap: %s',k,opts_psg.refseq_labels{k}));
    end
    opts_psg.refseq=getinp('choice','d',[1 length(opts_psg.refseq_labels)],opts_psg.refseq); 
    [sessions_withaug,sessions_withaug_sorted,opts_used,psg_desc]=psg_sessconfig_make(opts_psg);
    opts_psg.cond_desc=psg_desc;
    s_all{isetup}.if_frozen_psg=if_frozen_psg;
    %
    %accumulate and display statistics of the augmented configuration
    %
    disp(sprintf('Analyzing the session configuration %s prior to replacement',psg_desc));
    session_stats_withaug=psg_session_stats(sessions_withaug,setfield(opts_psg,'if_log',1));
    s_all{isetup}.sessions_withaug=sessions_withaug;
    s_all{isetup}.session_stats_withaug=session_stats_withaug;
    %
    %replace the extra stimuli and do stats
    [sessions,opts_psg,warnings]=psg_sessconfig_replace(sessions_withaug,opts_psg);
    if ~isempty(warnings)
        disp(warnings);
    end
    disp(sprintf('Analyzing the session configuration %s after replacement',psg_desc));
    session_stats=psg_session_stats(sessions,setfield(opts_psg,'if_log',1));
    s_all{isetup}.sessions=sessions;
    s_all{isetup}.session_stats=session_stats;
    s_all{isetup}.opts_psg=opts_psg;
    %
    [session_cells,perms_used,examps_used]=psg_cond_create(sessions,fms.typenames,opts_psg);
    s_all{isetup}.session_cells=session_cells;
    s_all{isetup}.perms_used=perms_used;
    s_all{isetup}.examps_used=examps_used;
    %       
    disp('key variables')
    disp(s_all{isetup})
    %
    s_all{isetup}.creation_time=datestr(now);
    filename_base=getinp('file name base, _sess[#].csv will be appended for cond files','s',[],cat(2,cond_file_prefix,fms.name_brief));
    filename_mat=cat(2,pathname,filename_base,'.mat');
    s=s_all{isetup};
    save(filename_mat,'s');
	disp(sprintf('key variables saved in %s',filename_mat));
    %
    for igyfc=1:2 %gray level or full color
        disp(sprintf(' creating files for igyfc=%s (%s)',gyfc_infix{igyfc},gyfc_labels{igyfc}));
            %
        for isess=1:opts_psg.cond_nsess
            filename=cat(2,pathname,filename_base,'_',gyfc_infix{igyfc},'_sess',zpad(isess,opts_psg.sess_zpad));
            psg_cond_write(filename,session_cells{isess},setfields(opts_psg,{'if_log','cond_write_prefix'},...
                {1,cat(2,gyfc_infix{igyfc},'_')}));
        end
        %create and save stimulus files
        nexamps_expt=length(unique(examps_used(:)));
        disp(sprintf(' a maximum of %3.0f unique examples are needed for each of %2.0f stimuli (%s)',...
            nexamps_expt,fms.nstims,gyfc_labels{igyfc}));
        for istim=1:fms.nstims
            switch igyfc
                case 1
                    image_data=rgb2gray(imgs{isetup}{istim});
                case 2
                    image_data=imgs{isetup}{istim};
            end
            filename_stim=cat(2,pathname,gyfc_infix{igyfc},'_',fms.typenames{istim});
            filename_stim_full=cat(2,filename_stim,'.',opts_psg.stim_filetype);
            imwrite(image_data,filename_stim_full);
            disp(sprintf(' wrote stimulus %2.0f: %s',istim,filename_stim_full));
        end
    end
end

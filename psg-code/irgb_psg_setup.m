%irgb_psg_setup: sets up perceptual space geometry with independently distributed rgb stimuli (irgb)
%
%  Derived from faces_mpi_psg_setup amd psg_spokes_setup, but customized for irgb.  Some
%  features of psg_spokes_setup borrowed, as each stimulus typically has multiple examples
%
% stimulus file names: irgb_[paradigmID]_sXX_YYY.png, where XX is stimulus
% ID (e.g., 01 to 25), and YYY is stimulus example (000 to ?)
% 
% opts_spec, if specified, determines the stimulus set, otherwise defaults
% are taken from irgb_spec_make
%
% See also:  PSG_SPOKES_SETUP, FACES_MPI_PSG_SETUP, PSG_DEFOPTS, IRGB_SPEC_MAKE, IRGB_STIM_MAKE.
% PSG_COND_CREATE, PSG_COND_WRITE, PSG_SESSCONFIG_MAKE, PSG_SESSION_STATS.
%
nrgb=3;
%psg defaults
%
if ~exist('opts_psg') opts_psg=struct; end
opts_psg=psg_defopts(opts_psg);
if ~exist('opts_spec') opts_spec=struct; end
if ~exist('opts_stim') opts_stim=struct; end
if ~exist('nchecks') nchecks=16; end
if ~exist('nrep_display') nrep_display=20; end %block replication for display
%
%defaults for faces
if ~exist('cond_file_prefix') cond_file_prefix='irgb';end
%
ifok=0;
while ifok==0
    opts_stim_used=cell(0);
    [s,opts_spec_used]=irgb_spec_make(opts_spec);
    s.nchecks=nchecks;
    nchecks=getinp('nchecks','d',[4 128],nchecks);
    %
    nmean_dirs=size(s.opts_irgb.mean_dirs,1);
    nmean_mults=length(s.opts_irgb.mean_mults);
    nmeans=nmean_dirs*nmean_mults+s.opts_irgb.mean_include_zero;
    mean_include_zero=s.opts_irgb.mean_include_zero;
    ncovs=length(s.opts_irgb.cov_mults);
    nstims=ncovs*nmeans;
    stim_example=zeros(nchecks,nchecks,nrgb,nstims);
    %
    %display images and create session files
    %
    tstring=sprintf('name %s, type %s: cov_mode=%s, ifz=%1.0f offset=[%5.3f %.3f %5.3f]',...
        s.paradigm_name,s.paradigm_type,s.opts_irgb.cov_mode,s.opts_irgb.mean_include_zero,s.opts_irgb.mean_offset);
    figure;
    set(gcf,'Position',[50 150 1200 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',tstring);
    nc=nmean_mults*ncovs;
    nr=nmean_dirs+mean_include_zero;
    for istim=1:nstims
        switch s.paradigm_type
            case 'spokes'
                %determine row and column placement:
                %each row is a direction, columns are covariances cycled within mean multipliers 
                icov=ceil(istim/nmeans);
                imean=mod(istim-1,nmeans)+1;
                imult=mod(imean-1,nmean_mults)+1;
                idir=ceil(imean/nmean_mults);
                irow=idir;
                icol=icov+(imult-1)*ncovs;
                tshort=s.spec_labels{istim};
                tshort=strrep(tshort,'meandir','md');
                tshort=strrep(tshort,'meanmult','mm');
        end
        subplot(nr,nc,icol+(irow-1)*nc);
        [stim_example(:,:,:,istim),opts_stim_used{istim}]=irgb_stim_make(s.specs{istim},nchecks,1,opts_stim);
        imshow((1+repblk(stim_example(:,:,:,istim),[nrep_display,nrep_display,1,1]))/2);
        axis equal;
        axis tight;
        title(tshort);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        xlabel(sprintf('stim %1.0f',istim));
    end
    axes('Position',[0.01,0.04,0.01,0.01]); %for text
    text(0,0,tstring,'Interpreter','none','FontSize',10);
    axis off;
    ifok=getinp('1 if ok','d',[0 1]);
end
%section modified from psg_spokes_setup
%
% create psg spoke_setup
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
    disp(sprintf('current psg spoke_setup: %3.0f stimuli, %3.0f comparison stimuli per trial, overlap %3.0f; %3.0f sessions',...
        opts_psg.cond_nstims,opts_psg.cond_ncompares,opts_psg.cond_novlp,opts_psg.cond_nsess));
    if (nstims==opts_psg.cond_nstims)
        if getinp('1 if ok','d',[0 1])
            ifok=1;
        end
    else
        disp(sprintf('current setup invalid: number of stimuli in spokes: %3.0f, in psg: %3.0f',...
            nstims,opts_psg.cond_nstims));
    end
    if (ifok==0)
        opts_psg.cond_nstims=getinp('nstims','d',[1 1000],opts_psg.cond_nstims);
        opts_psg.cond_ncompares=getinp('ncompares','d',[1 1000],opts_psg.cond_ncompares);
        opts_psg.cond_novlp=getinp('novlp','d',[1 1000],opts_psg.cond_novlp);
        opts_psg.cond_nsess=getinp('number of sessions','d',[1 100],opts_psg.cond_nsess);
    end
end
for k=1:length(opts_psg.refseq_labels)
    disp(sprintf('%1.0f->method for choosing stimuli in overlap: %s',k,opts_psg.refseq_labels{k}));
end
opts_psg.refseq=getinp('choice','d',[1 length(opts_psg.refseq_labels)],opts_psg.refseq); 
[sessions,sessions_sorted,opts_used,psg_desc]=psg_sessconfig_make(opts_psg);
opts_psg.cond_desc=psg_desc;
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
[session_cells,perms_used,examps_used]=psg_cond_create(sessions,typenames,opts_psg);
%
s=struct;
s.nstims=nstims;
s.nchecks=nchecks;
s.nsubsamp=nsubsamp;
s.specs=specs;
s.spec_labels=spec_labels;
%
s.opts_psg=opts_psg;
s.typenames=typenames;
s.session_stats=session_stats;
s.sessions=sessions;
s.session_cells=session_cells;
s.perms_used=perms_used;
s.examps_used=examps_used;
%
s.btc_dict=dict;
s.btc_aug_opts=aug_opts;
s.btc_augcoords=augcoords;
s.btc_methods=methods;
s.if_frozen_btc=if_frozen_btc;
s.if_frozen_psg=if_frozen_psg;
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
filename_base=getinp('file name base (and path), _sess[#].csv will be appended for cond files','s',[]);
filename_mat=cat(2,pathname,filename_base,'.mat');
save(filename_mat,'s');
disp(sprintf('key variables saved in %s',filename_mat));
%
for isess=1:opts_psg.cond_nsess
    filename=cat(2,pathname,filename_base,'_sess',zpad(isess,opts_psg.sess_zpad));
    psg_cond_write(filename,session_cells{isess},setfield(opts_psg,'if_log',1));
end
%create and save stimulus files
nexamps_expt=length(unique(examps_used(:)));
disp(sprintf(' a maximum of %3.0f unique examples are needed for each stimulus',nexamps_expt));
for istim=1:nstims
    examples_reduced=btc_makemaps(methods{istim},setfields([],{'area','nmaps'},{nchecks,nexamps_expt}));
    for iexamp=1:nexamps_expt
        filename_stim=cat(2,typenames{istim},opts_psg.example_infix_string,zpad(iexamp-1,opts_psg.example_infix_zpad));
        filename_stim_full=cat(2,filename_stim,'.',opts_psg.stim_filetype);
        imwrite(pxlrep(examples_reduced(:,:,iexamp),nsubsamp),cat(2,pathname,filename_stim_full));
        if (iexamp==1) fn_first=filename_stim_full; end
        if (iexamp==nexamps_expt) fn_last=filename_stim_full; end
    end
    disp(sprintf(' stimulus type %2.0f: created %3.0f examples of %25s and wrote files %s to %s in %s',...
        istim,nexamps_expt,spec_labels{istim},fn_first,fn_last,pathname));
end

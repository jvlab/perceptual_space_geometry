%psg_spokes_setup: sets up perceptual space geometry with btc stimuli along spokes
%
%  Derived from spokes_layout_demo, but 
%   * has modified values of cmax
%   * does not do phase randomization or whitening or windowing; contrast fixed to 1
%   * invokes psg* modules to create cond files for Waraich's distance comparison protocol
%   * creates and writes png files for texture examples
%   * creates all files needed for a session (csv cond files, png image files, and
%     a mat file of parameters) saved in a user-entered directory
%
% 17Nov22:  allow for frozen randomization or not, for texture generation and for session config
% 08Dec22:  create max contrast lists in spokes_setup_create
%
% See also:  SPOKES_LAYOUT_DEMO, BTC_DEFINE, BTC_AUGCOORDS, BTC_MAKEMAPS, REPBLK, REPPXL, 
% PSG_DEFOPTS, PSG_COND_CREATE, PSG_COND_WRITE, PSG_SESSCONFIG_MAKE,
% PSG_SPEC2FILENAME, SPOKES_SETUP_CREATE, PSG_SHOWTRIAL, FACES_MPI_PSG_SETUP, IRGB_PSG_SETUP.
%
if_frozen_btc=getinp('1 for frozen random numbers, 0 for new random numbers each time for texture generation, <0 for a specific seed','d',[-10000 1]);
%
if (if_frozen_btc~=0)
    rng('default');
    if (if_frozen_btc<0)
        rand(1,abs(if_frozen_btc));
    end
else
    rng('shuffle');
end
%btc defaults
dict=btc_define;
nbtc=length(dict.codel);
aug_opts=struct;
aug_opts.ifstd=1;
%psg defaults
if ~exist('opts_psg')
    opts_psg=struct;
end
opts_psg=psg_defopts(opts_psg);
%
if ~exist('opts_stn')
    opts_stn=struct;
    opts_stn=filldefault(opts_stn,'decdig',3); %decimal digits in corr val in stimulus file name 3: cval=1->1000 
    opts_stn=filldefault(opts_stn,'sign_chars',{'m','z','p'}); %prefix characters for negative, zero, and positive cvals
    opts_stn=filldefault(opts_stn,'base',''); %start of stimulus file name
    opts_stn=filldefault(opts_stn,'rand','rand'); %name for random stimulus
end
%
if ~exist('spoke_setups')
    spokes_setup_create;
end
%
if ~exist('spoke_setup_choice') spoke_setup_choice=2; end
if ~exist('cmax_choice') cmax_choice=2; end
% already done in spokes_setup_create
% if ~exist('cmax_sets')
%     cmax_sets=cell(1);
%     cmax_sets{1}.desc='standard for thresholds';
%     cmax_sets{1}.vals=[0.2 0.4 0.4 0.5 0.5 1.0 1.0 1.0 1.0 0.8];
%     cmax_sets{2}.desc='standard for thresholds';
%     cmax_sets{2}.vals=[0.4 0.6 0.6 0.7 0.7 1.0 1.0 1.0 1.0 1.0];
% end
disp(cat(2,'      ',sprintf('%12s',dict.name_order_aug{:})));
for ic=1:length(cmax_sets)
    disp(sprintf(' set %2.0f: %s',ic,cmax_sets{ic}.desc));
    disp(cat(2,'      ',sprintf('%12.3f',cmax_sets{ic}.vals)));
end
cmax_choice=getinp('choice','d',[1 length(cmax_sets)],cmax_choice);
cmax=struct;
for ibtc=1:nbtc
    cmax.(dict.codel(ibtc))=cmax_sets{cmax_choice}.vals(ibtc);
end
%
disp('maximum values of stimuli')
disp(cmax)
%
if ~exist('nchecks') nchecks=16; end
if ~exist('nexamps') nexamps=1; end %number of examples of stimulus maps
if ~exist('nsubsamp') nsubsamp=16; end %subsampling of checks for display and phase reandomization
if ~exist('nclevs_max') nclevs_max=6; end
%figure layout options
if ~exist('axis_exspace_frac') axis_exspace_frac=0.2; end %extra spacing between axes, in units of axis length
if ~exist('axis_margin_frac') axis_margin_frac=0.5; end %margin, in units of axis length
%
for isetup=1:length(spoke_setups)
    disp(sprintf('%2.0f->%s',isetup,spoke_setups{isetup}.name));
end
%choose setup and btc params
spoke_setup_choice=getinp('choice','d',[1 length(spoke_setups)],spoke_setup_choice);
spoke_setup=spoke_setups{spoke_setup_choice};
if length(spoke_setup.btc_choices)>1
    for btc_choice=1:length(spoke_setup.btc_choices)
        btc_string=[];
        for ibtc=1:length(spoke_setup.btc_choices{btc_choice})
            btc_string=cat(2,btc_string,spoke_setup.btc_choices{btc_choice}{ibtc},' ');
        end
        disp(sprintf('%2.0f->%s',btc_choice,btc_string));
    end
    btc_choice=getinp('choice','d',[1 length(spoke_setup.btc_choices)]);
else
    btc_choice=1;
end
btc_choices=spoke_setup.btc_choices{btc_choice};
%
nclevs=getinp('number of c-levels','d',[1 nclevs_max],spoke_setup.nclevs);
clev_fracvals=getinp('c-levels as fraction of max along each spoke','f',[0 1],[1:nclevs]/nclevs); %multipliers for correlation values along each spoke
posit_fracvals=getinp('positions for display as fractions of max along each spoke','f',[0 1],clev_fracvals); %multipliers for positions along each spoke for display
nchecks=getinp('number of checks','d',[4 64],nchecks);
%
nstims=nclevs*spoke_setup.nspokes+1; %random is the last
stimlocs=zeros(nstims,2);
spokes=zeros(nstims,1);
cmults=zeros(nstims,1);
for ispoke=1:spoke_setup.nspokes
    for iclev=1:nclevs
        istim=iclev+(ispoke-1)*nclevs;
        spokes(istim)=ispoke;
        cmults(istim)=clev_fracvals(iclev);
        stimlocs(istim,:)=spoke_setup.endpoints(ispoke,:)*posit_fracvals(iclev);
    end
end
%
%compute closest approach to determine spacing
%
nearest=Inf;
for istim=2:nstims
    for jstim=1:istim-1
        nearest=min(nearest,max(abs(stimlocs(istim,:)-stimlocs(jstim,:))));
    end
end
axis_space_frac=axis_exspace_frac+max(0,1/(nclevs*nearest)-1); %nclevs since default spacing is 1/nclevs
%
%determine texture specifications and labels
%
specs=cell(nstims,1);
spec_labels=cell(nstims,1);
for istim=1:nstims
    specs{istim}=struct;
    if (istim<nstims)
        for ibtc=1:size(spoke_setup.mixing,2)
            letcode=btc_choices{ibtc};
            mixval=spoke_setup.mixing(spokes(istim),ibtc);
            if (mixval~=0)
                cval=mixval*cmax.(letcode)*cmults(istim);
                specs{istim}.(letcode)=cval;
                spec_labels{istim}=cat(2,spec_labels{istim},sprintf('%s=%5.2f ',letcode,cval));
            end
        end
        spec_labels{istim}=deblank(spec_labels{istim});
    else
        spec_labels{istim}='random';
    end
end
%
%make texture samples for display
%
nexamps_disp=1;
btc_samples_reduced=zeros(nchecks,nchecks,nexamps_disp,nstims); %the binary texture
btc_samples_display=zeros(nchecks*nsubsamp,nchecks*nsubsamp,nexamps_disp,nstims); %scaled to [-1 1],subsampled, and windowed 
augcoords=cell(nstims,1);
methods=cell(nstims,1);
typenames=cell(nstims,1);
for istim=1:nstims
    typenames{istim}=psg_spec2filename(specs{istim},opts_stn);
    disp(sprintf(' creating %4.0f stimuli of stimulus type %2.0f, specification %25s typename %s',...
        nexamps_disp,istim,spec_labels{istim},typenames{istim}))
    augcoords{istim}=btc_augcoords(specs{istim},dict,aug_opts);
    methods{istim}=augcoords{istim}.method{1};
    opts=struct;
    btc_samples_reduced(:,:,:,istim)=btc_makemaps(methods{istim},setfields([],{'area','nmaps'},{nchecks,nexamps}));
end
btc_samples_display=repblk(2*btc_samples_reduced-1,[nsubsamp nsubsamp 1 1]); %subsample
%
%determine figure layout
axis_n1d=(1+2*nclevs);
axis_size=1/(axis_n1d+2*axis_margin_frac+(axis_n1d-1)*axis_space_frac);
haxes=cell(nstims,1);
positions=zeros(nstims,4);
for istim=1:nstims
    offsets=stimlocs(istim,:)*(axis_n1d-1)*(1+axis_space_frac);
    positions(istim,:)=[(1+offsets*axis_size-axis_size)/2,axis_size,axis_size];
end
spokes_label=sprintf('nchecks %3.0f nsubsamp %2.0f',nchecks,nsubsamp);
%
%show sample textures
%
figure;
set(gcf,'NumberTitle','off');
tstring=spoke_setup.name;
set(gcf,'Name',tstring);
set(gcf,'Position',[100 100 900 800]);
for istim=1:nstims
    axes('Position',positions(istim,:));
    stim=btc_samples_display(:,:,1,istim);
    imagesc(stim,[-1 1]);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    axis equal;
    axis tight;
    colormap gray;
    title(sprintf('s %2.0f: %s',istim,spec_labels{istim}));
    haxes{istim}=gca;
end
axes('Position',[0.02,0.02,0.01,0.01]); %for text
text(0,0,spokes_label,'Interpreter','none');
axis off;
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

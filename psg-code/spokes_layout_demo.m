%spokes_layout_demo:  demonstrate some stimulus layouts along spokes
%
% Optionally randomizes the phases (_pr), obtained by Fourier-transforming a larger sample of the texture, with cosine bell
% windowing and pad factor 1, and also computes a whitened (_pr) texture by setting all Fourier amplitudes to 1
%
% See also:  BTC_AUGCOORDS, BTC_MAKEMAPS, REPBLK, REPPXL,
% MLIS_WINDOW_SETUP, MLIS_RUN_RESTOREMV, SPOKES_PSG_SETUP, SPOKES_SETUP_CREATE.
%
if ~exist('cmax')
    cmax=struct; %default maximum values of texture stats
end
cmax=filldefault(cmax,'g',0.2);
cmax=filldefault(cmax,'b',0.4);
cmax=filldefault(cmax,'c',0.4);
cmax=filldefault(cmax,'d',0.5);
cmax=filldefault(cmax,'e',0.5);
cmax=filldefault(cmax,'t',1.0);
cmax=filldefault(cmax,'u',1.0);
cmax=filldefault(cmax,'v',1.0);
cmax=filldefault(cmax,'w',1.0);
cmax=filldefault(cmax,'a',0.8);
%
if ~exist('nchecks') nchecks=16; end
if ~exist('nexamps') nexamps=1; end %number of examples of stimulus maps
if ~exist('nsubsamp') nsubsamp=16; end %subsampling of checks for display and phase reandomization
if ~exist('nclevs_max') nclevs_max=6; end
if ~exist('nchecks_pr') nchecks_pr=128; end % texture size for phase randomization, note that this is larger than nchecks
%   this size is also used for whitening
if ~exist('contrast_stdv_pr') contrst_stdv_pr=0.5; end %contrast if phase-randomized stimuli are used
if ~exist('fft2_padfactor') fft2_padfactor=2; end %no padding
%
if ~exist('if_prw') if_prw=0; end %default is not to make phase-randomized and whitened examples
%figure layout options
if ~exist('axis_exspace_frac') axis_exspace_frac=0.2; end %extra spacing between axes, in units of axis length
if ~exist('axis_margin_frac') axis_margin_frac=0.5; end %margin, in units of axis length
%
aug_opts=struct;
aug_opts.ifstd=1;
dict=btc_define;
%
nbtc=length(dict.codel);
%
if ~exist('setups')
    spokes_setup_create;
end
%
for isetup=1:length(spoke_setups)
    disp(sprintf('%2.0f->%s',isetup,spoke_setups{isetup}.name));
end
%choose setup and btc params
isetup=getinp('choice','d',[1 length(spoke_setups)]);
spoke_setup=spoke_setups{isetup};
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
nexamps=getinp('number of examples','d',[1 64],nexamps);
window_opts={'hard','cb1','cb2','gau'};
for iwindow=1:length(window_opts)
    disp(sprintf('%1.0f->%s',iwindow,window_opts{iwindow}));
end
iwindow=getinp('choice','d',[1 length(window_opts)],1);
window=mlis_window_setup(nchecks*nsubsamp,window_opts{iwindow});
%
if_prw=getinp('1 to make phase-randomized and whitened stimuli','d',[0 1],if_prw);
if (if_prw)
    nchecks_pr=getinp('nchecks for computing sample fft','d',[nchecks*2, 1024],nchecks_pr);
    fft2_padfactor=getinp('fft2 padfactor','d',[1 8],fft2_padfactor);
    nsubsamp=getinp('nsubsamp','d',[1 32],nsubsamp);
    contrast_stdv=getinp('contrast (stdv)','f',[.01 1],contrst_stdv_pr);
else
    contrast_stdv=getinp('contrast (stdv)','f',[.01 1],1);
end

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
%make texture samples
%
btc_samples_reduced=zeros(nchecks,nchecks,nexamps,nstims); %the binary texture
btc_samples_display=zeros(nchecks*nsubsamp,nchecks*nsubsamp,nexamps,nstims); %scaled to [-1 1],subsampled, and windowed 
if (if_prw)
    btc_samples_pr=zeros(nchecks*nsubsamp,nchecks*nsubsamp,nexamps,nstims);
    btc_samples_wh=zeros(nchecks*nsubsamp,nchecks*nsubsamp,nexamps,nstims);
    window_big=mlis_window_setup(nchecks_pr*nsubsamp,'cb2');
    center_select=[1:nchecks*nsubsamp]+nsubsamp*(nchecks_pr-nchecks)/2;%offset into middle of sample
    fft_size=nchecks_pr*nsubsamp*fft2_padfactor;
    fft_half=fft_size/2;
    max_imag_pr=0;
    max_real_pr=0;
    max_imag_wh=0;
    max_real_wh=0;
end
augcoords=cell(nstims,1);
methods=cell(nstims,1);
for istim=1:nstims
    disp(sprintf(' creating %4.0f stimuli of stimulus type %2.0f, specification %s',nexamps,istim,spec_labels{istim}))
    augcoords{istim}=btc_augcoords(specs{istim},dict,aug_opts);
    methods{istim}=augcoords{istim}.method{1};
    opts=struct;
    btc_samples_reduced(:,:,:,istim)=btc_makemaps(methods{istim},setfields([],{'area','nmaps'},{nchecks,nexamps}));
    if (if_prw)
        btc_working_all=btc_makemaps(methods{istim},setfields([],{'area','nmaps'},{nchecks_pr,nexamps}));
        for iexamp=1:nexamps
            btc_working=contrast_stdv*(2*btc_working_all(:,:,iexamp)-1); %set contrast and map to [-1 1]
            btc_working_mean=mean(btc_working(:));
            btc_working_var=var(btc_working(:));
            btc_working=window_big.*reppxl(btc_working-mean(btc_working(:)),nsubsamp); %subtract mean and replicate pixels
            btc_working_fft=fft2(btc_working,fft_size,fft_size);
            %
            %phase-randomized
            %
            %create random phases
            rand_phases=exp(2*pi*i*rand(fft_half+1,fft_size));
            rand_phases(1+[0 fft_half],1+[0 fft_half])=sign(real(rand_phases(1+[0 fft_half],1+[0 fft_half]))); %these components must be real
            rand_phases(1,1+fft_half+[1:fft_half-1])=conj(fliplr(rand_phases(1,[2:fft_half]))); %these components have to be complex_congjuages
            rand_phases=[rand_phases;conj(flipud(rand_phases(2:fft_half,:)))];
            rand_phases(fft_half+[1:fft_half],2:end)=fliplr(rand_phases(fft_half+[1:fft_half],2:end));
            rand_phases(1+fft_half,1+fft_half+[1:fft_half-1])=conj(fliplr(rand_phases(1+fft_half,[2:fft_half]))); %these components have to be complex_congjuages
            %apply them
            btc_working_ifft_pr=ifft2(btc_working_fft.*rand_phases); %apply phase randomization and invert the FFT
            max_real_pr=max(max_real_pr,max(abs(real(btc_working_ifft_pr(:)))));
            max_imag_pr=max(max_imag_pr,max(abs(imag(btc_working_ifft_pr(:)))));
            btc_working_center_pr=real(btc_working_ifft_pr(center_select,center_select)); %subsampled but not yet windowed
            btc_working_mv_pr=mlis_run_restoremv(btc_working_center_pr,btc_working_mean,btc_working_var); %restore mean and variance
            btc_samples_pr(:,:,iexamp,istim)=btc_working_mv_pr;
            %
            %whitened(based on each image's spectrum)
            %
            btc_working_mag=abs(btc_working_fft);
            btc_working_mag(btc_working_mag==0)=1;
            btc_working_ifft_wh=ifft2(btc_working_fft./btc_working_mag);
            max_real_wh=max(max_real_wh,max(abs(real(btc_working_ifft_wh(:)))));
            max_imag_wh=max(max_imag_wh,max(abs(imag(btc_working_ifft_wh(:)))));
            btc_working_center_wh=real(btc_working_ifft_wh(center_select,center_select));
            btc_working_mv_wh=mlis_run_restoremv(btc_working_center_wh,btc_working_mean,btc_working_var); %restore mean and variance
            btc_samples_wh(:,:,iexamp,istim)=btc_working_mv_wh;          
        end
    end
end
btc_samples_display=contrast_stdv*repblk(2*btc_samples_reduced-1,[nsubsamp nsubsamp 1 1]).*repmat(window,[1,1,nexamps nstims]); %scale to [-1 1],subsample, and window
if (if_prw)
    btc_samples_pr=btc_samples_pr.*repmat(window,[1,1,nexamps nstims]); %window
    btc_samples_wh=btc_samples_wh.*repmat(window,[1,1,nexamps nstims]); %window
    disp(sprintf(' maximum real and imaginary parts for phase-randomization: %15.12f %15.12f',max_real_pr,max_imag_pr));
    disp(sprintf(' maximum real and imaginary parts for           whitening: %15.12f %15.12f',max_real_wh,max_imag_wh));
end
%
%determine layout
axis_n1d=(1+2*nclevs);
axis_size=1/(axis_n1d+2*axis_margin_frac+(axis_n1d-1)*axis_space_frac);
haxes=cell(nstims,1);
positions=zeros(nstims,4);
for istim=1:nstims
    offsets=stimlocs(istim,:)*(axis_n1d-1)*(1+axis_space_frac);
    positions(istim,:)=[(1+offsets*axis_size-axis_size)/2,axis_size,axis_size];
end
spokes_label=sprintf('nchecks %3.0f contrast %4.2f pad factor %2.0f window %s nsubsamp %2.0f nchecks(pr&wh) %3.0f',...
    nchecks,contrast_stdv,fft2_padfactor,window_opts{iwindow},nsubsamp,nchecks_pr);
%
%show sample textures
%
for isample_type=1:1+2*if_prw
    figure;
    set(gcf,'NumberTitle','off');
    switch isample_type
        case 1
            tstring=spoke_setup.name;
        case 2
            tstring=cat(2,'phase-randomized ',spoke_setup.name);
        case 3
            tstring=cat(2,'whitened ',spoke_setup.name);
    end
    set(gcf,'Name',tstring);
    set(gcf,'Position',[100 100 900 800]);
    for istim=1:nstims
        axes('Position',positions(istim,:));
        switch isample_type
            case 1
                stim=btc_samples_display(:,:,1,istim);
            case 2
                stim=btc_samples_pr(:,:,1,istim);
            case 3
                stim=btc_samples_wh(:,:,1,istim);
        end
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
end

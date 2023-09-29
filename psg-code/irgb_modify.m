function [imgs,stats,opts_used]=irgb_modify(img_orig,opts)
%
% [imgs,stats,opts_used]=irgb_modify(img_orig,opts) modifies an original image
% via whitening, phase_scrambling, etc.
%
% img_orig: a square image, [npts npts 1 or 3]
% opts: options
%   opts.if_show: 1 to show image and modifications
%   opts.label: name for display
%   opts.range: min and max range for image, defaults to [0 1]
%      should be [0 65535] for png
%   opts.make_gray: 1 to convert to gray level (if size(img_orig,3)=3)
%   opts.figh: figure handle, can be empty
%   opts.nreplicates: number of replicates prior to FFT for filtering (defaults to 1)
%   opts.nrandphase: number of phase-randomized examples (defaults to 1)
%
% imgs: modified images.  Always rescaled to [0 1]
%   img.rescaled
%   img.gray: gray-level version
%   img.whitened: whitened
%   img.randphase: phase-randomized
% stats: statistics
% opts_used: options used
%
%   See also:  MLIS_RUN_RESTOREMV, IRGB_MODIFY_DEMO.
%
if (nargin<=1)
    opts=struct;
end
opts=filldefault(opts,'if_show',1);
opts=filldefault(opts,'label','image');
opts=filldefault(opts,'make_gray',1);
opts=filldefault(opts,'range',[0 1]);
opts=filldefault(opts,'figh',[]);
opts=filldefault(opts,'nreplicates',1);
opts=filldefault(opts,'nrandphase',1);
%
npts=size(img_orig,1);
%set up filter
freqs=[0:2*npts*opts.nreplicates-1];
cosbell=(1+cos(pi*freqs/npts/opts.nreplicates))/2;
cosbell2=cosbell.*cosbell';
%
imgs=struct;
%
imgs_use=double(img_orig-opts.range(1))/diff(opts.range);
imgs.rescaled=imgs_use;
%
imgs.gray=rgb2gray(imgs_use);
orig_mean=mean(imgs.gray(:));
orig_var=var(imgs.gray(:));
%
mirrored_h=[imgs.gray,fliplr(imgs.gray)];
%replicate
imgs.mirrored=repmat([mirrored_h;flipud(mirrored_h)],opts.nreplicates,opts.nreplicates);
%
img_fft2=fft2(imgs.mirrored);
imgs.fft_fftinv=real(ifft2(img_fft2));
img_filt=real(ifft2(cosbell2.*img_fft2));
imgs.filt=img_filt([1:npts],[1:npts]);
%
filt_mean=mean(imgs.filt(:));
filt_var=var(imgs.filt(:));
%
stats.orig_mean=orig_mean;
stats.orig_var=orig_var;
stats.filt_mean=filt_mean;
stats.filt_var=filt_var;
%
%flatten the power spectrum
fft_mag=sqrt(img_fft2.*conj(img_fft2));
fft_mag(fft_mag==0)=1;
img_whitened_mirrored=real(ifft2(img_fft2./fft_mag));
img_whitened=mlis_run_restoremv(img_whitened_mirrored(1:npts,1:npts),orig_mean,orig_var);
imgs.whitened=max(min(img_whitened,1),0);
%
img_whitened_filt_mirrored=real(ifft2(cosbell2.*img_fft2./fft_mag));
img_whitened_filt=mlis_run_restoremv(img_whitened_filt_mirrored(1:npts,1:npts),filt_mean,filt_var);
imgs.whitened_filt=max(min(img_whitened_filt,1),0);
%
%randomize the phases: ignore phase relationships needed for real, and then
%just take real part at the end
%
imgs.randphase=zeros(npts,npts,opts.nrandphase);
imgs.randphase_filt=zeros(npts,npts,opts.nrandphase);
%
for ir=1:opts.nrandphase
    rand_phases=exp(2*pi*i*rand(2*npts*opts.nreplicates,2*npts*opts.nreplicates));
    img_randphase_mirrored=real(ifft2(img_fft2.*rand_phases));
    img_randphase=mlis_run_restoremv(img_randphase_mirrored(1:npts,1:npts),orig_mean,orig_var);
    imgs.randphase(:,:,ir)=max(min(img_randphase,1),0);
    %
    img_randphase_filt_mirrored=real(ifft2(cosbell2.*img_fft2.*rand_phases));
    img_randphase_filt=mlis_run_restoremv(img_randphase_filt_mirrored(1:npts,1:npts),filt_mean,filt_var);
    imgs.randphase_filt(:,:,ir)=max(min(img_randphase_filt,1),0);
end
%
ncol=4+double(opts.nrandphase>1);
if (opts.if_show)
    if isempty(opts.figh)
        opts.figh=figure;
    end
    figure(opts.figh);
    set(gcf,'Position',[100 100 1200 800]);
    set(gcf,'Name',opts.label);
    set(gcf,'NumberTitle','off');
    irgb_modify_show(imgs.rescaled,1,opts.label,ncol);
    irgb_modify_show(imgs.gray,2,cat(2,opts.label,' gray'),ncol);
    irgb_modify_show(imgs.whitened,3,cat(2,opts.label,' whitened'),ncol);
    irgb_modify_show(imgs.randphase(:,:,1),4,cat(2,opts.label,' rand phases'),ncol);
    %
    irgb_modify_show(imgs.filt,ncol+2,cat(2,opts.label,' filt'),ncol);
    irgb_modify_show(imgs.whitened_filt,ncol+3,cat(2,opts.label,' whitened, filt'),ncol);
    irgb_modify_show(imgs.randphase_filt(:,:,1),ncol+4,cat(2,opts.label,' rand phases, filt'),ncol);
    %
    if (opts.nrandphase>1)
        irgb_modify_show(imgs.randphase(:,:,end),5,cat(2,opts.label,' rand phases'),ncol);
        irgb_modify_show(imgs.randphase_filt(:,:,end),ncol+5,cat(2,opts.label,' rand phases, filt'),ncol);
    end
end
%
opts_used=opts;

%
return
function irgb_modify_show(img_use,isub,label,ncol)
subplot(2,ncol,isub);
if size(img_use,3)==1 %convert to rgb so we can show in a figure with a rgb colormap
    img_use=repmat(img_use,[1 1 3]);
end
imagesc(img_use,[0 1]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
axis tight;
axis square;
ht=title(label,'Interpreter','none');
return

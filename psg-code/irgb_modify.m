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
%
npts=size(img_orig,1);
%set up filter
freqs=[0:2*npts-1];
cosbell=(1+cos(pi*freqs/npts))/2;
cosbell2=cosbell.*cosbell';
%
imgs=struct;
%
imgs.use=(img_orig-opts.range(1))/diff(opts.range);
imgs.rescaled=imgs.use;
%
orig_mean=mean(imgs.rescaled(:));
orig_var=var(imgs.rescaled(:));
%
imgs.gray=rgb2gray(imgs.use);
%
mirrored_h=[imgs.gray,fliplr(imgs.gray)];
imgs.mirrored=[mirrored_h;flipud(mirrored_h)];
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
rand_phases=exp(2*pi*i*rand(2*npts,2*npts));
img_randphase_mirrored=real(ifft2(img_fft2.*rand_phases));
img_randphase=mlis_run_restoremv(img_randphase_mirrored(1:npts,1:npts),orig_mean,orig_var);
imgs.randphase=max(min(img_randphase,1),0);
%
img_randphase_filt_mirrored=real(ifft2(cosbell2.*img_fft2.*rand_phases));
img_randphase_filt=mlis_run_restoremv(img_randphase_filt_mirrored(1:npts,1:npts),filt_mean,filt_var);
imgs.randphase_filt=max(min(img_randphase_filt,1),0);
%
orig_mean=mean(imgs.rescaled(:));
orig_var=var(imgs.rescaled(:));
%
if (opts.if_show)
    if isempty(opts.figh)
        opts.figh=figure;
    end
    figure(opts.figh);
    set(gcf,'Position',[100 100 1200 800]);
    set(gcf,'Name',opts.label);
    set(gcf,'NumberTitle','off');
    irgb_modify_show(imgs.rescaled,1,opts.label);
    irgb_modify_show(imgs.gray,2,cat(2,opts.label,' gray'));
    %irgb_modify_show(imgs.mirrored,3,cat(2,opts.label,' mirrored'));
    %irgb_modify_show(imgs.fft_fftinv,4,cat(2,opts.label,'fft and inv'));
    irgb_modify_show(imgs.whitened,3,cat(2,opts.label,' whitened'));
    irgb_modify_show(imgs.randphase,4,cat(2,opts.label,' rand phases'));
    %
    irgb_modify_show(imgs.filt,6,cat(2,opts.label,' filt'));
    irgb_modify_show(imgs.whitened_filt,7,cat(2,opts.label,' whitened, filt'));
    irgb_modify_show(imgs.randphase_filt,8,cat(2,opts.label,' rand phases, filt'));
end
%
opts_used=opts;

%
return
function irgb_modify_show(img_use,isub,label)
subplot(2,4,isub);
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


    




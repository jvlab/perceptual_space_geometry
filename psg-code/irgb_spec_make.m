function [s,ou]=irgb_spec_make(opts)
%[s,ou]=irgb_spec_make creates many of the metadata fields of a structure s
%for an iid color experiment.
%
% rgb values are specified with 0=global mean (gray), -1=lowest possible, +1=highest possible
%
% opts:options, can be omitted 
%  paradigm_name: paradigm name, no blanks or underscores, defaults to 'irgbtest'
%  mean_dirs: rgb triplets for maximal modulation along rays
%  mean_offset: rgb triplet to offset mean_dirs*mean_mults, defaults to [0 0 0]
%  mean_mults: scalars to multiply the rgb triplets
%  mean_incluce_zero: 1 to include zero mean (plus offset), 0 does not include.  defaults to 1
%  cov_mode: shape of covariance
%  cov_mults: magnitude of covariance
%  discrete: to specify discrete component (future)
%
%  See also: PSG_SPOKES_SETUP, IRGB_STIM_MAKE.
%
if nargin<=0
    opts=struct;
end
nrgb=3;
%
opts=filldefault(opts,'paradigm_name','irgbtest');
opts=filldefault(opts,'paradigm_type','spokes');
opts=filldefault(opts,'mean_dirs',[1 1 1;1 0 0;0 1 0;0 0 1]);
opts=filldefault(opts,'mean_offset',[0 0 0]);
opts=filldefault(opts,'mean_mults',[-0.75 -0.50 -0.25 0.25 0.50 0.75]);
opts=filldefault(opts,'mean_include_zero',1);
opts=filldefault(opts,'cov_mode','gaussian'); %gaussian=gaussian with specified covariance, 'ellipsoid'=ellipsoid with specified covariance
opts=filldefault(opts,'cov_mults',0.01); % covariance magnitude (square of std dev)
opts=filldefault(opts,'discrete','none'); %for discrete component
%
ou=opts;
%
s.paradigm_name=opts.paradigm_name;
s.paradigm_type=opts.paradigm_type;
s.opts_irgb=opts;
s.creation_time=datestr(now);
%
nmean_dirs=size(opts.mean_dirs,1);
nmean_mults=length(opts.mean_mults);
nmeans=nmean_dirs*nmean_mults+opts.mean_include_zero;
ncovs=length(opts.cov_mults);
nstims=ncovs*nmeans;
opts.nstims=nstims;
s.specs=cell(nstims,1);
s.typenames=cell(nstims,1);
s.spec_labels=cell(nstims,1);
%
istim=0;
for icov=1:ncovs
    cov_mult=opts.cov_mults(icov);
    for imean_dir=1:nmean_dirs
        mean_dir=opts.mean_dirs(imean_dir,:);
        for imean_mult=1:nmean_mults
            mean_mult=opts.mean_mults(imean_mult);
            mean_val=mean_dir*mean_mult+opts.mean_offset;
            istim=istim+1;
            s.specs{istim}.mean_val=mean_val;
            s.specs{istim}.cov=cov_mult*eye(nrgb);
            s.typenames{istim}=sprintf('cov%1.0f_meandir%1.0f_meanmult%1.0f',icov,imean_dir,imean_mult); %typenames cannot have a blank
            s.spec_labels{istim}=sprintf('cov=%1.0f meandir=%1.0f meanmult=%1.0f',icov,imean_dir,imean_mult);
        end
    end
    for iz=1:opts.mean_include_zero
        istim=istim+1;
        s.specs{istim}.mean_val=opts.mean_offset;
        s.specs{istim}.cov=cov_mult*eye(nrgb);
        s.typenames{istim}=sprintf('cov%1.0f_meandir%1.0f_meanmult%1.0f',icov,0,0);
        s.spec_labels{istim}=sprintf('cov=%1.0f meanzero',icov);
    end
end
for istim=1:length(s.specs)
    s.specs{istim}=filldefault(s.specs{istim},'cov_mode',opts.cov_mode);
    s.specs{istim}=filldefault(s.specs{istim},'discrete',opts.discrete);
end
return

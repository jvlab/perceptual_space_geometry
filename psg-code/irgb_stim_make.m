function [stims,ou]=irgb_stim_make(spec,nchecks,nexamps,opts)
% [stims,ou]=irgb_stim_make(spec,nchecks,nexamps,opts) creastes one or more examples of a stimulus for
% an iid color experiment
%
% rgb values are specified with 0=global mean (gray), -1=lowest possible, +1=highest possible
%
% spec: a strucrture, typically s.specs{k}, from irgb_stim_make. Fields include
%   mean_val: mean value, size [1 3]
%   cov: covariance size [3 3]
% nchecks: number of checks
% nexamps: number of examples, defaults to 1
% opts: options (for future)
%
% stims: array, size [nchecks nchecks 3 nexamps], entries are in range [-1 1]
% ou: options used
%  ou.ntrunc: number of values truncated on each channel, low and high (size [2 3])
%
%  See also: PSG_SPOKES_SETUP, IRGB_STIM_MAKE, GNORMCOR.
%
if nargin<=3
    opts=struct;
end
if nargin<=2
    nexamps=1;
end
nrgb=3;
%
npts=nchecks*nchecks*nexamps;
x=gnormcor(spec.cov,npts)';
rawvals=x+repmat(spec.mean_val,npts,1);
ntrunc=zeros(2,nrgb);
ntrunc(1,:)=sum(rawvals<-1,1);
ntrunc(2,:)=sum(rawvals>+1,1);
rawvals=min(+1,max(-1,rawvals));
stims=reshape(rawvals,[nchecks nchecks nexamps nrgb]);
stims=permute(stims,[1 2 4 3]);
ou.ntrunc=ntrunc;
return

function [stims,opts_stim_used]=irgb_stim_make(spec,nchecks,nexamps,opts_stim)
% [stims,opts_stim_used]=irgb_stim_make(spec,nchecks,nexamps,opts_stim) creastes one or more examples of a stimulus for
% an iid color experiment
%
% rgb values are specified with 0=global mean (gray), -1=lowest possible, +1=highest possible
%
% spec: a strucrture, typically s.specs{k}, from irgb_stim_make. Fields include
%   mean_val: mean value, size [1 3]
%   cov: covariance size [3 3]
%   transform2rgb: b,m,a describe how the values specified by mean, cov, discrete are transformed into rgb values
%    rgb_vals=(rawvals-a)*m+b
% nchecks: number of checks
% nexamps: number of examples, defaults to 1
% opts_stim: options (for future)
%
% stims: array, size [nchecks nchecks 3 nexamps], entries are in range [-1 1]
% opts_stim_used: options used
%  opts_stim_used.ntrunc: number of values truncated on each channel, low and high (size [2 3])
%
%  See also: PSG_SPOKES_SETUP, IRGB_STIM_MAKE, GNORMCOR, ELLIPCOR.
%
if nargin<=3
    opts_stim=struct;
end
if nargin<=2
    nexamps=1;
end
nrgb=3;
%
npts=nchecks*nchecks*nexamps;
switch spec.paradigm_type
    case 'spokes'
        switch spec.cov_mode
            case 'gaussian'
                x=gnormcor(spec.cov,npts)';
                rawvals=x+repmat(spec.mean_val,npts,1);
            case 'ellipsoid'
                x=ellipcor(spec.cov,npts)';
                rawvals=x+repmat(spec.mean_val,npts,1);
            otherwise
                rawvals=repmat(spec.mean_val,npts,1);
                warning(sprintf('unknown cov_mode )%s)',spec.cov_mode));
        end
    case 'distributions'
        weights=spec.distribution_weights;
        vals=spec.distribution_vals;
        nweights=length(weights);
        %
        rand_indices=1+sum(double(repmat(rand(npts,1),[1 nweights])>repmat(cumsum(weights)/sum(weights),[npts 1])),2);
        rawvals=vals(rand_indices,:);
    otherwise
        warning(sprintf('unknown paradigm type (%s)',spec.paradigm_type));
end
%transform
if isfield(spec,'transform2rgb')
    t=spec.transform2rgb;
    rawvals=(rawvals-t.a)*t.m+t.b;
end
ntrunc=zeros(2,nrgb);
ntrunc(1,:)=sum(rawvals<-1,1);
ntrunc(2,:)=sum(rawvals>+1,1);
rawvals=min(+1,max(-1,rawvals));

stims=reshape(rawvals,[nchecks nchecks nexamps nrgb]);
stims=permute(stims,[1 2 4 3]);
opts_stim_used.ntrunc=ntrunc;
return

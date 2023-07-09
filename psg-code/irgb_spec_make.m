function [s,spec_params_used]=irgb_spec_make(spec_params)
%[s,spec_params_used]=irgb_spec_make(spec_params) creates the stimlulus specifications in a structure s
% for an iid color experiment
%
% rgb values are specified with 0=global mean (gray), -1=lowest possible, +1=highest possible
%
% spec_params: parameters, can be omitted
%  *  parameters global to s
%  paradigm_name: paradigm name, no blanks or underscores, defaults to 'irgb_test', should begin with irgb_
%  paradigm_type: paradigm type, no blanks or underscores, defaults to 'spokes'
%  *  parameters copied into each substructure s.spec{k} but can be over-ridden and may depend on param_type
%  cov_mode: shape of covariance
%  discrete: to specify discrete component (future)
%  transform2rgb: b,m,a describe how the values specified by mean, cov, discrete are transformed into rgb values via rgb_vals=(rawvals-a)*m+b
%     transform2rgb.label is a free-text label
%  * params that do not appear in specs but are used to calculate mean and covariance of continuous specification
%  mean_dirs: triplets for maximal modulation along rays
%  mean_offset: triplet to offset mean_dirs*mean_mults, defaults to [0 0 0]
%  mean_mults: scalars to multiply the rgb triplets
%  mean_incluce_zero: 1 to include zero mean (plus offset), 0 does not include.  defaults to 1
%  cov_mults: magnitude of covariance
%
% s: a metadata structure; has fields for all the global parameters and also
%     s.typenames{istim}: used in csv files to specif stimulus image files, and for analysis, must be usable as a file name and cannot contain blanks
%     s.spec_labels{istim}: friendlier version of typenames, suitble for plot legends
%     s.creation_time
%     s.spec_params
%     s.specs{istim}, for each of the stimulus classes
%       has cov_mode, discrete, transform2rgb and
%       s.specs{istim}.mean (size [1 3]): mean value prior to transformation by transform2rgb
%       s.specs{istim}.cov  (size [3 3]): covariance, prior to transformation by transform2rgb, shape detemined by cov_mode
%  spec_params_used: spec_parameters used
%
%  See also: PSG_SPOKES_SETUP, IRGB_STIM_MAKE.
%
if nargin<=0
    spec_params=struct;
end
nrgb=3;
%
transform2rgb_def=struct;
transform2rgb_def.a=zeros(1,nrgb);
transform2rgb_def.m=eye(nrgb);
transform2rgb_def.b=zeros(1,nrgb);
transform2rgb_def.label='none';
%global params
spec_params=filldefault(spec_params,'paradigm_name','irgb_test');
spec_params=filldefault(spec_params,'paradigm_type','spokes');
%params used in each of specs{istim}
spec_params=filldefault(spec_params,'discrete','none'); %for discrete component
spec_params=filldefault(spec_params,'transform2rgb',transform2rgb_def); %transformation
spec_params=filldefault(spec_params,'cov_mode','gaussian'); %gaussian=gaussian with specified covariance, 'ellipsoid'=ellipsoid with specified covariance
%params used to calculate mean and covariance
spec_params=filldefault(spec_params,'mean_dirs',[1 1 1;1 0 0;0 1 0;0 0 1]);
spec_params=filldefault(spec_params,'mean_offset',[0 0 0]);
spec_params=filldefault(spec_params,'mean_mults',[-0.75 -0.50 -0.25 0.25 0.50 0.75]);
spec_params=filldefault(spec_params,'mean_include_zero',1);
spec_params=filldefault(spec_params,'cov_mults',0.01); % covariance magnitude (square of std dev)
%
s.paradigm_name=spec_params.paradigm_name;
s.paradigm_type=spec_params.paradigm_type;
s.creation_time=datestr(now);
%
switch spec_params.paradigm_type
    case 'spokes'
        nmean_dirs=size(spec_params.mean_dirs,1);
        nmean_mults=length(spec_params.mean_mults);
        nmeans=nmean_dirs*nmean_mults+spec_params.mean_include_zero;
        ncovs=length(spec_params.cov_mults);
        nstims=ncovs*nmeans;
    otherwise
        nstims=0;
end
s.specs=cell(nstims,1);
s.typenames=cell(nstims,1);
s.spec_labels=cell(nstims,1);
%
spec_params.nstims=nstims;
spec_params_used=spec_params;
s.spec_params=spec_params;
%
switch spec_params.paradigm_type
    case 'spokes'
        istim=0;
        for icov=1:ncovs
            cov_mult=spec_params.cov_mults(icov);
            for imean_dir=1:nmean_dirs
                mean_dir=spec_params.mean_dirs(imean_dir,:);
                for imean_mult=1:nmean_mults
                    mean_mult=spec_params.mean_mults(imean_mult);
                    mean_val=mean_dir*mean_mult+spec_params.mean_offset;
                    istim=istim+1;
                    s.specs{istim}.mean_val=mean_val;
                    s.specs{istim}.cov=cov_mult*eye(nrgb);
                    s.typenames{istim}=sprintf('cov%1.0f_meandir%1.0f_meanmult%1.0f',icov,imean_dir,imean_mult); %typenames cannot have a blank
                    s.spec_labels{istim}=sprintf('cov=%1.0f meandir=%1.0f meanmult=%1.0f',icov,imean_dir,imean_mult);
                end
            end
            for iz=1:spec_params.mean_include_zero
                istim=istim+1;
                s.specs{istim}.mean_val=spec_params.mean_offset;
                s.specs{istim}.cov=cov_mult*eye(nrgb);
                s.typenames{istim}=sprintf('cov%1.0f_meandir%1.0f_meanmult%1.0f',icov,0,0);
                s.spec_labels{istim}=sprintf('cov=%1.0f meanzero',icov);
            end
        end
    otherwise
        warning(sprintf('unknown paradigm type (%s)',spec_params.paradigm_type));
end
for istim=1:length(s.specs)
    s.specs{istim}=filldefault(s.specs{istim},'cov_mode',spec_params.cov_mode);
    s.specs{istim}=filldefault(s.specs{istim},'discrete',spec_params.discrete);
    s.specs{istim}=filldefault(s.specs{istim},'transform2rgb',spec_params.transform2rgb);
end
return

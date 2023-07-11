function spec_params=irgb_spec_defaults(spec_init)
% spec_params=irgb_spec_defaults(spec_init) initializes a stimulus specifications in a structure s
% for an iid color experiment
%
% rgb values are specified with 0=global mean (gray), -1=lowest possible, +1=highest possible
%
% spec_init: starting structure with non-default entries, can be zero
%  *  parameters global to s
%  paradigm_name: paradigm name, no blanks or underscores, defaults to 'irgb_test', should begin with irgb_
%  paradigm_type: paradigm type, no blanks or underscores, defaults to 'spokes'
%  *  parameters copied into each substructure s.spec{k} but can be over-ridden and may depend on param_type
%  [spokes]
%  cov_mode: shape of covariance
%  [distributions]
%  distribution_[weights|type|count]
%  [all]
%  transform2rgb: b,m,a describe how the values specified by mean, cov, discrete are transformed into rgb values via rgb_vals=(rawvals-a)*m+b
%     transform2rgb.label is a free-text label
%  * params that do not appear in specs but are used to calculate mean and covariance of continuous specification
%  [spokes]
%  mean_dirs: size [ndirs 3]: triplets for maximal modulation of mean along rays
%  mean_offset: triplet to offset mean_dirs*mean_mults, defaults to [0 0 0]
%  mean_mults: scalars to multiply the mean along each ray prior to offset
%  mean_incluce_zero: 1 to include zero mean (plus offset), 0 does not include.  defaults to 1
%  cov_mults: magnitude of covariance
%
% spec_params: initialized structure
%
%  See also: PSG_SPOKES_SETUP, IRGB_STIM_MAKE, IRGB_SPEC_DEFAULTS.
%
if nargin<=0
    spec_init=struct;
end
nrgb=3;
%
spec_params=spec_init;
%
transform2rgb_def=struct;
transform2rgb_def.a=zeros(1,nrgb);
transform2rgb_def.m=eye(nrgb);
transform2rgb_def.b=zeros(1,nrgb);
transform2rgb_def.label='none';
%global params
spec_params=filldefault(spec_params,'paradigm_name','irgb_test');
spec_params=filldefault(spec_params,'paradigm_type','spokes');
spec_params=filldefault(spec_params,'paradigm_type_list',{'spokes','distributions'});
%
%params used in each of specs{istim}
%for 'spokes'
spec_params=filldefault(spec_params,'cov_mode','gaussian'); %gaussian=gaussian with specified covariance, 'ellipsoid'=ellipsoid with specified covariance
spec_params=filldefault(spec_params,'cov_mode_list',{'gaussian','ellipsoid'});
%for 'distributions'
spec_params=filldefault(spec_params,'distribution_weights',ones(1,2)/2); %for distributions
spec_params=filldefault(spec_params,'distribution_type','random');
spec_params=filldefault(spec_params,'distribution_type_list',{'random','specified'});
spec_params=filldefault(spec_params,'distribution_count',25); 
%for all
spec_params=filldefault(spec_params,'transform2rgb',transform2rgb_def); %transformation
%
%params used to calculate mean and covariance
%for 'spokes'
spec_params=filldefault(spec_params,'mean_dirs',[1 1 1;1 0 0;0 1 0;0 0 1]);
spec_params=filldefault(spec_params,'mean_offset',[0 0 0]);
spec_params=filldefault(spec_params,'mean_mults',[-0.75 -0.50 -0.25 0.25 0.50 0.75]);
spec_params=filldefault(spec_params,'mean_include_zero',1);
spec_params=filldefault(spec_params,'cov_mults',0.01); % covariance magnitude (square of std dev)
%
return

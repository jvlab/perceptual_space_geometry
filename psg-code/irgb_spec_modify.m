function params_new=irgb_spec_modify(spec_params)
% params_new=irgb_spec_modify(spec_params) modifies a stimulus specification structure for independently distributed rgb stimuli (irgb)
%
% See also:  IRGB_PSG_SETUP, IRGB_SPEC_MAKE, IRGB_STIM_MAKE, IRGB_SPEC_DEFAULTS.
%
nrgb=3;
%
if (nargin<=0)
    spec_params=struct;
end
% fill defaults
spec_params=irgb_spec_defaults(spec_params);
spec_params_orig=spec_params;
ifok=0;
while (ifok==0)
    switch spec_params.paradigm_type
        case 'spokes'
            spec_params.mean_mults=getinp('mean multipliers (number of points on a ray)','f',[-1 1],spec_params.mean_mults);
        otherwise
            warning(sprintf('unknown paradigm type (%s)',spec_params.paradigm_type));
    end
    ifok=getinp('1 if ok, -1 to revert (discard changes)','d',[-1 1]);
    if (ifok==-1)
        spec_params=spec_params_orig;
    end
end
params_new=spec_params;
return

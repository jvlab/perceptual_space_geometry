function opts_fit=btc_soidfg_define(opts)
% opts_fit=btc_soidfg_define sets up an option s structure for figure-ground fits to psychophysical data
%
% modeled after structures used in btc_soid_fit and bt_soid_demo.
%
% opts: nondefault overrides
%
% opts_fit: options used for fitting
%   verbose:  1 for verbose terminal output
%   coords:  the coordinate set to include in the model ,e.g.,'gbcdetuvwa' or 'bcde'
%   nsurr:  number of surrogates to analyze (surrogates generated based on supplied error bars) defaults to 0
%   need_eb: 1 to only use thresholds with valid error bars, 0 to use all, defaults to 0
%   need_coords: 1 to use only the thresholds that are from planes whose statistics are included in opts_fit.coords (default)
%   ebpval:  error-bar p-value, defaults to 0.05, used to create surrogates
%   minvec: minimum absolute value for a vector to be used as a regressor
%   maxvec: maximum absolute value for a vector to be used as a regressor
%   th_min: minimum value for a non-exceptional threshold
%   th_max: maximum value for a non-exceptional threshold
%   eb_min: minimum value for a non-exceptional error bar
%   eb_max: maximum value for a non-exceptional error bar
%   eb_fill: how to fill in exceptional error bars, defaults to 1
%      0->do not fill
%      n->fill with a fraction of the threshold that is equal to the mean of
%      the n largest fractions (Inf means use all valid error bars in that plane)
%   ifaug: 1 to augment the coords (as is done to create the textures), 0 to not augment
%   sym_type: one of sym_type_avail
%   sym_type_avail: {'none','hflip','vflip','hvflip','rot90','rot180','diag','full'}
%   model_type: one of model_type_avail, defaults to 'f-g'
%   model_type_avail: {'f-g','f-g,f','f-g,g','f-g,f,g(linked)','f-g,f,g(indep)','allquad'}
%     options to restrict the model parameters:  all conditions to hold for terms to be kept. 
%       model_axes_probed and model_planes_probed are checked for each of the equivlent-by-symmetry names; if any are OK then they are kept..
%   model_exclude: cell array (1,n) of model parameters to exclude by name, defaults to cell(0), an entry of e.g., '_d_e' to exclude (d,e) cross-terms
%   model_axes_probed: cell array  (1,n) of axes probed, defaults to cell(0); an entry {'b'} {'c'} {'d'} means that only model terms that use these axes are kept
%   model_planes_probed: cell array (1,n) of planes probed, defaults to cell(0); an entry {'bc'} {'bd'} {'de'} means that only model terms that use these planes are kept
%  
% 27Apr20:  Add options to exclude based on axes and planes probed
%
%    See also:  BTC_SOIDFG_DEMO, BTC_SOID_FIT, BTC_SOID_DEMO, FIGGND_PROC, FIGGND_DBASE_PLOT,
%    BTC_SOIDFG_FIT, BTC_SOIDFG_MODEL, BTC_SYMCLASSES, BTC_SOIDFG_VALIDATE.
%
if (nargin==0)
    opts_fit=[];
else
    opts_fit=opts;
end
opts_fit=filldefault(opts_fit,'verbose',0);
opts_fit=filldefault(opts_fit,'coords','gbcdetuvwa'); %coordinates  to include in model
opts_fit=filldefault(opts_fit,'nsurr',0); % number of surrogates
opts_fit=filldefault(opts_fit,'need_eb',0); % 1 to require that error bars are present
opts_fit=filldefault(opts_fit,'need_coords',1); % 1 to require that thresholds use coordinates in opts_fit.coords
%
opts_fit=filldefault(opts_fit,'ebpval',0.05); %p-value corresponding to error bars, used for parameteric bootstrap
%
opts_fit=filldefault(opts_fit,'minvec',0.0001); %minimum threshold that can be used for a fit
opts_fit=filldefault(opts_fit,'maxvec',10); %maximum threshold that can be used for a fit
%
opts_fit=filldefault(opts_fit,'th_min',0.01); %minimum value for non-exceptional threshold
opts_fit=filldefault(opts_fit,'th_max',2); %maximum value for non-exceptional threshold
opts_fit=filldefault(opts_fit,'eb_min',0); %maximum value for non-exceptional error bar
opts_fit=filldefault(opts_fit,'eb_max',3); %maximum value for non-exceptional error bar
opts_fit=filldefault(opts_fit,'eb_fill',1); %how to fill in error bars
%
opts_fit=filldefault(opts_fit,'ifaug',1); %1 to fit based on augmented coordinate values (via btc_augcoords)
opts_fit=filldefault(opts_fit,'sym_type','full');
opts_fit=filldefault(opts_fit,'sym_type_avail',{'none','hflip','vflip','hvflip','rot90','rot180','diag','full'});
opts_fit=filldefault(opts_fit,'model_type','f-g');
opts_fit=filldefault(opts_fit,'model_type_avail',{'f-g','f-g,f','f-g,g','f-g,f,g(linked)','f-g,f,g(indep)','allquad'});
opts_fit=filldefault(opts_fit,'model_exclude',cell(0)); %list of name strings to exclude model parameters
opts_fit=filldefault(opts_fit,'model_axes_probed',cell(0)); %list of axes that must contain the model parameter coords
opts_fit=filldefault(opts_fit,'model_planes_probed',cell(0)); %list of planes (axis pairs) that must contain the model parameter coords
%
% used in btc_soid_fit but not btc_soidfg_*
%
%opts_fit=filldefault(opts_fit,'use_large',0); %set to 1 to allow use of out-of-range values
%opts_fit=filldefault(opts_fit,'symfull',1);
%opts_fit=filldefault(opts_fit,'consistency_tol',10^-5);
%opts_fit=filldefault(opts_fit,'consistency_chk',strvcat('uvecs','maxvecs','thresh_vecs'));
return

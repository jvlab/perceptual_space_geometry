function [specxs,avec,results,ou]=btc_soid_find_mix(specs,dsq,qform,dict_std,opts)
% [specxs,avec,results,ou]=btc_soid_find(specs,dsq,qform,dict_std,opts) finds the
% augmented vector in the direction of tangent vector spec that achieves
% a given value of a quadratic form dsq=avec*qform*avec'
%
% This is useful for finding vectors at a known distance from the origin,
% to create test datasets for btc_soid routines
%
% specs: a cell array of spec, a structure of up to 10 fields, chosen from the code letters gbcdetuvwa,
%    indicating the pairwise coordinates of the tangent vectors specifying the mixtures to be combined
%    Each of these is augmented by btc_augcoords, and the results added.
% dsq: the desired value of the quadratic form
% qform: a positive definite quadratic form, if empty, the identity is used
% dict_std: the dictionary, must be standard (fastest is to specify it and also to use 
%     opts.aug_opts.ifstd=1; if dict_std is empty, it is created with btc_define)
% opts: opts.aug_opts: options used for btc_augcoords, defaults to []
%
% specxs: the specification of the cell arrays of plane tangent vectors that, when augmented,
%    will achieve the desired value of the quadratic form.
% avec: the augmented vector that achieves the desired quadratic form
% results: outputs of the minimization routine
%    results.x is the multiplier: fields of specx = (fields of spec) *results.x 
%    results.fval= the fitted value, i.e., results.fval+dsq is the actual 
%      value of the quadratic form
% ou: options used
%
% Note that specxs and avec will have Inf components if the length element is zero or negative,
%  which can occur if the quadratic form qform is not positive-definite.
%
%    See also: BTC_SOID_FIND, BTC_DEFINE, BTC_AUGCOORDS_MIX, BTC_SOID_FIND_MIX_OF, FZERO, TEXRES_SETUP.
%
if (nargin<3)
    qform=[];
end
if isempty(qform)
    qform=eye(length(dict_std.codel));
end
if (nargin<4)
    dict_std=[];
end
if (isempty(dict_std))
    dict_std=btc_define;
end
if (nargin<5)
    opts=[];
end
opts=filldefault(opts,'aug_opts',setfield([],'ifstd',1));
opts.aug_opts.ifstd=1;
opts.aug_opts.nocheck=1;
opts.dict_std=dict_std;
aug_opts=opts.aug_opts;
%
results=[];
ou=opts;
%
vecsquared=btc_soid_find_mix_of(1,0,specs,qform,opts);
if (vecsquared>0)
    [x,fval,exitflag,output] = fzero(@(x) btc_soid_find_mix_of(x,dsq,specs,qform,opts),1);
else
    x=Inf;
    fval=0;
    exitflag=[];
    output=[];
end
if (x<0)
    x=Inf;
    fval=0;
    exitflag=[];
    output=[];
end
results.x=x; %this is the multipler for specs
results.fval=fval; %this should be near 0
results.exitflag=exitflag;
results.output=output;
%
specxs=cell(0);
for imix=1:length(specs)
    specxs{imix}=[];
    vars=fieldnames(specs{imix});
    for ivar=1:length(vars)
        specxs{imix}=setfield(specxs{imix},vars{ivar},x*getfield(specs{imix},vars{ivar}));
    end
end %imix
avec=btc_augcoords_mix(specxs,opts.dict_std,opts.aug_opts);
%
return

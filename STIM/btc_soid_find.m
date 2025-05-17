function [specx,avec,results,ou]=btc_soid_find(spec,dsq,qform,dict_std,opts)
% [specx,avec,results,ou]=btc_soid_find(spec,dsq,qform,dict_std,opts) finds the
% augmented vector in the direction of tangent vector spec that achieves
% a given value of a quadratic form dsq=avec*qform*avec'
%
% This is useful for finding vectors at a known distance from the origin,
% to create test datasets for btc_soid routines
%
% spec: a tangent vector, specified as a structure, with field names
%   given by the nonzero btc coordinate letters
% dsq: the desired value of the quadratic form
% qform: a quadratic form, typically positive-definite, if empty, the identity is used
% dict_std: the dictionary, must be standard (fastest is to specify it and also to use 
%     opts.aug_opts.ifstd=1; if dict_std is empty, it is created with btc_define)
% opts: opts.aug_opts: options used for btc_augcoords, defaults to []
%
% specx: the specification of the tangent vector that, when augmented,
%    will achieve the desired value of the quadratic form.
% avec: the augmented vector that achieves the desired quadratic form
% results: outputs of the minimization routine
%    results.x is the multiplier: fields of specx = (fields of spec) *results.x 
%    results.fval= the fitted value, i.e., results.fval+dsq is the actual 
%      value of the quadratic form
% ou: options used
%
% Note that specx and avec will have Inf components if the length element is zero or negative,
%  which can occur if the quadratic form qform is not positive-definite (8-Dec-11).
%
%    See also: BTC_DEFINE, BTC_AUGCOORDS, BTC_SOID_FIND_OF, FZERO, BTC_SOID_FIND_MIX, BTC_SOIDFG_FIND.
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
% use [x,fval,exitflag,output] = fzero(@(x) btc_soid_find_of(x,dsq,spec,qform,opts),1)
% conditional on length element being positive (8-Dec-11)
%
vecsquared=btc_soid_find_of(1,0,spec,qform,opts);
if (vecsquared>0)
    [x,fval,exitflag,output] = fzero(@(x) btc_soid_find_of(x,dsq,spec,qform,opts),1);
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
results.x=x; %this is the multipler for spec
results.fval=fval; %this should be near 0
results.exitflag=exitflag;
results.output=output;
%
specx=[];
vars=fieldnames(spec);
for ivar=1:length(vars)
    specx=setfield(specx,vars{ivar},x*getfield(spec,vars{ivar}));
end
augcoords=btc_augcoords(specx,opts.dict_std,opts.aug_opts);
avec=augcoords.method{1}.vec;
%
return

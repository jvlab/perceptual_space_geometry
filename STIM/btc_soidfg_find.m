function [x,specx,avec,results,ou]=btc_soidfg_find(specs,dsq,qform,dict_std,opts_fit,m0)
% [x,specx,avec,results,ou]=btc_soidfg_find(specs,dsq,qform,dict_std,opts_fit,m0) finds the
% augmented figure,ground vector in the direction of tangent figure, ground vector given by specs
% that achieves a given value of a quadratic form dsq=avec*qform*avec'
%
% Structure is similar to btc_soid_find, but arguments are quite different.
% * If no solution is found, x is Inf, specx is null, and avec is NaN. With btc_soid_find, specx and avec may have Inf components
% * As with btc_soid_find, this augments the specification vector (via btc_augcoords_mults), and should therefore
%    not be called if opts_fit.ifaug=0
%
% specs: specs{1} is the figure specification, specs{2} is groud specification
% dsq: the desired value of the quadratic form
% qform: a quadratic form, typically positive definite, if empty, the identity is used
% dict_std: the dictionary, must be standard (fastest is to specify it and also to use 
%     opts_fit.aug_opts.ifstd=1; if dict_std is empty, it is created with btc_define)
% opts_fit: from btc_soidfg_define, but also 
%   opts_fit.aug_opts: options used for btc_augcoords, defaults to []
% m0: initial guess for multiplier, calcated from non-augmented (figure,ground) vectors if not supplied.
%
% x: multiplier (schematically, specx=x*specs) 
% specx: specifications of (figure,ground) vector that, when augmented,
%    will achieve the desired value of the quadratic form.
%    specx{1} for figure, specx{2} for ground
% avec: the augmented vector that achieves the desired quadratic form
% results: outputs of the minimization routine
%    results.x is the multiplier: fields of specx = (fields of spec) *results.x 
%    results.fval= the fitted value, i.e., results.fval+dsq is the actual 
%      value of the quadratic form
% ou: options used
%
%    See also: BTC_DEFINE, BTC_SOIDFG_PREDICT, BTC_AUGCOORDS_MULTS, BTC_SOIDFG_FIND_OF, FZERO, BTC_SOID_FIND.
%
if (nargin<3)
    qform=[];
end
if isempty(qform)
    qform=eye(2*length(dict_std.codel));
end
if (nargin<4)
    dict_std=[];
end
if (isempty(dict_std))
    dict_std=btc_define;
end
if (nargin<5)
    opts_fit=[];
end
if (nargin<6)
    m0=[];
end
opts_fit=filldefault(opts_fit,'aug_opts',setfield([],'ifstd',1));
opts_fit.aug_opts.ifstd=1;
opts_fit.aug_opts.nocheck=1;
opts_fit.dict_std=dict_std;
aug_opts=opts_fit.aug_opts;
btc_n=length(dict_std.codel);
%
ou=opts_fit;
%
% use [x,fval,exitflag,output] = fzero(@(x) btc_soidfg_find_of(x,dsq,specs,qform,opts_fit),1)
% conditional on length element being positive (8-Dec-11)
%
vecsquared=btc_soidfg_find_of(1,0,specs,qform,opts_fit);
if (vecsquared>0)
    if isempty(m0)
        m0=sqrt(dsq/vecsquared);
    end
    [x,fval,exitflag,output] = fzero(@(x) btc_soidfg_find_of(x,dsq,specs,qform,opts_fit),m0);
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
% recover augmented vector and specifications for results of fzero
%
avec=NaN(1,2*btc_n);
specx=cell(1,2);
if isfinite(x)
    for ifg=1:2 %augment figure and ground separately
        specx{ifg}=btc_vecnan2letcode(x*btc_letcode2vec(specs{ifg},dict_std),dict_std);
        avec(btc_n*(ifg-1)+[1:btc_n])=btc_augcoords_mults(specs{ifg},x,dict_std,opts_fit.aug_opts);
    end
end
return

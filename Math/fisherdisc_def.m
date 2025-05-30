function optsdef=fisherdisc_def(opts)
% optsdef=fisherdisc_def(opts) sets up the options for fisherdisc
%
% opts: an options structure
% optsdef: options structure, with defaults filled in
%
% See fisherdisc.m for description of options
%
% See also:
%  FISHERDISC_TEST, FISHERDISC_DEF, FISHERDISC, FINDSEGFLIPS.
%
if (nargin<1)
    opts=[];
end
opts=filldefault(opts,'nshuffle',0);
opts=filldefault(opts,'nshuffle_max',200);
opts=filldefault(opts,'xvalid',0);
opts=filldefault(opts,'sub1',0);
opts=filldefault(opts,'classifiers',[]);
opts=filldefault(opts,'condmax',Inf);
opts=filldefault(opts,'segs',[]);
opts=filldefault(opts,'delrad',Inf);
opts=filldefault(opts,'nflips',200);
opts=filldefault(opts,'nflips_maxtries',10000);
opts=filldefault(opts,'nflips_tol',1);
optsdef=opts;
return

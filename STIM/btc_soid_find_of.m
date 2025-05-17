function d=btc_soid_find_of(x,dsq,spec,qform,opts)
% d=btc_soid_find_of(x,dsq,spec,qform,opts) is the objective function for btc_soid_find
%   x:scalar multipler
%   dsq: the desired squared distance
%   spec: a tangent-vector specification, as a structure of 10 btc parameters
%   qform: a quadratic form
%   opts:  options argument
%       opts.dict_std must be the standard dictionary
%       opts.aug_opts must be present
%       opts.aug_opts.ifstd=1
%
%   See also:  BTC_SOID_FIND, BTC_AUGCOORDS.
%
specx=[];
opts=filldefault(opts,'verbose',0);
vars=fieldnames(spec);
for ivar=1:length(vars)
    specx=setfield(specx,vars{ivar},x*getfield(spec,vars{ivar}));
end
augcoords=btc_augcoords(specx,opts.dict_std,opts.aug_opts);
vec=augcoords.method{1}.vec;
d=vec*qform*vec'-dsq;
if (opts.verbose)
    disp(sprintf(' btc_soid_find_of called with x=%10.5f and finds vec=',x));
    disp(sprintf(repmat(' %10.5f',1,10),vec)); 
end
return
 
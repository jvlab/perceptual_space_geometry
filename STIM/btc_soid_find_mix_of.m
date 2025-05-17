function d=btc_soid_find_mix_of(x,dsq,specs,qform,opts)
% d=btc_soid_find_of(x,dsq,specs,qform,opts) is the objective function for
% btc_soid_find_mix
%   x:scalar multipler
%   dsq: the desired squared distance
%   specs: a cell array of tangent-vector specifications, as a structure of 10 btc parameters
%   qform: a quadratic form
%   opts:  options argument
%       opts.dict_std must be the standard dictionary
%       opts.aug_opts must be present
%       opts.aug_opts.ifstd=1
%
%   See also:  BTC_SOID_FIND_MIX, BTC_AUGCOORDS_MIX, BTC_SOID_FIND_OF.
%
opts=filldefault(opts,'verbose',0);
specxs=cell(0);
for imix=1:length(specs)
    specxs{imix}=[];
    vars=fieldnames(specs{imix});
    for ivar=1:length(vars)
        specxs{imix}=setfield(specxs{imix},vars{ivar},x*getfield(specs{imix},vars{ivar}));
    end
end
vec=btc_augcoords_mix(specxs,opts.dict_std,opts.aug_opts);
d=vec*qform*vec'-dsq;
if (opts.verbose)
    disp(sprintf(' btc_soid_find_mix_of called with x=%10.5f and finds vec=',x));
    disp(sprintf(repmat(' %10.5f',1,10),vec)); 
end
return
 
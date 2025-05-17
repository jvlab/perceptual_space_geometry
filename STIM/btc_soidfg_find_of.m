function d=btc_soidfg_find_of(x,dsq,specs,qform,opts)
% d=btc_soid_findfg_of(x,dsq,spec,qform,opts) is the objective function for btc_soid_find
%   x:scalar multipler
%   dsq: the desired squared distance
%   specs: cell array of tangent vector specifications, specs{1} 10 btc parameters for figure, specs{2} for ground
%   qform: a quadratic form
%   opts:  options argument
%       opts.dict_std must be the standard dictionary
%       opts.aug_opts must be present
%       opts.aug_opts.ifstd=1
%
%   See also:  BTC_SOIDFG_FIND, BTC_AUGCOORDS, BTC_AUGCOORDS_MULTS.
%
specx=[];
opts=filldefault(opts,'verbose',0);
dict=opts.dict_std;
btc_n=length(dict.codel);
vec_fg=zeros(1,2*btc_n);
for ifg=1:2 %augment figure and ground separately
    vec_fg(btc_n*(ifg-1)+[1:btc_n])=btc_augcoords_mults(specs{ifg},x,dict,opts.aug_opts);
end
d=vec_fg*qform*vec_fg'-dsq;
if (opts.verbose)
    disp(sprintf(' btc_soidfg_find_of called with x=%10.5f and finds vec_fg=',x));
    disp(reshape(vec_fg',[btc_n 2])');
end
return
 
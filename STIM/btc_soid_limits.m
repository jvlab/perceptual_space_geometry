function opts_limits=btc_soid_limits(opts)
% opts_limits=btc_soid_limits(opts) sets up a structure of low and high
% limits for psychophysical data
%
% opts: input structure, can be omitted
%    opts.limits_X: [low high] limits for parameter X (X=ndirs,namps, ntrials, fraccor, thresh_mags, bexpon_mags)
%    opts.nan_X: 0 to cause a "not OK" message if limit is exceeded, 1 to replace by NaN
% opts_limits: opts, with limiting values added
%    opts_limits is opts, but with values omitted filled in
%
%   See also:  BTC_SOID_XLSREAD, BTC_SOID_ARCH2DS, BTC_SOID_GETDATA.
%
% sets up limits for sane data in reading psychophysical data from spreadsheets
if (nargin<1) opts=[]; end
opts=filldefault(opts,'limits_ndirs',[3 24]); %inserted and broadened, Sept 3 2017
opts=filldefault(opts,'limits_namps',[1 100]);
opts=filldefault(opts,'limits_ntrials',[1 10000]);
opts=filldefault(opts,'limits_fraccor',[0 1]);
opts=filldefault(opts,'limits_thresh_mags',[0 10]); %don't check error bars for sanity, as this is done in btc_soid_getdata
opts=filldefault(opts,'limits_bexpon_mags',[0 100]);
opts=filldefault(opts,'nan_namps',0);
opts=filldefault(opts,'nan_ntrials',0);
opts=filldefault(opts,'nan_fraccor',0);
opts=filldefault(opts,'nan_thresh_mags',1);
opts=filldefault(opts,'nan_bexpon_mags',1);
opts_limits=opts;
return

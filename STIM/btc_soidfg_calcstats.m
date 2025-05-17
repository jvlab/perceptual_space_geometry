function stats=btc_soidfg_calcstats(weib_a_data,weib_a_fit,nparams,opts_fit)
% stats=btc_soidfg_calcstats(weib_a_data,weib_a_fit,nparams,opts_fit) calculates statistics of a figure-ground fit
%
% weib_a_data: triplets (threshold, cl_lo, cl_hi) of threshold data
% weib_a_fit: a column vector of fitted thresholds, typically from btc_soidfg_predict
% nparams: number of model params
% opts_fit: options structure from btc_soidfg_define, filled with defaults if omitted
%  
% stats: statistics structure, with following fields
%
% zscore: rms of (model error/data s.d.)
% chi2: corresponding chi-squared value
% chi2_dof: number of degrees of freedom
% chi2_pval: p-value for chi2
% chi2reduced: reduced chi-square, equal to zscore*(number of values fit)/(nvals_fit-nparams)
% chi2reduced_dof: degrees of freedom in reduced chi-squared (number of values fit-nparams)
% chi2reduced_pval p-value for reduced chi-squared
% ss_[total,aboutmean,unexplained]: sum of squares of [thresholds,thresholds-mean thresholds, thresholds-fit]
% ss_ratio: ratio of unexplained ss to about-mean ss
% var_measurement: estimated variances of measuremnts
% var_unexpained: variances unexplained
% fratio: ratio of variance unexplained to variance of measurement
% fratio_pval: p-value for f-ratio
% fratio_dof: degrees of freedom for fratio
%
% 24May20:  Add ss_* var_*, fratio_* fields, logic from figgnd_plane_fit
%
% See also:
% BTC_SOIDFG_DEMO, BTC_SOIDFG_MODEL, BTC_SOIDFG_FIT, BTC_SOIDFG_PREDICT,
% BTC_SOIDFG_DEFINE, NORMINV, CHI2CDF, FIGGND_PLANE_FIT.
%
if nargin<=3
    opts_fit=[];
end
opts_fit=btc_soidfg_define(opts_fit);
clim2sd=2*norminv(1-opts_fit.ebpval/2); %divide by this to convert full confidence limit to 1 standard deviation
dmf=weib_a_fit-weib_a_data(:,1);
rmse=sqrt(mean((dmf(~isnan(dmf)).^2))); %root mean squared error
mabe=median(abs(dmf)); %median average error
stats.rmse=rmse;
stats.mabe=mabe;
sd_data=(weib_a_data(:,3)-weib_a_data(:,2))/clim2sd;
%
zval=dmf./sd_data;
stats.zscore=sqrt(mean(zval(~isnan(zval).^2)));
nvals_fit=length(~isnan(zval)); %number of values fit
stats.chi2=stats.zscore*nvals_fit;
stats.chi2_dof=nvals_fit;
stats.chi2_pval=1-chi2cdf(stats.chi2,nvals_fit);
stats.chi2reduced=stats.zscore*nvals_fit/(nvals_fit-nparams);
stats.chi2reduced_dof=nvals_fit-nparams;
stats.chi2reduced_pval=1-chi2cdf(stats.chi2,nvals_fit-nparams);
%
%ftest logic from figgnd_plane_fit
%
stats.ss_total=sum(weib_a_data(:,1).^2);
stats.ss_aboutmean=sum((weib_a_data(:,1)-mean(weib_a_data(:,1))).^2);
stats.ss_unexplained=sum((weib_a_data(:,1)-weib_a_fit).^2);
stats.ss_ratio=stats.ss_unexplained/stats.ss_aboutmean;
stats.var_measurement=sum(sd_data.^2)/nvals_fit;
stats.var_unexplained=stats.ss_unexplained/(nvals_fit-nparams);
stats.fratio=stats.var_unexplained/stats.var_measurement;
stats.fratio_dof=[nvals_fit-nparams,nvals_fit];
stats.fratio_pval=1-fcdf(stats.fratio,stats.fratio_dof(1),stats.fratio_dof(2)); %compare variance on basis of misfit, vs variance from confidence limits
return

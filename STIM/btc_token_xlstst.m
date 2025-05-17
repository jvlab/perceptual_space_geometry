% btc_token_xlstst tests the format of an xls sheet with token data, by checking data for
% all models
%
%  See also: BTC_TOKEN_ARCH2DS, BTC_SOID_XLSTST, BTC_SOID_GETDATA, BTC_SOID_XLSREAD,
%  BTC_SOID_LIMITS, BTC_SOID_XLSFMT, BTC_SOID_ARCH2DS.
%
if ~exist('xls_name') xls_name='..\iso\token\Token_weibull_analysis.xls'; end
if ~exist('xls_sheet') xls_sheet='MC_data'; end
% 
modelnames={'allraysfixedb'};
dataul_rows=[113];
opts_list.plane_spacing=182; %override a default in  BTC_SOID_XLSFMT.
opts_list.nplanes=3; %number of planes
%
ds_list=cell(0);
ou_list=cell(0);
%
pr=[];
for imodel=1:length(modelnames)
    disp(sprintf('reading model %s',modelnames{imodel}));
    opts_list=btc_soid_xlsfmt(opts_list);
    opts_list=setfield(opts_list,'fit_type',modelnames{imodel});
    opts_list=setfield(opts_list,'dataul_row',dataul_rows(imodel));
    [ds_list.(modelnames{imodel}),ou_list.(modelnames{imodel}),pr]=...
        btc_soid_arch2ds(xls_name,xls_sheet,opts_list,pr);
end

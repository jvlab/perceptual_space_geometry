% btc_token_xls2ds creates mat-files with dataset structures that
% can be read by btc_soid_getdata.
%
%  See also: BTC_TOKEN_XLSTST, BTC_SOID_XLSTST, BTC_SOID_GETDATA, BTC_SOID_XLSREAD,
%  BTC_SOID_LIMITS, BTC_SOID_XLSFMT, BTC_SOID_GETDATA, BTC_TOKEN_EDIRS.
%

if ~exist('xls_name') xls_name='..\iso\token\Token_weibull_analysis.xls'; end
if ~exist('xls_sheet') xls_sheet='MC_data'; end
if ~exist('modelnames') modelnames={'allraysfixedb'}; end
xls_sheet=getinp('xls sheet name','s',[0 1],xls_sheet);
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
    [ds,ou_list.(modelnames{imodel}),pr]=...
        btc_soid_arch2ds(xls_name,xls_sheet,opts_list,pr);
end
ifsave=getinp('1 to save file','d',[0 1],1);
if (ifsave==1)
    fn=strrep('token_allraysfixedb_*','*',lower(xls_sheet(1:2)));
    fn=getinp('file name','s',[0 1],fn);
    save(fn,'ds');
    disp(sprintf('%s can be read by [ed,eb_avail,ds,condno,ou]=btc_soid_getdata([])',fn));
    disp('caution: USE btc_token_edirs to retrieve stimulus directions.');
else
    disp('file not saved.');
end

function [ds,ou,pr]=btc_soid_arch2ds(xls_name,xls_sheet,opts,preread)
% [ds,ou,pr]=btc_soid_arch2ds(xls_name,xls_sheet,opts,preread) creates a matlab structure from 
% an xls archive of psychophysical data
%
%  xls_name: name of an xls file
%  xls_sheet: name of a sheet in the xls file
%  opts: options
%     opts.nplanes: number of planes to read
%     opts.fit_type is a prefix for the description
%     opts.ifaug is 1 to also read in values
%        of the exponent (b), fraction correct,  number of trials, number of amplitudes
%     opts.dataul_row controls which kind of fit to read
%     opts.limits_X: limit values for parameters
% preread: if present, preread.x* gives values of xnum, xtxt, xraw from the sheet that has already been read
%
%  ds:  a structure of one or more edirs, to be saved into a matlab file readable by btc_soid_getdata.
%  it has the following fields:
%  ds.subject
%  ds.desc
%  ds.condition{*}, a cell array for each background (or other) condition
%  ds.condition{*}.figgr (0, 1, or 2)
%  ds.condition{*}.desc
%  ds.condition{*}.edirs (has edirs fields: ds.condition{*}.edirs.yx.thresh_mags, returned as ed)
%
%  ou: options used
%
%  pr: values obtained by reading, can be used as preread next time through
%
%  See also:  BTC_SOID_GETDATA, BTC_SOID_XLSREAD, BTC_SOID_LIMITS, BTC_SOID_XLSFMT.
%
if (nargin<=2) opts=[]; end
opts=filldefault(opts,'fit_type','AllRaysFixedB'); %will go into the descr
opts=filldefault(opts,'ifaug',1);
opts=filldefault(opts,'noread',0); %set to 1 to skip reading
opts=btc_soid_xlsfmt(opts);
opts=btc_soid_limits(opts);
ou=opts;
if (nargin<=3) preread=[]; end
%
ds=[];
rowoff=[[0:opts.nplanes-1]'*opts.plane_spacing,zeros(opts.nplanes,1)]; %row offsets
xls_rc_plane=rowoff+repmat([opts.name_row,opts.name_col],opts.nplanes,1);
ou.xls_rc_plane=xls_rc_plane;
xls_rc_ndirs=rowoff+repmat([opts.ndir_row,opts.ndir_col],opts.nplanes,1);
ou.xls_rc_ndirs=xls_rc_ndirs;
xls_rc_dataul=rowoff+repmat([opts.dataul_row,opts.dataul_col],opts.nplanes,1);
ou.xls_rc_dataul=xls_rc_dataul;
xls_rc_figgr_spacing=repmat([opts.figgr_spacing 0],opts.nplanes,1);
ou.xls_rc_figgr_spacing=xls_rc_figgr_spacing;
%
if (opts.noread==1)
    return
end
%
ds.subject=xls_sheet;
ds.desc=sprintf('%s from %s',opts.fit_type,xls_name);
pr=preread;
for figgr=0:2
    c=[];
    c.figgr=figgr;
    if (figgr==0) c.desc='bkgd structured'; end
    if (figgr==1) c.desc='bkgd random'; end
    if (figgr==2) c.desc='bkgd combined'; end
    disp(c.desc);
    [c.edirs,ou,pr]=btc_soid_xlsread(xls_name,xls_sheet,xls_rc_plane,xls_rc_ndirs,...
        xls_rc_dataul+figgr*xls_rc_figgr_spacing,opts,pr);
    ds.condition{figgr+1}=c;
end
%
return


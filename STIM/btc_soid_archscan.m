function [scanned,expected,ou]=btc_soid_archscan(xls_name,xls_sheet,opts)
% [scanned,expected,ou]=btc_soid_archscan(xls_name,xls_sheet,opts) scans an xls archive of psychophysical data
%  to determine whether the basic formatting is correct
%
%  xls_name: name of an xls file
%  xls_sheet: name of a sheet in the xls file
%  opts: options (see btc_soid_limits, btc_soid_xlsfmt)
%
%  scanned. xnum, xtxt, xraw:  numeric, text, and raw contents of the spreadsheet
%  scanned.name_rows: rows that may have a plane name
%  scanned.ndir_rows: rows that may have an ndir value
%  expected.*: expected values of above params
%  ou: options used
%
%  See also:  BTC_SOID_ARCH2DS, BTC_SOID_XLSREAD, BTC_SOID_LIMITS, BTC_SOID_XLSFMT.
%
if (nargin<=2) opts=[]; end
opts=filldefault(opts,'fit_type','AllRaysFixedB'); %will go into the descr
opts=btc_soid_xlsfmt(opts);
opts=btc_soid_limits(opts);
ou=opts;
%
[xnum,xtxt,xraw]=xlsread(xls_name,xls_sheet);
scanned.xnum=xnum;
scanned.xtxt=xtxt;
scanned.xraw=xraw;
disp(sprintf('reading spreadsheet %s, sheet %s, for all planes of data',xls_name,xls_sheet));
disp(sprintf('size of read-in region is %4.0f x %4.0f',size(xraw)));
%
scanned.nplanes=floor(size(xraw,1)/opts.plane_spacing);
expected.nplanes=opts.nplanes;
%
% find the rows at the start of each plane; expected values calculated from
% plane spacing
%
name_rows=[];
ndir_rows=[];
for irow=1:size(xraw,1)
    if (ischar(xraw{irow,1}))
        name_rows=[name_rows,irow];
    end
    if (isnumeric(xraw{irow,1})) & ~(isnan(xraw{irow,1}))
        ndir_rows=[ndir_rows,irow];
    end
end
scanned.name_rows=name_rows;
scanned.ndir_rows=ndir_rows;
expected.name_rows=opts.name_row+[0:scanned.nplanes-1]*opts.plane_spacing;
expected.ndir_rows=opts.ndir_row+[0:scanned.nplanes-1]*opts.plane_spacing;
%
% find the rows at the start of a data block; expected values calculated
% from plane_spacing and first group
%
ray_rows=[];
for irow=1:size(xraw,1)
    qray=xraw{irow,opts.dataul_col-6};
    if (ischar(qray))
        if strcmp('ray',deblank(qray))
            ray_rows=[ray_rows,irow];
        end
    end
end
scanned.ray_rows=ray_rows;
ray_rows_plane1=ray_rows(find(ray_rows<=opts.plane_spacing));
expected.ray_rows=[];
for iplane=1:scanned.nplanes
    expected.ray_rows=[expected.ray_rows,opts.plane_spacing*(iplane-1)+ray_rows_plane1];
end
%
disp('comparing scanned and expected values')
%
scf=fieldnames(scanned);
for iscf=1:length(scf)
    fn=scf{iscf};
    if (isfield(expected,scf{iscf}))
        ifmatch=1;
        es=setdiff(expected.(fn),scanned.(fn));
        if ~isempty(es)
            disp(sprintf('parameter %20s has expected values that are not found in scan:',fn));
            disp(es);
            ifmatch=0;
        end
        se=setdiff(scanned.(fn),expected.(fn));
        if ~isempty(se)
            disp(sprintf('parameter %20s has scanned values that are not expected:',fn));
            disp(se);
            ifmatch=0;
        end
        if (ifmatch==1)
            disp(sprintf('parameter %20s has complete match of scanned and expected values.',fn));
        end
    end
end
return


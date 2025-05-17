% btc_soid_symtest: test symmetrizing of a btc psychophysical dataset
%
%   See also:  BTC_SOID_GETDATA, BTC_SOID_DEMO.
%
if ~exist('dict') dict=btc_define([]); end
nbtc=length(dict.codel); %10
%
% get the data from a named mat file
%
if ~exist('pdata_fn') pdata_fn='btc_allraysfixedb_mc'; end
pdata_fn=getinp('name of file with psychophysical data','s',[],pdata_fn);
%
[edata,eb_avail,ds,condno]=btc_soid_getdata(pdata_fn,setfield([],'if_choose',1));
%
planes=fieldnames(edata);
nplanes=length(planes);
%
edirs=[];
edirs_data=[];
ifok=1;
for iplane=1:nplanes
    pn=planes{iplane};
    edir=btc_edirs(pn);
    edirs=setfield(edirs,pn,edir);
    edirs_data=setfield(edirs_data,pn,edir);
    if isfield(edata,pn)
        if size(edata.(pn).thresh_mags,1)==edir.ndirs
            edirs_data.(pn).thresh_mags=edata.(pn).thresh_mags;
            edirs_data.(pn).thresh_vecs=edirs.(pn).uvecs.*repmat(edata.(pn).thresh_mags,1,2);
            if (eb_avail>0)
                if size(edata.(pn).thresh_mags_eblo,1)==edir.ndirs & size(edata.(pn).thresh_mags_eblo,1)==edir.ndirs
                    edirs_data.(pn).thresh_mags_eblo=edata.(pn).thresh_mags_eblo;
                    edirs_data.(pn).thresh_mags_ebhi=edata.(pn).thresh_mags_ebhi;
                    edirs_data.(pn).thresh_vecs_eblo=edirs.(pn).uvecs.*repmat(edata.(pn).thresh_mags_eblo,1,2);
                    edirs_data.(pn).thresh_vecs_ebhi=edirs.(pn).uvecs.*repmat(edata.(pn).thresh_mags_ebhi,1,2);
                else
                    disp(sprintf('Error bars for plane %2s from dataset in %s have wrong number of directions.',...
                        pn,pdata_fn))
                    ifok=0;
                end
            end
        else
            disp(sprintf('Threshold data for plane %2s from dataset in %s have wrong number of directions.',...
                pn,pdata_fn))
            ifok=0;
        end
    else
        disp(sprintf('Threshold data for plane %2s missing from dataset in %s.',...
            pn,pdata_fn))
        ifok=0;
    end
end
%
%
if ~exist('if_symm') if_symm=0; end
if_symm=getinp('1 to symmetrize the measured thresholds around the origin','d',[0 1],if_symm);
if ~exist('if_axes') if_axes=0; end
if_axes=getinp('1 to uniformize along axes','d',[-1 1],if_axes);
%
symopts=[];
symopts.if_symm=if_symm;
symopts.if_axes=if_axes;
[edirs_sym,symopts_used]=btc_soid_sym(edirs_data,symopts);
%
figure;
opts_plot=[];
opts_plot.datafield=strvcat('thresh_vecs','thresh_vecs_eblo','thresh_vecs_ebhi');
opts_plot.marker=strvcat('k','b','b');
opts_plot.cyclic=1;
set(gcf,'Position',[50 50 1200 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name','thresholds');
btc_soid_plot(edirs_data,opts_plot);
%
figure;
opts_plot=[];
opts_plot.datafield=strvcat('thresh_vecs','thresh_vecs_eblo','thresh_vecs_ebhi');
opts_plot.marker=strvcat('m','r','r');
opts_plot.cyclic=1;
set(gcf,'Position',[50 50 1200 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name','thresholds');
btc_soid_plot(edirs_sym,opts_plot);

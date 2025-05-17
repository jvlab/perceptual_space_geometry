% btc_soid_stats:  analyze the best-fitting ellipsoid to psychophysical discrimination data
%
% derived from btc_soid_demo.m, April 2014
%
%changes:
%  * defaults set for standard analysis (10 coords, madj=1, model variants just with and without symmetrizing
%  * plots by default are made but not saved
%  * further stats added
%  * some fields in r{:} output have slightly different names for uniformity
%  * display of number of unique points corrected for whether symmetry in positive and negative directions is assumed
%
%  modifications Nov 30 2011 to step through multiple fitting options
%  and save a cell array r{k}, with r{k}.results and other fields
%
%  modifications Nov 15 2012 to allow for fitting of models even if some thresholds are large
% (lines here and in btc_soid_fit involving use_large)
% 
%  modifications Nov 24 2013 to invoke 'madj' option for btc_edirs, to adjust diagonal directions
%
%   See also:  BTC_DEFINE, BTC_PAIRSNEEDED, BTC_SOID_PLOT, BTC_AUGCOORDS,
%   BTC_SOID_FIND, BTC_SOID_PLOT3D, BTC_PRED_DEMO, BTC_SOID_DEMO.
%
%
if ~exist('dict') dict=btc_define([]); end
nbtc=length(dict.codel); %10
if ~exist('coordsets'); coordsets={'atg','atbg','bcdetuvw','gbcdetuvwa','abcg'}; end
if ~exist('triplets');
    triplets{1}={'atg'};
    triplets{2}={'abg','atg','abt'};
    triplets{3}={'bcd','bde','tuv','bct'};
    triplets{4}={'atg','abg','adg','abd','atu','atb','tuv','tub','bcd','bde','bct','bdt'};
    triplets{5}={'abg','acg','abc','bcg'};
end
if ~exist('csno'); csno=4; end
if ~exist('madj') madj=1; end
madj=getinp('1 to adjust diagonals to match BigConds ("madj")','d',[0 1],madj);
for ics=1:length(coordsets)
    disp(sprintf(' coord set %2.0f->%10s',ics,coordsets{ics}));
end
csno=getinp('choice','d',[1 length(coordsets)],csno);
coords=coordsets{csno};
ncoords=length(coords);
if ~exist('symfull') symfull=1; end
symfull=getinp('1 for symmetric setup, 2 for full setup, 3 for symmetric setup without subs','d',[1 3],1);
%
[planes_sym,planes_full,planes_sym_unsub]=btc_pairsneeded(coords,dict);
if (symfull==1) planes=planes_sym; end
if (symfull==2) planes=planes_full; end
if (symfull==3) planes=planes_sym_unsub; end
disp(sprintf(' %2.0f planes needed with symmetry (%2.0f for full set)',...
    size(planes_sym,1),size(planes_full,1)));
nplanes=size(planes,1);
%
% get the data from a named mat file
%
if ~exist('pdata_fn') pdata_fn='btc_soid_demo'; end
pdata_fn=getinp('name of file with psychophysical data (e.g., btc_allraysfixedb_xx)','s',[],pdata_fn);
%
[edata,eb_avail,ds,condno]=btc_soid_getdata(pdata_fn,setfield([],'if_choose',1));
%
%set up the experiment directions
%transfer the experimental data
%and check that there are the right number of directions
%
edirs=[];
edirs_data=[];
ifok=1;
for iplane=1:nplanes
    pn=planes(iplane,:);
    edir=btc_edirs(pn,setfield([],'madj',madj));
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
if (ifok==1)
    if (eb_avail>0)
        disp(sprintf('Threshold data and error bars successfully transferred for %2.0f planes.',nplanes))
    else
        disp(sprintf('Threshold data without error bars successfully transferred for %2.0f planes.',nplanes))
    end
else
    disp('Cannot proceed.');
end
%
if ~exist('symopts') symopts=[]; end
symopts=filldefault(symopts,'if_symm',0);
symopts=filldefault(symopts,'if_axes',0);
%
if ~exist('setups') setups=[]; end
setups=filldefault(setups,'if_symm_list',[0 1]); %basic: use [0 1]
setups=filldefault(setups,'if_axes_list',[-1 0 1]); %basic: use [0 1]
setups=filldefault(setups,'ifaug_list',[0 1]); %basic: use 1
setups.if_symm_list=getinp('list of (0 for no symmetrization around origin, 1 to symmetrize)','d',[0 1],setups.if_symm_list);
setups.if_axes_list=getinp('list of (1 to uniformize the thresholds on axes, (0 for not, -1 to keep separate within [bc],[de],[tuvw])','d',[-1 1],1);
setups.ifaug_list=getinp('list of (1 to augment coordinates to match stimuli [default], 0 to ignore augmentation)','d',[0 1],1);
%
setups.variants=cell(0);
nvariants=0;
for ptr_symm=1:length(setups.if_symm_list)
    for ptr_axes=1:length(setups.if_axes_list)
        for ptr_aug=1:length(setups.ifaug_list)
            nvariants=nvariants+1;
            setups.variants{nvariants}.if_symm=setups.if_symm_list(ptr_symm);
            setups.variants{nvariants}.if_axes=setups.if_axes_list(ptr_axes);
            setups.variants{nvariants}.ifaug=setups.ifaug_list(ptr_aug);
            setups.variants{nvariants}.label=sprintf('sy=%1.0f ax=%1.0f au=%1.0f',...
                setups.variants{nvariants}.if_symm,...
                setups.variants{nvariants}.if_axes,...
                setups.variants{nvariants}.ifaug);
            disp(sprintf('setup %2.0f->%s',nvariants,setups.variants{nvariants}.label))
        end
    end
end
if ~exist('nsurr') nsurr=10; end
if (eb_avail>0)
    nsurr=getinp('number of surrogates to run','d',[0 10000],nsurr);
else
    disp('Number of surrogates set to 0 as error bars are not present.')
    nsurr=0;
end
use_large=getinp('1 to use out-of-range values for fitting','d',[0 1],0);
%
if_savefigs=getinp('1 to save figs as fig file','d',[0 1],0);
if_psfigs=getinp('1 to save figs as ps','d',[0 1],0);
if_closefigs=getinp('1 to close figs after making','d',[0 1],0);
if ~exist('fig_infix')
    if (madj==0)
        fig_infix=sprintf('_%1.0fsurrs_stats',nsurr);
    else
        fig_infix=sprintf('_%1.0fsurrs_stats_madj',nsurr);
    end
end

%if (if_savefigs>0) | (if_psfigs>0)
if (0==0)
    fig_infix=getinp('figure and ps infix and suggested saved name, with leading underscore','s',[],fig_infix);
    if ~(fig_infix(1)=='_');
        fig_infix=cat(2,'_',fig_infix);
    end
    if (strcmp(fig_infix,'_'))
        fig_infix='';
    end
end
r=cell(0);
%now step through each variant
for ivariant=1:nvariants
    fig_basename=cat(2,strrep(pdata_fn,'.mat',''),fig_infix,'_var',zpad(ivariant,3));
    r{ivariant}.setup=setups.variants{ivariant};
    vlabel=setups.variants{ivariant}.label;
    disp(' ');
    disp(sprintf('Quadratic form analysis for file %s',pdata_fn'));
    disp(sprintf('setup %2.0f->%s',ivariant,setups.variants{ivariant}.label))
    symopts.if_symm=setups.variants{ivariant}.if_symm;
    symopts.if_axes=setups.variants{ivariant}.if_axes;
    ifaug=setups.variants{ivariant}.ifaug;
    %
    edirs_data_orig=edirs_data;
    [edirs_data,symopts_used]=btc_soid_sym(edirs_data_orig,symopts);
    unique_onaxis=symopts_used.onaxis_list(symopts_used.onaxis_coord_table_ptrs);
    unique_offaxis=setdiff([1:size(symopts_used.all_coord_table,1)],symopts_used.onaxis_list);
    %
    disp(rmfield(ds,'condition'));
    disp(sprintf('Condition %2.0f',condno));
    disp(rmfield(ds.condition{condno},'edirs'));
    disp(' ');
    %
    opts_fit=[];
    opts_fit.verbose=0;
    opts_fit.coords=coords;
    opts_fit.symfull=symfull;
    opts_fit.nsurr=nsurr;
    opts_fit.ifaug=ifaug;
    opts_fit.use_large=use_large;
    disp('fitting options');
    disp(opts_fit);
    disp('running btc_soid_fit');
    [results,ou_fit,ou_dict]=btc_soid_fit(edirs,edirs_data,opts_fit,dict);
    r{ivariant}.results=results;
    r{ivariant}.ou_fit=ou_fit;
    r{ivariant}.ou_dict=ou_dict;
    r{ivariant}.edirs=edirs;
    %
    qform_fit=results.qfit; %for later plotting
    qfa=[];
    qfa(:,:,1)=qform_fit;
    %
    % augment qform_fit to qform_fit_aug with 0's in missing locations
    %
    qform_fit_aug=zeros(length(dict.codel));
    qform_fit_aug(results.which_axes,results.which_axes)=qform_fit;
    if (nsurr>0)
        qform_med=median(results.surrogates.qfit,3);
        qform_med_aug=zeros(length(dict.codel));
        qform_med_aug(results.which_axes,results.which_axes)=qform_med;
    end
    %
    disp('symmetrizing options:');
    disp(symopts);
    %
    %here, calculate differences between fitted and measured thresholds as in-plane magnitudes
    %and the prediction of qform_fit
    %
    allmag=[];
    %
    disp('using fitted qform')
    [edirs_fits,allmag,allmag_fits,allmag_std]=btc_soid_predict(edirs,edirs_data,qform_fit_aug,dict);
    nonan=find(~isnan(allmag));
    % axes and diagonal pointers already determined by btc_soid_sym
    nonan_onaxis=intersect(nonan,unique_onaxis);
    nonan_offaxis=intersect(nonan,unique_offaxis);
    %
    nthresh=length(nonan);
    disp(sprintf('There are %3.0f threshold values (%3.0f are out of range)',nthresh,length(allmag)-nthresh));
    dmf=allmag-allmag_fits;
    r{ivariant}.mean_data_minus_fit=mean(dmf(nonan));
    r{ivariant}.mean_data_minus_fit_onaxis=mean(dmf(nonan_onaxis));
    r{ivariant}.mean_data_minus_fit_offaxis=mean(dmf(nonan_offaxis));
    r{ivariant}.rmse=sqrt(mean(dmf(nonan).^2));
    r{ivariant}.rmse_onaxis=sqrt(mean(dmf(nonan_onaxis).^2));
    r{ivariant}.rmse_offaxis=sqrt(mean(dmf(nonan_offaxis).^2));
    r{ivariant}.rmsstd=sqrt(mean(allmag_std(nonan).^2));
    r{ivariant}.rmsstd_onaxis=sqrt(mean(allmag_std(nonan_onaxis).^2));
    r{ivariant}.rmsstd_offaxis=sqrt(mean(allmag_std(nonan_offaxis).^2));
    r{ivariant}.rmsz=sqrt(mean((dmf(nonan)./allmag_std(nonan)).^2));
    r{ivariant}.rmsz_onaxis=sqrt(mean((dmf(nonan_onaxis)./allmag_std(nonan_onaxis)).^2));
    r{ivariant}.rmsz_offaxis=sqrt(mean((dmf(nonan_offaxis)./allmag_std(nonan_offaxis)).^2));
    disp('unique values        mean data-fit   rms(data-fit) rms(data std)  rms((data-fit)/data_std)');
    disp(sprintf(' on-axis:  %5.0f       %7.4f        %7.4f       %7.4f          %7.4f',...
        length(nonan_onaxis)/(1+symopts.if_symm),...
        r{ivariant}.mean_data_minus_fit_onaxis,r{ivariant}.rmse_onaxis,r{ivariant}.rmsstd_onaxis,r{ivariant}.rmsz_onaxis));
    disp(sprintf('off-axis:  %5.0f       %7.4f        %7.4f       %7.4f          %7.4f',...
        length(nonan_offaxis)/(1+symopts.if_symm),...
        r{ivariant}.mean_data_minus_fit_offaxis,r{ivariant}.rmse_offaxis,r{ivariant}.rmsstd_offaxis,r{ivariant}.rmsz_offaxis));
    disp(sprintf(' all pts:  %5.0f       %7.4f        %7.4f       %7.4f          %7.4f',...
        length(nonan)/(1+symopts.if_symm),...
        r{ivariant}.mean_data_minus_fit,r{ivariant}.rmse,r{ivariant}.rmsstd,r{ivariant}.rmsz));
%
    if (nsurr>0)
        disp('using median of surrogate-fitted qform')
        [edirs_fits_med,allmag_med,allmag_fits_med]=...
            btc_soid_predict(edirs,edirs_data,qform_med_aug,dict); %allmag_med should match allmag
        disp(sprintf('rmse= %7.4f',sqrt(mean((allmag_med(nonan)-allmag_fits_med(nonan)).^2))));
        r{ivariant}.rmse_med_surrogates=sqrt(mean((allmag_med(nonan)-allmag_fits_med(nonan)).^2));
        qfa(:,:,[2:4])=quantile(results.surrogates.qfit,[.025 .975 0.5],3);
    end
    % merge the data and fitted thresholds, and plot
    edirs_data_and_fits=edirs_data;
    for iplane=1:nplanes
        pn=planes(iplane,:);
        edirs_data_and_fits.(pn).thresh_vecs_fits=edirs_data.(pn).uvecs.*repmat(edirs_fits.(pn).thresh_mags,1,2);   
    end
    r{ivariant}.edirs_data_and_fits=edirs_data_and_fits;
    %
    figure;
    set(gcf,'Position',[50 50 1200 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',cat(2,'thresholds ',vlabel));
    opts_plot=[];
    if (nsurr>0)
        opts_plot.datafield=strvcat('thresh_vecs','thresh_vecs_eblo','thresh_vecs_ebhi','thresh_vecs_fits');
        opts_plot.marker=strvcat('k','b','b','r');
    else
        opts_plot.datafield=strvcat('thresh_vecs','thresh_vecs_fits');
        opts_plot.marker=strvcat('k','r');
    end
    opts_plot.cyclic=1;
    opts_plot.tstring=vlabel;
    btc_soid_plot(edirs_data_and_fits,opts_plot);
    if (if_psfigs)
        orient landscape;
        print(gcf,'-dwinc','-append',cat(2,fig_basename,'.ps'));
    end
    if (if_savefigs)
        fig_savename=cat(2,fig_basename,'_thresh');
        saveas(gcf,fig_savename);
    end
    if (if_closefigs)
        close(gcf);
    end
    %
    %plot a scattergram of thresholds and fitted thresholds
    %
    figure;
    set(gcf,'Position',[50 50 1200 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',cat(2,'error summary ',vlabel));
    ilmax=max(allmag(nonan))*1.2;
    for ilinlog=0:1
        for iscat=1:1+2*double(nsurr>0)
            switch iscat
                case 1
                    xvals=allmag;
                    xlabel_string='measured';
                    yvals=allmag_fits;
                    ylabel_string='fitted';
                case 2
                    xvals=allmag;
                    xlabel_string='measured';
                    yvals=allmag_fits_med;
                    ylabel_string='fitted (median surr)';
                case 3
                    xvals=allmag_fits;
                    xlabel_string='fitted';
                    yvals=allmag_fits_med;
                    ylabel_string='fitted (median surr)';
            end
            subplot(2,3,iscat+ilinlog*3);
            plot(xvals(nonan_onaxis),yvals(nonan_onaxis),'k.');hold on;
            plot(xvals(nonan_offaxis),yvals(nonan_offaxis),'m.');hold on;
            title(vlabel);
            xlabel(xlabel_string);
            ylabel(ylabel_string);
            if (ilinlog==0)
                ilmin=0;
            end
            if (ilinlog==1)
                ilmin=0.1;
                set(gca,'XScale','log');
                set(gca,'YScale','log');
            end
            plot([ilmin ilmax],[ilmin ilmax],'k--');hold on;
            set(gca,'XLim',[ilmin ilmax]);
            set(gca,'YLim',[ilmin ilmax]);
            axis equal;
            axis tight;
            legend('axes','diags','Location','SouthEast');    
        end %iscat
    end %ilinlog
    if (if_psfigs)
        orient landscape;
        print(gcf,'-dwinc','-append',cat(2,fig_basename,'.ps'));
    end
    if (if_savefigs)
        fig_savename=cat(2,fig_basename,'_errsum');
        saveas(gcf,fig_savename);
    end
    if (if_closefigs)
        close(gcf);
    end
    %
    %plot quadratic form values and matrix square roots
    %
    figure;
    set(gcf,'Position',[50 50 1200 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',cat(2,'qf values ',vlabel));
    opts_plotqf=[];
    opts_plotqf.coords=coords;
    opts_plotqf.view=[-60 50];
    subplot(2,2,1);
    btc_soid_plotqf(qform_fit,opts_plotqf);
    title(cat(2,'qform ',vlabel));
    subplot(2,2,3);
    btc_soid_plotqf(real(sqrtm(qform_fit)),opts_plotqf);
    title(cat(2,'qform matrix sqrt ',vlabel));
    results.qfit_sqrtm=sqrtm(qform_fit);
    if (nsurr>0)
        subplot(2,2,2);
        btc_soid_plotqf(qfa,opts_plotqf);
        title(cat(2,'qform from surrogates ',vlabel));
        %matrix square roots from surrogates   
        msqrts=zeros(size(results.surrogates.qfit));
        for isurr=1:nsurr
            msqrts(:,:,isurr)=sqrtm(results.surrogates.qfit(:,:,isurr));
        end
        subplot(2,2,4);
        btc_soid_plotqf(cat(3,real(sqrtm(qform_fit)),quantile(real(msqrts),[.025 .975 0.5],3)),opts_plotqf);
        title(cat(2,'qform matrix sqrt from surrogates ',vlabel));
        results.surrogates.qfit_sqrtm=msqrts;
    end
    if (if_psfigs)
        orient landscape;
        print(gcf,'-dwinc','-append',cat(2,fig_basename,'.ps'));
    end
    if (if_savefigs)
        fig_savename=cat(2,fig_basename,'_qfvals');
        saveas(gcf,fig_savename);
    end
    if (if_closefigs)
        close(gcf);
    end
    %
    %plot quadratic form (and confidence limits) as colormaps
    %
    figure;
    set(gcf,'Position',[100 50 1000 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',cat(2,'qf colormap ',vlabel));
    zscale=[-1 1]*max(qfa(:));
    for iqf=1:size(qfa,3)
        if (iqf==1) tstring='qform fitted'; sp=1; end
        if (iqf==2) tstring='qform lower CL'; sp=3; end
        if (iqf==3) tstring='qform upper CL'; sp=4; end
        if (iqf==4) tstring='qform median'; sp=2; end
        tstring=cat(2,tstring,' ',vlabel);
        subplot(2,2,sp);
        set(gca,'Fontsize',7);
        imagesc(qfa(:,:,iqf),zscale);
        set(gca,'XLim',0.5+[0 length(coords)]);
        set(gca,'YLim',0.5+[0 length(coords)]);
        set(gca,'XTick',[1:length(coords)]);
        set(gca,'XTickLabel',coords');
        set(gca,'XAxisLocation','top');
        set(gca,'YTick',[1:length(coords)]);
        set(gca,'YTickLabel',coords');
        title(tstring);
        axis square;
        colorbar;
    end
    if (if_psfigs)
        orient landscape;
        print(gcf,'-dwinc','-append',cat(2,fig_basename,'.ps'));
    end
    if (if_savefigs)
        fig_savename=cat(2,fig_basename,'_qfmap');
        saveas(gcf,fig_savename);
    end
    if (if_closefigs)
        close(gcf);
    end
    %
    % plot triplets
    %
    if ~isempty(triplets{csno})
        figure;
        set(gcf,'Position',[100 50 1000 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,'triplets ',vlabel));
        if ismember(symfull,[1 3])
            opts_plot.use_sym=1;
        else
            opts_plot.use_sym=0;
        end
        btc_soid_plot3d(edirs_data_and_fits,triplets{csno},setfields(opts_plot,{'dperm','tstring'},{[3 2 1],vlabel}));
    end
    if (if_psfigs)
        orient landscape;
        print(gcf,'-dwinc','-append',cat(2,fig_basename,'.ps'));
    end
    if (if_savefigs)
        fig_savename=cat(2,fig_basename,'_triplets');
        saveas(gcf,fig_savename);
    end
    if (if_closefigs)
        close(gcf);
    end
end %ivariant
disp(' ');
disp(' Summary');
for ivariant=1:nvariants
    disp(sprintf('setup %2.0f: %s, rmse=%7.5f',ivariant,r{ivariant}.setup.label,r{ivariant}.rmse));
end
savefile=cat(2,strrep(pdata_fn,'.mat',''),fig_infix);
if getinp(sprintf('1 to save r in %s',savefile),'d',[0 1],1);
    save(savefile,'r');
end
%
clear ivariant vlabel
clear ifn fn fns
clear iplane ics icond ix spec specx vecs_inplane
clear ou iqf qf tstring sp

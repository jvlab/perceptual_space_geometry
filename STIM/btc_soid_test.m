% btc_soid_test:  test routines for finding the best-fitting ellipsoid
% to psychophysical discrimination data; for use see btc_soid_demo
%
%   See also:  BTC_DEFINE, BTC_PAIRSNEEDED, BTC_SOID_PLOT, BTC_AUGCOORDS, BTC_SOID_FIND.
%
if ~exist('dict') dict=btc_define([]); end
nbtc=length(dict.codel); %10
if ~exist('coordsets'); coordsets={'atg','atbg','bcdetuvw','gbcdetuvwa'}; end
if ~exist('csno'); csno=1; end
if ~exist('qno') qno=1; end
if ~exist('qparams') qparams=[0.5 1]; end
if ~exist('symfull') symfull=1; end
%
if ~exist('dthr') dthr=0.5; end %dthr=distance from origin used to determine thresholds
%
for ics=1:length(coordsets)
    disp(sprintf(' coord set %2.0f->%10s',ics,coordsets{ics}));
end
csno=getinp('choice','d',[1 length(coordsets)],csno);
coords=coordsets{csno};
ncoords=length(coords);
%
disp('quadratic form options for the just noticeable difference')
disp('0->already created as qform')
disp('1->identity');
disp('2->ramp of values on diagonal');
disp('3->one off-diagonal value, one on-diagonal value');
disp('4->hand-tuned set with some off-diagonal values');
qno=getinp('choice','d',[0 4],qno);
if (qno>0)
    dthr=getinp('distance required to reach threshold','f',[0 1],dthr);
end
if (qno==1)
    qform=eye(nbtc);
    qform=qform/dthr^2;
end
if (qno==2)
    qparams=getinp('low and high diagonal value','f',[0 Inf],qparams);
    qform=diag(qparams(1)+(qparams(2)-qparams(1))*[0:(nbtc-1)]/(nbtc-1));
    qform=qform/dthr^2;
end
if (qno==3)
    qparams=getinp('off and on diagonal value','f',[0 Inf],qparams);
    qform=qparams(1)*ones(nbtc)+(qparams(2)-qparams(1))*eye(nbtc);
    qform=qform/dthr^2;
end
if (qno==4)
    qform=diag([10 2.5 2.5 1.5 1.5 0.4 0.4 0.4 0.4 0.625]);
    qform(2,3)=1; %beta, beta interaction
    qform(3,2)=1;
    qform(2:3,10)=1.3; %beta, alpha interaction
    qform(3:2,10)=1.3;
    qform(1,[6:9])=0.5; %gamma, theta interacton
    qform([6:9],1)=0.5;
    qform=qform/dthr^2;
end
qeigs=eig(qform);
disp(sprintf(' eigenvalues range from %7.4f to %7.4f',min(qeigs),max(qeigs)));
if (min(qeigs)<=0)
    disp('WARNING: the quadratic form is not positive definite');
end
%
symfull=getinp('1 for symmetric setup, 2 for full setup, 3 for symmetric setup without subs','d',[1 3],symfull);
%
[planes_sym,planes_full,planes_sym_unsub]=btc_pairsneeded(coords,dict);
if (symfull==1) planes=planes_sym; end
if (symfull==2) planes=planes_full; end
if (symfull==3) planes=planes_sym_unsub; end
disp(sprintf(' %2.0f planes needed with symmetry (%2.0f for full set)',...
    size(planes_sym,1),size(planes_full,1)));
nplanes=size(planes,1);
%
if ~exist('noise_frac') noise_frac=.25; end
if ~exist('noise_ebamt') noise_ebamt=.2; end
if ~exist('if_fit') if_fit=1; end
if ~exist('nsurr') nsurr=10; end
noise_frac=getinp('noise fraction for simulated data','f',[0 0.5],noise_frac);
noise_ebamt=getinp('error bar fractional size for simulated data','f',[0 1],noise_ebamt);
if_fit=getinp('1 to run btc_soid_fit','d',[0 1],if_fit);
if (if_fit==1)
    nsurr=getinp('number of surrogates to run','d',[0 1000],nsurr);
end
%
%now set up the experiment directions
%
edirs=[];
for iplane=1:nplanes
    edirs=setfield(edirs,planes(iplane,:),btc_edirs(planes(iplane,:)));
end
%plot the values, just to verify
%
figure;
set(gcf,'Position',[50 50 1200 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name','experiment setup');
hplanes=btc_soid_plot(edirs,setfield([],'datafield','allvecs_concat'));
%
% create concatenated arrays of coordinates with and without augmentation
% at a fixed (small) multiple of unity -- these will be used to establish
% thresholds for a noisy dataset
%
vecs_thr_aug=[];
edirs_thr=[]; %an edirs-structure with exact thresholds corresponding to qform=1
%
%determine the experimental design and create a structure
%with threshold data
for iplane=1:nplanes
   disp(sprintf('plane %s: creating ideal dataset',planes(iplane,:)))
   edir=getfield(edirs,planes(iplane,:));
   %set up a subfield to hold the thresholds
   edir_thr=[];
   edir_thr.plane=edir.plane;
   edir_thr.opts=edir.opts;
   edir_thr.thresh_vecs=[];
   edir_thr.maxvecs=edir.maxvecs; %to plot also
   edir_thr.thresh_mags=[];
   %
   vecs_inplane=getfield(edirs,planes(iplane,:),'uvecs');
   spec=[];
   for icond=1:edir.ndirs
       for ix=1:2
            %the "3-ix" is so that the SECOND coordinate of vecs_inplane 
            %is the FIRST btc coordinate (in planes)
            spec=setfield(spec,planes(iplane,ix),vecs_inplane(icond,3-ix));
       end
       [specx,avec,results,ou]=btc_soid_find(spec,1,qform,dict);
       % Here is where we make the second coordinate of plane into the first
       % coordinate that is plotted
       thr_inplane=[getfield(specx,planes(iplane,2)),getfield(specx,planes(iplane,1))];
       edir_thr.thresh_vecs=[edir_thr.thresh_vecs;thr_inplane];
       vecs_thr_aug=[vecs_thr_aug;avec];
       edir_thr.thresh_mags(icond,1)=sqrt(sum(thr_inplane.^2));
       edirs_thr=setfield(edirs_thr,planes(iplane,:),edir_thr);
   end
end
% check whether the ideal thresholds are correct
distsq=vecs_thr_aug*qform*vecs_thr_aug';
disp(sprintf(' maximum deviation of squared distance from 1: %12.5f',max(abs(diag(distsq)-1))));
%
edirs_noise_thr=[]; %an edirs-structure with noisy thresholds
for iplane=1:nplanes
   disp(sprintf('plane %s: creating noisy dataset',planes(iplane,:)))
   edir=getfield(edirs,planes(iplane,:));
   edir_thr=getfield(edirs_thr,planes(iplane,:));
   noise_mult=max(0,1+noise_frac*randn(edir.ndirs,1));
   edir_noise_thr=[];
   edir_noise_thr.maxvecs=edir_thr.maxvecs;
   edir_noise_thr.plane=edir_thr.plane;
   edir_noise_thr.opts=edir_thr.opts;
   edir_noise_thr.thresh_mags     =noise_mult.*edir_thr.thresh_mags;
   edir_noise_thr.thresh_mags_eblo=noise_mult.*edir_thr.thresh_mags./(1+noise_ebamt);
   edir_noise_thr.thresh_mags_ebhi=noise_mult.*edir_thr.thresh_mags.*(1+noise_ebamt);
   edir_noise_thr.thresh_vecs     =repmat(edir_noise_thr.thresh_mags,1,2).*edir.uvecs;
   edir_noise_thr.thresh_vecs_eblo=repmat(edir_noise_thr.thresh_mags_eblo,1,2).*edir.uvecs;
   edir_noise_thr.thresh_vecs_ebhi=repmat(edir_noise_thr.thresh_mags_ebhi,1,2).*edir.uvecs;
   %
   edirs_noise_thr=setfield(edirs_noise_thr,planes(iplane,:),edir_noise_thr);
   %
end
%
%
disp(sprintf('edirs_thr or edirs_noise_thr, qform, and coords now set up to test ellipsoid-finding routines'));
%
figure;
set(gcf,'Position',[50 50 1200 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name','thresholds');
opts_thr=[];
opts_thr.datafield=strvcat('maxvecs','thresh_vecs');
opts_thr.marker=strvcat('k.','r-');
opts_thr.cyclic=[0 1];
opts_hplanes_thr=btc_soid_plot(edirs_thr,opts_thr);
%
figure;
set(gcf,'Position',[50 50 1200 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name','noisy thresholds');
opts_thr=[];
opts_thr.datafield=strvcat('maxvecs','thresh_vecs','thresh_vecs_eblo','thresh_vecs_ebhi');
opts_thr.marker=strvcat('k.','r-','g','b');
opts_thr.cyclic=[0 1];
opts_hplanes_thr=btc_soid_plot(edirs_noise_thr,opts_thr);
%
if (if_fit)
    disp('running btc_soid_fit');
    opts_fit=[];
    opts_fit.verbose=1;
    opts_fit.coords=coords;
    opts_fit.symfull=symfull;
    opts_fit.nsurr=nsurr;
    %
    % do the regression
    %
    [results,ou_fit,ou_dict]=btc_soid_fit(edirs,edirs_noise_thr,opts_fit,dict);
    qform_fit=results.qfit;
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
    %here, calculate differnces between predicted and measured thresholds as in-plane magnitudes
    %first argument  edirs.yx.uvecs
    %second argument edirs[_noise].yx.thresh_mags)
    %and the prediction of qform_fit
    disp('using ideal dataset and known qform')
    [edirs_ideal     ,allmag_ideal,allmag_ideal_pred]      =btc_soid_predict(edirs,edirs_thr,qform,dict);
    disp(sprintf('rmse= %7.4f',sqrt(mean((allmag_ideal-allmag_ideal_pred).^2))));
    disp('using ideal dataset and fitted qform')
    [edirs_pred      ,allmag      ,allmag_pred]      =btc_soid_predict(edirs,edirs_thr,qform_fit_aug,dict);
    disp(sprintf('rmse= %7.4f',sqrt(mean((allmag-allmag_pred).^2))));
    disp('using noisy dataset and fitted qform')
    [edirs_noise_pred,allmag_noise,allmag_noise_pred]=btc_soid_predict(edirs,edirs_noise_thr,qform_fit_aug,dict);
    disp(sprintf('rmse= %7.4f',sqrt(mean((allmag_noise-allmag_noise_pred).^2))));
    if (nsurr>0)
        disp('using noisy dataset and median of surrogate-fitted qform')
        [edirs_noise_pred_med,allmag_noise_med,allmag_noise_pred_med]=...
            btc_soid_predict(edirs,edirs_noise_thr,qform_med_aug,dict);
        disp(sprintf('rmse= %7.4f',sqrt(mean((allmag_noise_med-allmag_noise_pred_med).^2))));
    end
    %
    figure;
    set(gcf,'Position',[100 50 1000 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name','qform and fit');
    zscale=[-1 1]*max(max(abs(qform(:))),max(abs(qform_fit(:))));
    for iqf=1:4
        if (iqf==1) qf=qform(results.which_axes,results.which_axes); tstring='qform'; end
        if (iqf==2) qf=qform_fit; tstring='qform fitted'; end
        if (iqf==3) qf=(qform(results.which_axes,results.which_axes)-qform_fit); tstring='difference'; end
        if (iqf==4) qf=(qform(results.which_axes,results.which_axes)-qform_fit)*10; tstring='difference x 10'; end
        subplot(2,2,iqf);
        set(gca,'Fontsize',7);
        imagesc(qf,zscale);
        set(gca,'XLim',0.5+[0 length(coords)]);
        set(gca,'YLim',0.5+[0 length(coords)]);
        set(gca,'XTick',[1:length(coords)]);
        set(gca,'XTickLabel',coords');
        set(gca,'XAxisLocation','top');
        set(gca,'YTick',[1:length(coords)]);
        set(gca,'YTickLabel',coords');
        title(tstring);
        axis square;
    end
    %
    figure;
    set(gcf,'Position',[50 50 1200 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name','regression setup');
    fns=fieldnames(results.qsetup);
    [nr,nc]=nicesubp(length(fns),0.7);
    for ifn=1:length(fns)
        fn=fns{ifn};
        subplot(nr,nc,ifn,'align');
        set(gca,'Fontsize',7);
        qones=getfield(getfield(results.qsetup,fn),'qones');
        spy(qones); hold on;
        plot(0.5+[0 length(coords)],0.5+[0 length(coords)],'k');
        title(fn);
        set(gca,'XLim',0.5+[0 length(coords)]);
        set(gca,'XTick',[1:length(coords)]);
        set(gca,'XTickLabel',coords');
        set(gca,'XAxisLocation','top');
        set(gca,'YLim',0.5+[0 length(coords)]);
        set(gca,'YTick',[1:length(coords)]);
        set(gca,'YTickLabel',coords');
        xlabel('');
    end
    clear ifn fn fns qones
end
clear iplane ics icond ix spec specx vecs_inplane
clear edir edir_thr edir_noise_thr noise_mult
clear ou iqf qf tstring

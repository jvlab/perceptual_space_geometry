function [img,stats,optsused,errs,metro]=btc_nopick(opts,area)
%[img,stats,optsused,errs,metro]=btc_nopick(opts,area) creates a binary noise image
%  via Markov random field followed by Metropolis for 2x2 blocks with two 
%  adjacent nonzero theta parameters
% This handles the method NoPickTT and NoPickBT in btc_makemaps.  Many maps can be
% made at once, to facilitate Metropolis stage.  So the internal stats
% routine is different from that of genmrfm.
%
% Initial stage based on genmrfm and genmrfmg, final stage based on metro*.m
% Metropolis stage params passed as opts.metro_opts and opts.metro_show.
%
% The checks A, B, C, and D referred to by the recursion relationship are
% A B
% C D
%
% This is used row by row, with (in standard orientation) nonzero theta params for ACD and BCD.
%
% opts.recurtype: 'NoPickTT' or 'NoPickBT'
%  [for NoPickTT]
%   opts.mACD:  theta-ACD (when rotated into standard position)
%   opts.mBCD:  theta-BCD (when rotated into standard position)
%   opts.mABCD: alpha-ABCD 
%  [for NoPickBT]
%   opts.mABC:  theta-ABC (when rotated into standard position)
%   opts.mAD:   beta-AD   (when rotated into standard position)
%     T-rule is defined on
%     tA tB tC
%        tD
%   t(BCD) comes from theta-ABC; t(AD) comes from beta-AD, and t(ABC)=t(AD)^2*t(BCD) must be added
%   to keep stationarity
% opts.show: set to 1 to draw the image
% opts.statsub:  number of statistical subdivisions on each axis (defaults to 2, must be >0)
% opts.err_rept: 0 to abort on errors, 1 to report them into errs
% opts.rotvar:  rotation variant. 1->ACD/BCD, 2->BCD/ABD, 3->ABD/ABC, 4->ABC/ACD.
% opts.nmaps: number of maps to make (defaults to 1)
% opts.burnin:  number of "burn-in" rows
% opts.minarea: (for future) minimum area to make, to allow for adequate Metropolis scrambling
% opts.metro_opts:  Metropolis options
% opts.metro_show: set to show Metropolis details
% opts.verbose: 1 for a log, 2 for a detailed log
% opts.stability_tol: tolerance for verifying stability of 2x2 rule
% opts.onaxis: 1 if on axis (i.e., at most one param value is nonzero); this bypasses Metropolis
% opts.method:  method, as determined by btc_augcoords or mtc_augcoords
% opts.mtcopts: options structure from mtc_defopts
%
% area: size of img generated, defaults to [256 256]
%
% uses recursive nonvectorized algorithm, and is therefore slow
%
% img: the image; size(img) = area
% stats: some image statistics (luminance bias and eo bias, for consistency
%   with genmrfm and btc_dimrf)
% optsused: the options used
% errs:  string to report errors
% metro: results of call to donut_metro
%
% 5 Oct 2015:  but fixed that resulted in Metropolis outputs NOT being
%    used, see ptbxcvt_btc_notes.txt
%
%   See also GENMRFM, BTC_AUGCOORDS, BTC_MAKEMAPS, NICESUBP, BTC_3X3ANAL, BTC_TRULE_ITERATE, GENMRFM_TEE.
%
metro=[];
ard=[256 256];
if (nargin >= 2)
   if (length(area) >= 2)
      ard=[area(1) area(2)];
   else
      ard=[area(1) area(1)];
   end
end
opts=filldefault(opts,'show',0);
opts=filldefault(opts,'statsub',2);
opts=filldefault(opts,'err_rept',0);
opts=filldefault(opts,'rotvar',1);
opts=filldefault(opts,'nmaps',1);
opts=filldefault(opts,'burnin',8);
opts=filldefault(opts,'minarea',2^20);
opts=filldefault(opts,'verbose',0);
opts=filldefault(opts,'stability_tol',10^-6);
opts=filldefault(opts,'metro',[]);
opts=filldefault(opts,strvcat('mACD','mBCD','mABCD','mABC','mAD'),0);
opts=filldefault(opts,'method',[]);
opts=filldefault(opts,'mtcopts',[]);
mtcopts=mtc_defopts(opts.mtcopts);
opts.mtcopts=mtcopts;
%
optsused=opts;
%
img=zeros([ard opts.nmaps]);
stats=[];
errs=[];
%
% set up calls to genmrfm or genmrfmg
%
og=[];
og.show=0;
og.statsub=0;
og.firstcol=[];
og.firstrow=[];
og.err_rept=opts.err_rept;
%
%determine number of gray levels
if isempty(opts.method)
    ng=2;
else
    p2x2=opts.method.p2x2;
    ng=size(p2x2,1);
end
%
if strcmp(opts.recurtype,'NoPickTT') %check stability
    if (ng==2)
        c=[];
        c.gamma=0;
        c.beta=[0 0 0 0];
        c.theta=[0 0 opts.mACD opts.mBCD];
        c.alpha=opts.mABCD;
        og.pblocks=getp2x2_corrs(c);
    else
        %here have to rotate p2x2 so that it is a "tu" configuration
        switch opts.method.variant_num
            case 1 %tu
                prot=[1 2 3 4];
            case 2 %tw
                prot=[3 1 4 2];
            case 3 %vw
                prot=[4 3 2 1];
            case 4 %uv
                prot=[2 4 1 3];
        end
        og.pblocks=permute(p2x2,prot);
    end
    optsused.pblock_2x2=og.pblocks;
    stabil_ok=1; %assume stable; if unstable will be set to 0, if unchecked, to -1
    p3x3_opts.ifshow=double((opts.verbose>1));
    if (ng<=mtcopts.maxng_3x3stabcheck)
        [p3x3,margs]=btc_3x3anal(og.pblocks,p3x3_opts); % calculate 3x3 block probabilities from 2x2 probs
        for k=1:size(margs.q2x2,5);
            if (max(max(max(max(abs(margs.q2x2(:,:,:,:,k)-og.pblocks)))))>opts.stability_tol)
                stabil_ok=0;
            end
        end
    else
        stabil_ok=-1;
        if (opts.verbose>0)
            disp(sprintf('btc_nopick: 3x3 stability check not done as ng=%2.0f, exceeding limit (%2.0f).',...
                ng,mtcopts.maxng_3x3stabcheck));
        end
    end
end
if strcmp(opts.recurtype,'NoPickBT')
    stabil_ok=1;
    if (ng==2)
        %use getp2x2_corrs to make T-region probabilities
        c=[];
        c.gamma=0;
        c.beta=[0 0 opts.mAD 0];
        c.theta=[opts.mAD.^2*opts.mABC 0 0 opts.mABC]; %change the opposite theta to induce stability
        c.alpha=0;
        og.pblocks=getp2x2_corrs(c);
    %
        [tc_new,tp_new]=btc_trule_iterate(btc_trule_getcorrs(og.pblocks));
        if (max(max(max(max(abs(tp_new-og.pblocks)))))>opts.stability_tol)
            stabil_ok=0;
        end
    else
        og.pblocks=opts.method.tee_probs_r;
        if (ng<=mtcopts.maxng_teestabcheck)
            [ptee_new,tee_margs]=mtc_trule_531(og.pblocks);
            for k=1:size(tee_margs,5);
                if (max(max(max(max(abs(tee_margs(:,:,:,:,k)-og.pblocks)))))>opts.stability_tol)
                    stabil_ok=0;
                end
            end
            if (max(abs(ptee_new(:)-og.pblocks(:)))>opts.stability_tol)
                stabil_ok=0;
            end
        else
            stabil_ok=-1;
            if (opts.verbose>0)
            disp(sprintf('btc_nopick: tee stability check not done as ng (%2.0f) exceeds limit (%2.0f).',...
                ng,mtcopts.maxng_teestabcheck));
            end
        end
    end
    optsused.pblock_tee=og.pblocks;
end
%
if (min(og.pblocks(:))<0)
    estring='requested correlations require negative probabilities';
    if (opts.err_rept==0)
        error(estring);
    else
        errs=strvcat(errs,estring);
        return
    end
end
%
%verify stability after iteration
%
if (stabil_ok==0)
    estring=sprintf('requested correlations result in unstable probabilities (%s) with tolerance %10.9f',opts.recurtype,opts.stability_tol);
    if (opts.verbose>0)
        disp(estring);
    end
    if (opts.err_rept==0)
        error(estring);
    else
        errs=strvcat(errs,estring);
        return
    end
else
    if (opts.verbose>0) & (stabil_ok==1)
        disp(sprintf('probabilities are stable (%s) with tolerance %10.9f',opts.recurtype,opts.stability_tol));
    end
end
% swap dimensions in case rotation is by 90 or 270 deg
if (mod(opts.rotvar,2)==0)
    ardg=[ard(2) ard(1)];
else
    ardg=ard;
end
[nr,nc]=nicesubp(opts.nmaps,ardg(2)/ardg(1));
layout=[nr*ardg(1) nc*ardg(2)];
ard_super=layout+opts.burnin;
if (opts.verbose)
    disp(sprintf('rotated map size  : %4.0f x %4.0f pixels',ardg));
    disp(sprintf('layout in supermap: %4.0f x %4.0f maps   (%5.0f requested)',nr,nc,opts.nmaps));
    disp(sprintf('supermap to make  : %4.0f x %4.0f pixels (including burnin of %4.0f pixels)',ard_super,opts.burnin));
end
% invoke the selected recursion rule
if strcmp(opts.recurtype,'NoPickTT')
    if (ng<=2)
        [imsuper,ss,ou,er1]=genmrfm(og,ard_super);
    else
        [imsuper,ss,ou]=genmrfmg(og,ard_super);
        er1=[];
    end
end
if strcmp(opts.recurtype,'NoPickBT')
     [imsuper,ss,ou]=genmrfmg_tee(og,ard_super);
     er1=[];
end
if ~isempty(er1)
    errs=strvcat(errs,cat(2,sprintf(' in map %5.0f',imap),': ',er1));
    if (opts.err_rept==0)
        return
    end
end
%
if opts.metro_show>0
    if opts.metro_show>1
        disp(sprintf('in btc_nopick ready for Metropolis algorithm, %s, with metro_opts:',opts.recurtype));
        disp(opts.metro_opts)
        disp(sprintf('map size: %5.0f x %5.0f',size(imsuper)));
    else
        disp(sprintf('in btc_nopick ready for Metropolis algorithm, %s.',opts.recurtype));
    end
end
donut=[];
donut.name='donut';
donut.matrix=[1 1 1;1 0 1; 1 1 1];
donut=glider_addcoords(donut);
if opts.onaxis==1
    disp('Metropolis algorithm not required, onaxis=1');
else
    if ~opts.metro_opts.numiters==0
        disp(sprintf('Calling Metropolis algorithm donut_metro with numiters=%7.0f',opts.metro_opts.numiters));
        [imsuper_new,metro_samp,metro_optsused]=donut_metro(ng,donut,imsuper,setfield(opts.metro_opts,'map_show_name',opts.recurtype));
        metro=[];
        metro.samp=metro_samp;
        metro.optsused=metro_optsused;
        metro.frac_changed=mean(abs(imsuper_new(:)-imsuper(:)));
        metro.numiters=opts.metro_opts.numiters;
        if opts.metro_show>0
            disp(sprintf('iterations: %5.0f, fraction of checks changed: %8.5f',...
                opts.metro_opts.numiters,metro.frac_changed));
        end
        imsuper=imsuper_new; %added 5 Oct 2015
    else
        disp('Metropolis algorithm not called, numiters=0');
    end
end
%cut up and rotate the results and put it into img
for imap=1:opts.nmaps
    ic=mod(imap-1,nc)+1;
    ir=floor((imap-1)/nc)+1;
    img(:,:,imap)=rot90(imsuper(opts.burnin+ardg(1)*(ir-1)+[1:ardg(1)],opts.burnin+ardg(2)*(ic-1)+[1:ardg(2)]),opts.rotvar-1);
end
%
%calculate stats
upx=ard;
if (opts.statsub==0) | (ng>2)
    stats=[];
else
    stats.bias_lum=0;
    stats.bias_eo=0;
    stats.bias_lum_sub=zeros(opts.statsub);
    stats.bias_eo_sub=zeros(opts.statsub);
    for imap=1:opts.nmaps
        [blum,beo]=getstats(img(:,:,imap));
        stats.bias_lum=stats.bias_lum+blum;
        stats.bias_eo=stats.bias_eo+beo;
        for isub=1:opts.statsub
            ilo=max(1,round(upx(1)*(isub-1)/opts.statsub));
            ihi=round(upx(1)*isub/opts.statsub);
            for jsub=1:opts.statsub
                jlo=max(1,round(upx(2)*(jsub-1)/opts.statsub));
                jhi=round(upx(2)*jsub/opts.statsub);
                [blum,beo]=getstats(img(ilo:ihi,jlo:jhi,imap));
                stats.bias_lum_sub(isub,jsub)=stats.bias_lum_sub(isub,jsub)+blum;
                stats.bias_eo_sub(isub,jsub)=stats.bias_eo_sub(isub,jsub)+beo;
            end
        end
    end %imap
    stats.bias_lum=stats.bias_lum/opts.nmaps;
    stats.bias_eo=stats.bias_eo/opts.nmaps;
    stats.bias_lum_sub=stats.bias_lum_sub/opts.nmaps;
    stats.bias_eo_sub=stats.bias_eo_sub/opts.nmaps;
end
%plot if requested
if (opts.show>0)
   figure;
   imagesc(img(:,:,1),[0 ng-1]);
   axis equal;axis off;colormap('gray');
end
return
%
function [bias_lum,bias_eo]=getstats(img)
bias_lum=2*mean(mean(img))-1;
img4=mod(img(2:end,2:end)+img(1:end-1,2:end)+img(2:end,1:end-1)+img(1:end-1,1:end-1),2);
bias_eo=1-2*mean(mean(img4));
return

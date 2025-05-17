function [thrs,specx,specs,opts_used]=btc_soidfg_predict(qform,vecs,dict,opts_fit)
% [thrs_pred,specx,specs,opts_used]=btc_soidfg_predict(qform,vecs,dict,opts_fit) calculates predicted thresholds from a quadratic 
% figure-ground model
% 
% This is similar in function to btc_soid_fit but arguments are quite different.
%
%input:
% qform: a quadratic form, size 2*length(dict.codel x 2*length(dict.codel),
%    i.e., 20 x 20; upper left block is (fig,fig), lower right is (gnd,gnd)
% vecs: one or more rows of length 2*dict.codel, NaN's where coordinates are to be filled in by btc_augcoords
% dict: output of btc_define
% opts_fit: output of btc_soidfg_define, can be omitted
%
%output:
% thrs: a multiplier for each element of vecs, so that vecs(k,:)'*qform*vecs=1, after augmentation of vecs
%    it opts_fit.ifaug=1, a nonlinear fit is used, via btc_soidfg_find.
%    if opts_fit.ifaug=0, threhold is calculated analytically
%  if qform is not positive definite, then thrs can be imaginary if if_aug=0, or NaN if if_aug=1
% specx: specx{ithr,ifg} is the structure of btc coords for figure (ifg=1) or ground (ifg=2) corresponding to threshold prediction
% specs: specs{ithr,ifg} is the structure of btc coords for figure (ifg=1) or ground (ifg=2) corresponding to original vecs(ithr,:)
% opts_used: options used
%
%   See also:  BTC_SOIDFG_FIT, BTC_SOIDFG_FIND, BTC_AUGCOORDS, BTC_AUGCOORDS_MULTS,
%   BTC_SOIDFG_DEFINE, BTC_SOIDFG_MODEL, BTC_DEFINE, BTC_SOID_PREDICT, BTC_VECNAN2LETCODE.
%
if (nargin<=3)
    opts_fit=[];
end
opts_fit=btc_soidfg_define(opts_fit);
opts_used=opts_fit;
btc_n=length(dict.codel);
thrs=NaN(size(vecs,1),1);
specs=cell(size(vecs,1),2);
for ithr=1:size(vecs,1)
    specs_fg=cell(1,2);
    for ifg=1:2
        stimvec=vecs(ithr,(ifg-1)*btc_n+[1:btc_n]);
        specs_fg{ifg}=btc_vecnan2letcode(stimvec,dict);
        specs{ithr,ifg}=specs_fg{ifg};
    end
    %determine a threshold without augmenting coordinates
    vec_unaug=vecs(ithr,:);
    vec_unaug(isnan(vec_unaug))=0;
    qval=vec_unaug*qform*vec_unaug';
    thr_unaug=1/sqrt(qval); %unaugmented threshold and initial guess for augmented threshold
    if (opts_fit.ifaug)
       if (qval>0) %if we get a reasonable threshold guess, use it -- otherwise force btc_soidfg_find to augment the coords
           thrs(ithr)=btc_soidfg_find(specs_fg,1,qform,dict,opts_fit,thr_unaug); 
       else
           thrs(ithr)=btc_soidfg_find(specs_fg,1,qform,dict,opts_fit,[]);
       end
    else
        thrs(ithr)=thr_unaug;
    end
    for ifg=1:2 %multiply each starting specification by the threshold
        specx{ithr,ifg}=btc_vecnan2letcode(thrs(ithr)*btc_letcode2vec(specs{ithr,ifg},dict),dict);
    end
end
return

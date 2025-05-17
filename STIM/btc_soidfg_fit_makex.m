function [x,regressors,augvecs_all,augvecs_sel]=btc_soidfg_fit_makex(model,thrs,coords,opts_fit,dict,isurr,aug_opts)
% [x,regressors,augvecs_all,augvecs_sel]=btc_soidfg_fit_makex(model,thrs,coords,opts_fit,dict,isurr,aug_opts)
% sets up variables for a regrssion for fitting figure-ground quadratic models
%
% This is modeled after btc_soid_fit_makex, captive in btc_soid_fit
%
%input:
% model: structure from btc_soidfg_model
% thrs: threshold data, one row for each threshold, [column of thresthresholds; [thr, eblo, ebhi]
%     if surrogates, then thr should be replaced by surrogate values
%     eblo and ebhi are optional and used to determine if error bars are present 
% coords: each row is [btc coords for figure btc coords for gnd] corresopnding to threshold of 1.  NaN's for values that need to be augmented.
% opts_fit: structure from btc_soidfg_define
% dict: structure from btc_define
% isurr: surrogate number, used for warning messages
% aug_opts: options for btc_augcoords_mults, can be omitted
%
%output:
% x: the design matrix
%    each row of x corresponds to a threshold measurement in augvecs_sel
%    each column of x corresponds to a degree of freedom in the quadratic form, i.e., a regressor of model
% regresssors: the names of the regressors, from model
% augvecs_all: coords, optionally augmented, for all thresholds that are non-NaN
% augvecs_sel: augvecs_all, selected for use in the model (will differ from augvecs_all if opts_fit.eb_need=1)
%
% See also:  BTC_SOIDFG_DEFINE, BTC_AUGCOORDS_MULTS, BTC_SOIDFG_MODEL, BTC_SOIDFG_FIT, BTC_VECNAN2LETCODE.
%
if (nargin<=6)
    aug_opts=[];
end
opts_fit=btc_soidfg_define(opts_fit);
btc_n=length(dict.codel);
regressors=model.param_name'; %make columm for consistency with btc_soid_fit_makex
aug_fails=zeros(1,2);
for ithr=1:size(thrs,1)
    mult=thrs(ithr,1);
    if ~isnan(mult)
        augvec=zeros(2,btc_n);
        for ifg=1:2
            stimvec=coords(ithr,btc_n*(ifg-1)+[1:btc_n]);
            if opts_fit.ifaug %augvec is augmented where stimulus coords are not specified
                spec=btc_vecnan2letcode(stimvec,dict);
                [augvec(ifg,:),augcoords_method,strats]=btc_augcoords_mults(spec,mult,dict,aug_opts);
                aug_fails(ifg)=aug_fails(ifg)+sum(double(strats<0));
            else %augvec has zeros where stimulus coords are not specified
                lets=find(~isnan(stimvec));
                augvec(ifg,lets)=mult*stimvec(lets);
            end
        end %ifg
    else
        augvec=NaN(2,btc_n);
    end
    augvecs_all(ithr,:)=[augvec(1,:) augvec(2,:)];
    %
end
if opts_fit.need_eb %threshold and error bars need to be non-NaN
    thr_sel=find(all(~isnan(thrs),2));
else %only threshold must be a non-NaN
    thr_sel=find(~isnan(thrs(:,1)));
end
augvecs_sel=augvecs_all(thr_sel,:);
if (opts_fit.verbose) | (isurr==0) | any(aug_fails>0)
    disp(sprintf(' surrogate %4.0f set up, %3.0f regressors, %3.0f thresholds total, %3.0f thresholds selected, augmentation failures (fig, gnd): %4.0f %4.0f',...
        isurr,size(regressors,1),size(augvecs_all,1),size(augvecs_sel,1),aug_fails));
end
%
%create the design matrix
%adapted from from btc_soid_fit_makex
x=zeros(size(augvecs_sel,1),size(regressors,1));
for nx=1:size(augvecs_sel,1)
    vec_outer=augvecs_sel(nx,:)'*augvecs_sel(nx,:); %outer product of augvec
    for ir=1:size(regressors,1)
        qones=model.qforms(:,:,ir);
        x(nx,ir)=sum(sum(vec_outer.*qones));
    end
end
%
return

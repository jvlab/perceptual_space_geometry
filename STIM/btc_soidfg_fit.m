function [results,opts_used]=btc_soidfg_fit(model,weib_a,coords,opts_fit,dict,aug_opts)
% [results,opts_used]=btc_soidfg_fit(model,weib_a,coords,opts_fit,dict,aug_opts)
% carries out the quadratic fit to figure-ground models
%
%input:
% model: model structure, from btv_soidfg_model
% weib_a: threshold data, one row for each threshold, [column of thresthresholds; [thr, eblo, ebhi]
% coords: each row is [btc coords for figure btc coords for gnd] corresopnding to threshold of 1.  NaN's for values that need to be augmented.
% opts_fit: structure from btc_soidfg_define
% dict: structure from btc_define
% aug_opts: options for btc_augcoords_mults, can be omitted
%
%output:
% results: structure containing regression matrices and results
%  and substructure surrogates, a stripped-down version for each surrogate
% opts_used: options used
%
% 15May20: add calculation of regression residuals, returned as rms
%
% See also:
% BTC_DEFINE, BTC_SOIDFG_DEFINE, BTC_SOIDFG_MODEL, BTC_AUGCOORDS_MULTS, BTC_SYMCLASSES,
% BTC_SOIDFG_FIT_MAKEX, BTC_SOIDFG_DEMO, BTC_SOID_FIT, REGRESS.
%
opts_fit=btc_soidfg_define(opts_fit);
opts_used=opts_fit;
if (nargin<=5)
    aug_opts=[];
end
btc_n=length(dict.codel);
%
eb_to_std_fixed=2*norminv(1-opts_fit.ebpval/2); %this CORRECTLY converts an error bar range (low to high) into a std dev;
%eb_to_std=2*norminv(1-opts_fit.ebpval); %INCORRECT from btc_soid_fit
%
results=struct();
surrogates=struct();
%calculate augmented coordinates and regressors
for isurr=0:opts_fit.nsurr
    if (isurr==0)
        rstring='';
        [x,regressors,augvecs_all,augvecs_sel]=btc_soidfg_fit_makex(model,weib_a,coords,opts_fit,dict,0,aug_opts); 
        results.x=x;
        results.regressors=regressors;
        results.augvecs_all=augvecs_all;
        results.augvecs_sel=augvecs_sel;
        %
        %condition numbers (from btc_soid_fit, only for real data)
        %
        condn_raw=cond(results.x);
        av=results.x;
        av_nzc=av(:,any(av~=0,1)); %remove any allzero-columns
        avc=av_nzc./repmat(sqrt(diag(av_nzc'*av_nzc))',size(av_nzc,1),1);
        condn_adjcol=cond(avc);
        av_nzr=av(any(av~=0,2),:); %remove any allzero-rows
        avr=av_nzr./repmat(sqrt(diag(av_nzr*av_nzr')),1,size(av_nzr,2));
        condn_adjrow=cond(avr);
        if (opts_fit.verbose>0)
            disp(sprintf('regression condition numbers: %8.3f (raw) %8.3f (column adjusted) %8.3f (row adjusted)',...
                condn_raw,condn_adjcol,condn_adjrow));
        end
        results.condn_raw=condn_raw;
        results.condn_adjcol=condn_adjcol;
        results.condn_adjrow=condn_adjrow;
        results.allzero_rows=size(av,1)-size(av_nzr,1);
        results.allzero_cols=size(av,2)-size(av_nzc,2);
    else %surrogates:
        if (isurr==1)
            disp(sprintf(' doing %4.0f surrogates',opts_fit.nsurr));
        end
        rstring=sprintf(', surrogate %4.0f',isurr);
        weib_a_surr=weib_a;
        ebars=weib_a_surr(:,3)-weib_a_surr(:,2);
        stdvs=abs(ebars)./eb_to_std_fixed;
        means=(weib_a_surr(:,2)+weib_a_surr(:,3))/2;
        weib_a_surr(:,1)=means+stdvs.*randn(size(stdvs,1),1);
        %sample from a Gaussian based on error bars
        [x,regressors,augvecs_all,augvecs_sel]=btc_soidfg_fit_makex(model,weib_a_surr,coords,opts_fit,dict,isurr,aug_opts); 
        surrogates.x(:,:,isurr)=x;
        % surrogates.regressors(:,:,isurr)=regressors; not needed since these are all the same
        surrogates.augvecs_all(:,:,isurr)=augvecs_all;
        surrogates.augvecs_sel(:,:,isurr)=augvecs_sel;
        %
    end
    %
    %check if any rows are all zeros (a threshold that can't be predicted)
    %and if any columns are all zeros (a parameter that doesn't matter)
    %
    if (opts_fit.verbose>0)
        disp(sprintf(' design matrix is %3.0f x %3.0f',size(av)));
    end
    count_zr=size(av,1)-size(av_nzr,1);
    count_zc=size(av,2)-size(av_nzc,2);
    if (opts_fit.verbose>0) | ((count_zr>0) | (count_zc>0))
        disp(sprintf(' design matrix has %3.0f all-zero rows (threshold can''t be reached) and %3.0f all-zero columns (parameter can''t be fit)',...
            count_zr,count_zc));
    end
    %
    % do the regression
    %
    if (isurr==1);hwait=waitbar(0,sprintf('Doing %2.0f surrogate regressions',opts_fit.nsurr));end
    if (isurr>0);waitbar(isurr/opts_fit.nsurr,hwait);end
    if (opts_fit.verbose>0) disp(sprintf('doing regression%s',rstring));end
    %
    [b,bint,resids]=regress(ones(size(x,1),1),x);
    rms_resids=sqrt(mean(resids.^2));
    %
    if (isurr==opts_fit.nsurr) & (isurr>0);close(hwait);end
    %
    % reconstitute the quadratic form from the regression matrix
    %
    qfit=zeros(btc_n*2,btc_n*2);
    for ic=1:length(results.regressors)
        qfit=qfit+b(ic)*model.qforms(:,:,ic);
    end
    if (isurr==0)
        results.b=b;
        results.qfit=qfit;
        results.rms_resids=rms_resids;
    else
        surrogates.b(:,isurr)=b;
        surrogates.qfit(:,:,isurr)=qfit;
        surrogates.rms_resids(:,isurr)=rms_resids;
    end
end %isurr
%
results.surrogates=surrogates;
return


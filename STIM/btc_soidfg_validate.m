function [ds_valid,counts,labels,eb_msg,opts_used]=btc_soidfg_validate(ds,opts_fit);
% [ds_valid,counts,labels,eb_msg,opts_used]=btc_soidfg_validate(ds,opts_fit) validates a plane of psychophysical data for
% figure-ground experiments, and fills in error bars when possible
%
% Error-bar filling uses same strategy as btc_soid_ebfix.
%
% ds:  a cell of d{:}, a database of figure-ground data, e.g., created by FIGGND_PROC
%      subj_id: 'MC'
%      files: {1×32 cell}
%      ndirs: 8
%       dirs: [1×1 struct]
%     weib_a: [8×3 double]
%     weib_b: [8×3 double]
% opts_fit: options, set up by btc_soidfg_define.  
%  Relevant fields:
%   th_min: minimum value for a non-exceptional threshold
%   th_max: maximum value for a non-exceptional threshold
%   eb_min: minimum value for a non-exceptional error bar
%   eb_max: maximum value for a non-exceptional error bar
%   eb_fill: how to fill in exceptional error bars, defaults to 1
%      0->do not fill
%      n->fill with a fraction of the threshold that is equal to the mean of
%      the n largest fractions (Inf means use all valid error bars in that plane)
%
% ds_valid: ds, with weib_a possibly adjusted: invalid thresholds or error
%   bars replaced by NaN's, and error bars filled in where possible
% counts: array in which each row corresponds to a threshold, indicating valid data, valid error bars, etc
% labels: cell array, describing columns of counts
% eb_msg: message from btc_soid_ebfix about how error bars were fixed
% opts_used: options used
%
%    See also:  BTC_SOIDFG_DEMO, BTC_SOID_EBFIX, BTC_SOIDFG_DEFINE.
%
if (nargin<=1)
    opts_fit=[];
end
opts_fit=btc_soidfg_define(opts_fit);
opts_used=opts_fit;
labels={'valid thresholds','out-of-range thresholds replaced by NaNs','valid error bars','filled-in error bars','error bars replaced by NaNs'};
counts=zeros(ds.ndirs,length(labels));
%assume all error bars and thresholds are valid
eb_msg='all valid';
ds_valid=ds;
th_valid=zeros(ds.ndirs,1);
%check that threhsolds are within range, and replace out-of-range thresholds by NaN's
th_valid=and(and(ds.weib_a(:,1)>=0,ds.weib_a(:,1)<=opts_fit.th_max),ds.weib_a(:,1)>=opts_fit.th_min);
counts(:,1)=th_valid;
counts(:,2)=1-th_valid;
ds.weib_a(th_valid==0,:)=NaN;
%find out-of-range error bars
eb_valid_lh=zeros(ds.ndirs,2);
for ilh=1:2 %1: low eb, 2: high eb
    eb_valid_lh(:,ilh)=and(and(ds.weib_a(:,1+ilh)>=0,ds.weib_a(:,1+ilh)<=opts_fit.eb_max),ds.weib_a(:,1+ilh)>=opts_fit.eb_min);
end
eb_valid=all(eb_valid_lh,2);
counts(:,3)=eb_valid;
%replace with Inf, to signal bt_soid_ebfix that error bars are invalid
ds.weib_a(eb_valid==0,[2:3])=Inf; %use repmat(ds.weib_a(eb_valid==0,1),1 2) to replace by 0
%convert weib_a to format needed by btc_soid_ebfix
edin=struct();
edin.plane.thresh_mags=ds.weib_a(:,1);
edin.plane.thresh_mags_eblo=ds.weib_a(:,2);
edin.plane.thresh_mags_ebhi=ds.weib_a(:,3);
%need to convert out-of-range error bars to zeros, so btc_soid_ebfix recognizes them as bad
[edfix,eb_avail,eb_msg]=btc_soid_ebfix(edin,opts_fit);
%convert back to weib_a
ds_valid.weib_a=[edfix.plane.thresh_mags,edfix.plane.thresh_mags_eblo,edfix.plane.thresh_mags_ebhi];
%identify filled-in error bars 
eb_ok=all(isfinite(ds_valid.weib_a(:,[2 3])),2);
counts(:,4)=eb_ok-eb_valid; %find the eb's that are present after btc_soid_ebfix but not before
counts(:,5)=1-eb_ok;
return

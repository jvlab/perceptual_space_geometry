function [qform_new,cust_factors]=btc_qform_customize(qform_orig,thr,dict)
% [qform_new,cust_factors]=btc_qform_customize(qform_orig,thr,dict) customizes a
% quadratic model tomatch a subject's on-axis thresholds
%
% qform_orig: a quadratic form threshold model
% thr: structure with fields g, b, d, t, a: subject's threhsolds on these axes
% dict: from btc_define.  Calculated if not provided.
%
%  See btc_qform_customize_notes.docx
%
%   See also: BTC_QFORM_CUSTOMIZE_TEST, BTC_QFORM_CUSTOMIZE_RUN,
%   BTC_DEFINE, BTC_AUGCOORDS, BTC_SOID_FIND.
%
if (nargin<3)
    dict=btc_define;
end
nbtc=length(dict.codel);
%
adj_order='atdbg'; %order to adjust gains
c0=0.5; %contrast for probing augmented coordinates
%
axes_to_fix=fieldnames(thr);
nfix=length(axes_to_fix);
%
aug_pwr=struct; %how each statistic induces correlations on other axes
sym_vec=struct; %how each statistic maps to a symmetry class, via dict.name_order_aug
for k=1:nfix
    ax=axes_to_fix{k};
    augcoords=btc_augcoords(setfield(struct(),ax,c0),dict);
    vec=augcoords.method{1}.vec;
    aug_pwr.(ax)=zeros(1,nbtc);
    aug_pwr.(ax)(vec~=0)=round(log(vec(vec~=0))./log(c0));
    name_order_match=dict.name_order_aug{find(dict.codel==ax)};
    sym_vec.(ax)=zeros(1,nbtc);
    sym_vec.(ax)(strmatch(name_order_match,dict.name_order_aug,'exact'))=1;
end
% disp('aug_pwr');
% disp(aug_pwr);
% disp('sym_vec');
% disp(sym_vec);
%adjust qform, and will need to do this for the surrogates
%work in direction of adj_order, successively modifying qform
qform_new=qform_orig;
cust_factors=struct;
for k=1:length(adj_order)
    ax=adj_order(k);
    pwr_ax=double(dict.codel==ax);
    pwr_aug=aug_pwr.(ax)-pwr_ax;
    vec_ax=thr.(ax)*pwr_ax';
    vec_aug=zeros(1,nbtc)';
    vec_aug(pwr_aug~=0)=thr.(ax).^(pwr_aug(pwr_aug~=0))';
    qa=vec_aug'*qform_new*vec_aug; %interaction of induced stats with themselves
    qb=vec_aug'*qform_new*vec_ax; %half of the interaction of ax with induced stats
    qc=vec_ax'*qform_new*vec_ax; %interaction ax with itself
%    disp(sprintf(' ax qa qb qc : %s, %12.7f, %12.7f, %12.7f',ax,qa,qb,qc));
    if (qb==0) & (qa==0)
        cust_factors.(ax)=1/sqrt(qc);
    else %solve qc*x^2+2(qb)x+qa=1
        cust_factors.(ax)=(-qb+sqrt((qb)^2-qc*(qa-1)))/qc;
    end
    %create a diagonal matrix with all of the stats of the same
    %symmetry class as ax using dict.name_order_aug
    %and then pre- and post-multiply qform_new by this
    d=ones(1,nbtc);
    d(sym_vec.(ax)==1)=cust_factors.(ax);
    d=diag(d);
    qform_new=d*qform_new*d;
end
return

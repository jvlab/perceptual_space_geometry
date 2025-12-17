function [amin,amax,p2x2_extremes]=btc_alpharange(p2x2)
% [amin,amax,p2x2_extremes]=btc_alpharange(p2x2) determines the minimum and
% maximum values of the 4-point statistic alpha that is consistent with
% the lower-order statistics implied by a 2x2 recursion matrix
%
% Compatibility of the recursion matrix with Pickard conditions is assumed,
% as is sum(p2x2(:))=1, and p2x2(:)>=0.
% 
% p2x2: Recursion matrix, size [2 2 2 2]
%
% amin, amax: min and max compatible alpha values
% p2x2_extremes: p2x2_extremes(:,:,:,:,im) is the 2x2 recursion for the minimum alpha (im=1) or maximum alpha (im=2)
%
%   See also:  BTC_DEFINE, BTC_ALPHARANGE_BG_DEMO.
% 
q2=[1 -1;-1 1];
q3=cat(3,q2,-q2);
q4=cat(4,q3,-q3);
qpos=find(q4>0);
qneg=find(q4<0);
%
m_max=min(min(1-p2x2(qpos)),min(p2x2(qneg)));
m_min=max(max(p2x2(qneg)-1),max(-p2x2(qpos)));
p2x2_extremes(:,:,:,:,1)=m_min*q4+p2x2;
p2x2_extremes(:,:,:,:,2)=m_max*q4+p2x2;
amin=sum(reshape(p2x2_extremes(:,:,:,:,1).*q4,16,1));
amax=sum(reshape(p2x2_extremes(:,:,:,:,2).*q4,16,1));
return
end

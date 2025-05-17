function [dsq,qh]=btc_soidfg_anal_of(h,q)
% [dsq,qh]=btc_soidfg_anal_of(h,q) is the objective function for btc_soidfg_anal
%   h: scalar
%   q: a model quadratic form, typically symmetric, size(q) is even
%
%   dsq: squared deviation of q-qh (see below)
%   qh: projection of q onto [h -1/2;-1/2 h] *ones(size(q)/2)
%
%   See also:  BTC_SOIDF_ANAL.
%
n2=length(q);
n=n2/2;
nf=[1:n];
ng=[n+1:n2];
%reorganize q into quadruples for each matrix position
qstack=cat(3,q(nf,nf),q(ng,nf),q(nf,ng),q(ng,ng));
qreorg=reshape(permute(qstack,[3 1 2]),[4,n^2]);
hvec=[h -1/2 -1/2 1-h];
hnorm_sq=sum(hvec.^2);
qh_reorg=hvec*qreorg/hnorm_sq;
qh_block=reshape(qh_reorg,[n n]);
qh=[qh_block*h -qh_block/2;-qh_block/2 qh_block*(1-h)];
dsq=sum((q(:)-qh(:)).^2);
return
 
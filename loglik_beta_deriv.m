function dloglik=loglik_beta_deriv(ab,obs)
% dloglik=loglik_beta_deriv(ab,obs) computes the derivatives of the log likelihood 
% (natural log) of a pair of parameters [a,b] of a beta distribution that
% would yield a set of observations
%
% also see .../jv/ey07977/psg_umi_notes.doc.
%
% a and b must be >0.
%
% ab: [a b], or, if scalar, the common value of a and b.  Note a=b=1 for flat prior.
% obs: a vector of observations where obs(:,1) is successes and obs(:,2) is total tries
%
% dloglik: derivs of the log likelihood w.r.t. a and b
%
% Note that if ab is a scalar and the derivative of the log likelhiood
% function with respect to the joint exponent (p^(a-1)*(1-p)^(a-1)) is
% desired, then this is d/da + d/db of the deriv of (p^(a-1)*(1-p)^(b-1)),
% i.e., the sum of dloglik
%
%   See also:  LOGLIK_BETA_DEMO, BETALN, FILLDEFAULT.
%
if length(ab)==1
    ab=[ab ab];
end
dloglik=zeros(1,2);
n=[obs(:,1),obs(:,2)-obs(:,1)];
tmax=max(obs(:,2));
if (tmax<1)
    return
end
recipsum_combined=[0 cumsum(1./(sum(ab)+[0:tmax-1]))]; %0-offset sum of recips
recipsum_each=cell(1,2);
for iv=1:2
    recipsum_each{iv}=[0 cumsum(1./(ab(iv)+[0:(max(n(:,iv))-1)]))];
end
for iv=1:2
    for k=1:size(n,1)
        dloglik(iv)=dloglik(iv)+recipsum_each{iv}(1+n(k,iv))-recipsum_combined(1+obs(k,2));
    end
end
return

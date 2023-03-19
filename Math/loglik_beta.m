function [loglik,opts_used]=loglik_beta(ab,obs,opts)
% [loglik,opts_used]=loglik_beta(ab,obs,opts) determines the log likelihood
% (natural log) of a pair of parameters [a,b] of a beta distribution (Dirichlet prior)
% that would yield a set of observations.
%
% Also has options for adding a set of point masses if obs is [successes trials]
%
% also see .../jv/ey07977/psg_umi_notes.doc.
%
% a and b must be >0.
% if any of obs are 0 or 1, results may be Inf, -Inf, or NaN 
%
% ab: [a b], or, if scalar, the common value of a and b.  Note a=b=1 for flat prior.
% obs: a column vector of probabilities
%    or, a vector of observations where obs(:,1) is successes and obs(:,2) is total tries
%    There are no exceptions in this case,  because of regularization due
%    to the beta-integration
% opts: options
%    opts.if_norm:
%        1: normalize by -log(beta(a,b))*length(non-exceptional obs)
%       -1: compute normalization and return in opts_used but do not include in loglik
%        0: do not compute normalization
%        Note that opts.if_norm should always be set to 1 if opts.hvec is nonempty or sums to >0
%    opts.hvec: vector of weights for point-mass probabilities
%       must be all >=0 and sum to <=1 (not checked for), defaults to [];
%    opts.qvec: vector of probabilities corresponding to the elements of hvec,
%       must all be strictly >0 and < 1 (not checked for)
%    Notes re hvec,qvec: 
%       * opts.qvec and opts.hvec ignored if obs is probabilities
%       * if hvec, qvec present, opts.if_norm ignored; normalization always computed and incorporated into loglik
% Re normalization (for size(obs,2)=2, obs is counts):  
%   the likelihood does not take into account other sequences of trials
%   that have the same counts.  So for a point-mass at 0.5 (qvec=.5,hvec=1),
%   log likelihood depends for obs=[k,n] depends only on n.
%   Same is true for continuous part only (hvec=0), with beta-parameters equal and large, which approximates
%   a point mass at 0.5.  See loglik_beta_test.
%
% loglik: the log likelihood
% opts_used: options used
%
% 09Jan23: add opts.hvec,opts.qvec
% 10Jan23: add checking to see if betaln arguments are <0 for finite-obs case
% 17Mar23: added documentation
%
%   See also:  LOGLIK_BETA_DEMO, LOGLIK_BETA_DERIV, LOGLIK_BETA_EXCEPT, BETALN, FILLDEFAULT, PBETABAYES_COMPARE,
%    LOGLIK_BETA_DEMO2, LOGLIK_BETA_TEST.
%
exception_vals.big=Inf;
exception_vals.small=-Inf;
if (nargin<3)
    opts=struct();
end
opts=filldefault(opts,'if_norm',1);
opts=filldefault(opts,'hvec',[]);
opts=filldefault(opts,'qvec',[]);
if size(obs,2)==1 & (~isempty(opts.hvec) | ~isempty(opts.qvec))
    warning('point-mass parameters (hvec and qvec) ignored for calculation based on probabilities');
    opts.hvec=[];
    opts.qvec=[];
end
if sum(opts.hvec)>0 & opts.if_norm~=1
    warning('un-normalized calculation not allowed with combined Dirichlet and point-mass prior, normalizatio turned on.');
    opts.if_norm=1;
end
%
if length(ab)==1
    ab=[ab ab];
end
a=ab(1);
b=ab(2);
nobs=size(obs,1);
if size(obs,2)==1 %obs are probabilities, hvec and qvec ignored
    p=obs;
    inside_ptr=find(p>0 & p<1);
    p_inside=p(inside_ptr); %can be non-exceptional if p=0 or 1 and a=b=1.
    nobs=length(p_inside);
    %
    exceptions=loglik_beta_except(a,b,p);
    if isempty(exceptions)
        if isempty(p_inside)
            loglik_beta=0;
        else
            loglik_beta=(a-1)*sum(log(p_inside))+(b-1)*sum(log(1-p_inside));
        end
    elseif length(exceptions)==2
        loglik_beta=NaN;
    else
        loglik_beta=exception_vals.(exceptions{1});
    end
else %obs is [successes tries]
    loglik_beta_each=zeros(nobs,1);
    exceptions=[];
    if nobs>0
        a_args=a+obs(:,1);
        b_args=b+obs(:,2)-obs(:,1);
        if any(a_args<=0) | any(b_args<=0)
            loglik_beta_each=Inf;
        else
            loglik_beta_each=betaln(a_args,b_args); %contribution of each non-exceptional beta term
        end
    end
    loglik_beta=sum(loglik_beta_each);
end
if opts.if_norm~=0
    if (a<=0) | (b<=0)
        lognorm_beta_each=-Inf;
    else
        lognorm_beta_each=-betaln(a,b);
    end
    lognorm_beta=lognorm_beta_each*nobs;
else
    lognorm_beta=[];
end
if sum(opts.hvec)>0
    lik_beta_each=exp(loglik_beta_each+lognorm_beta_each);% beta-component, with normalization
    lik_pointmass_each=ones(nobs,length(opts.hvec));
    for hptr=1:length(opts.hvec)
        q=opts.qvec(hptr);
        lik_pointmass_each(:,hptr)=q.^(obs(:,1)).*(1-q).^(obs(:,2)-obs(:,1)); %contribution of point-mass
    end
    lik_total_each=lik_beta_each*(1-sum(opts.hvec))+lik_pointmass_each*opts.hvec(:); %weighted sum
    loglik=sum(log(lik_total_each));
    lognorm=0;
else
    lognorm=lognorm_beta;
    loglik=loglik_beta;
end
if opts.if_norm>0
    loglik=loglik+lognorm;
end
opts.lognorm=lognorm;
opts.exceptions=exceptions;
opts_used=opts;
return

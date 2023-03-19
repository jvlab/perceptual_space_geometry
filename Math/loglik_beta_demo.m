% loglik_beta_demo demonstrates random draws from a symmetric beta distribution and
% fits to the beta-function parameters via maximum-likelihood
% by direct optimization, and by setting derivative to zero
%
% also see .../jv/ey07977/psg_umi_notes.doc.
%
%   See also: LOGLIK_BETA, LOGLIK_BETA_DERIV, LOGLIK_BETA_DEMO2.
%
rng('default');
if ~exist('limits') limits=[2^-10 2^3]; end
if ~exist('a_try') a_try=2.^[-10:.0625:10]; end
a_true=getinp('true value of beta parameter','f',[.1 10],.3);
nsamps=getinp('number of samples','d',[10 10^6],10^3);
if ~exist('ntries_list') ntries_list=[5 10]; end
ntries_list=getinp('ntries_list for finite samples','d',[1 1000],ntries_list);
%
pvals=betaincinv(rand(nsamps,1),a_true,a_true);
ll_true=loglik_beta(a_true,pvals(:));
%
figure;
set(gcf,'Position',[50 100 1200 700]);
subplot(3,1,1);
hist(pvals,100);
set(gca,'XLim',[0 1]);
xlabel('p');
title(sprintf('a_t_r_u_e=%6.3f',a_true));
%
ll_try=zeros(1,length(a_try));
for itry=1:length(a_try)
    ll_try(itry)=loglik_beta(a_try(itry),pvals(:));
end
ind_max=min(find(ll_try==max(ll_try)));
subplot(3,1,2);
semilogx(a_try,ll_try,'k');
hold on;
semilogx(a_true,ll_true,'ko');
semilogx(a_try(ind_max),ll_try(ind_max),'kx');
set(gca,'XLim',[min(a_try) max(a_try)]);
set(gca,'YLim',sort(ll_try(ind_max)*[.1 1.1]));
title(sprintf('a_m_a_x=%6.3f',a_try(ind_max)));
xlabel('a');
ylabel('ll');
%
[a_best,nll_best,exitflag]=fminbnd(@(x) -loglik_beta(x,pvals(:)),limits(1),limits(2));
semilogx(a_best,-nll_best,'ks');
legend({'all','true','max try','opt'});
%
%now optimize based on finite samples
%
tries=ntries_list(ceil(length(ntries_list)*rand(length(pvals),1)));
tries=tries(:);
successes=binornd(tries,pvals(:));
obs=[successes tries];
ll_true_samp=loglik_beta(a_true,obs);
ll_try_samp=zeros(1,length(a_try));
for itry=1:length(a_try)
    ll_try_samp(itry)=loglik_beta(a_try(itry),obs);
end
ind_max_samp=min(find(ll_try_samp==max(ll_try_samp)));
subplot(3,1,3);
semilogx(a_try,ll_try_samp,'k');
hold on;
semilogx(a_true,ll_true_samp,'ko');
semilogx(a_try(ind_max_samp),ll_try_samp(ind_max_samp),'kx');
set(gca,'XLim',[min(a_try) max(a_try)]);
set(gca,'YLim',sort(ll_try_samp(ind_max_samp)*[.1 1.1]));
title(sprintf('finite samples: a_m_a_x=%6.3f',a_try(ind_max_samp)));
xlabel('a');
ylabel('ll');
%
[a_best_samp,nll_best_samp,exitflag]=fminbnd(@(x) -loglik_beta(x,[successes tries]),limits(1),limits(2));
semilogx(a_best_samp,-nll_best_samp,'ks');
legend({'all','true','max try','opt'});
%
disp(sprintf('     true a                                   is %7.3f, with log likelihood %9.3f',a_true,ll_true));
disp(sprintf(' best-fit a via optimization                  is %7.3f, with log likelihood %9.3f',a_best,-nll_best));
disp(sprintf(' best-fit a via optimization (finite samples) is %7.3f, with log likelihood %9.3f',a_best_samp,-nll_best_samp));
%
%verify derivatives
%
a_best_tries=[0.9*a_best_samp,a_best_samp,0.1+0.9*a_best_samp];
for ia=1:3
    disp(sprintf(' for a=%12.5f, deriv of log likelihood is %12.5f',a_best_tries(ia),sum(loglik_beta_deriv(a_best_tries(ia),obs))));
end
%
%find minimum by derivative
%
disp('derivative method, finite samples');
a_deriv_samp=fzero(@(x) sum(loglik_beta_deriv(x,obs)),1);
ll_a_deriv_samp=loglik_beta(a_deriv_samp,obs);
disp(sprintf(' best-fit a via optimization (finite samples) is %7.3f, with log likelihood %9.3f',a_deriv_samp,ll_a_deriv_samp));
a_deriv_tries=[0.9*a_deriv_samp,a_deriv_samp,0.1+0.9*a_deriv_samp];
for ia=1:3
    disp(sprintf(' for a=%12.5f, deriv of log likelihood is %12.5f',a_deriv_tries(ia),sum(loglik_beta_deriv(a_deriv_tries(ia),obs))));
end





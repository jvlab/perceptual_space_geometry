% loglik_beta_demo2 demonstrates random draws from a symmetric beta distribution 
% with optional point mass, known to be at 0.5, and fits to the beta-function 
% and point-mass parameters via maximum-likelihood.
%
%  Compared with loglik_beta_demo:
%     * Only does finite-sample scenario.
%     * Uses fminbnd to find the beta parameter, then
%       fminsearch, starting with that value of the beta parameter, to find
%       best beta and mass
%
% also see .../jv/ey07977/psg_umi_notes.doc.
%
%   See also: LOGLIK_BETA, LOGLIK_BETA_DERIV, LOGLIK_BETA_DEMO.
%
rng('default');
if ~exist('a_limits') a_limits=[2^-10 2^3]; end
if ~exist('h_limits') h_limits=[0 1]; end
if ~exist('a_try') a_try=2.^[-10:.0625:10]; end
if ~exist('q_limits') q_limits=[0.001 0.999]; end
if ~exist('h_init') h_init=0; end %initial value for h
a_true=getinp('true value of beta parameter','f',[.1 10],0.3);
h_true=getinp('true value of weight parameter','f',[0 1],0.2);
q_true=getinp('true value of q','f',[0 1],0.5);
%
if ~exist('opts_disc')
    opts_disc=struct;
end
opts_disc.hvec=h_true;
opts_disc.qvec=q_true;
%
nsamps=getinp('number of samples','d',[10 10^6],10^3);
if ~exist('ntries_list') ntries_list=[5 10]; end
ntries_list=getinp('ntries_list for finite samples','d',[1 1000],ntries_list);
%
pvals_beta=betaincinv(rand(nsamps,1),a_true,a_true);
beta_vs_disc=find(double(rand(nsamps,1)>=h_true));
pvals=repmat(q_true,nsamps,1); %assume all discrete
pvals(beta_vs_disc)=pvals_beta(beta_vs_disc);
%
figure;
set(gcf,'Position',[50 100 1200 700]);
subplot(2,1,1);
hist(pvals,100);
set(gca,'XLim',[0 1]);
xlabel('p');
title(sprintf('a_t_r_u_e=%6.3f, h_t_r_u_e=%6.3f, q_t_r_u_e=%6.3f',a_true,h_true,q_true));
%
%create finite-sample dataset
%
tries=ntries_list(ceil(length(ntries_list)*rand(length(pvals),1)));
tries=tries(:);
successes=binornd(tries,pvals(:));
obs=[successes tries];
%
ll_true_samp_beta=loglik_beta(a_true,obs);
ll_true_samp_disc=loglik_beta(a_true,obs,opts_disc); %log likelihood calc with discrete part, with h and q
ll_try_samp_beta=zeros(1,length(a_try));
ll_try_samp_disc=zeros(1,length(a_try));
for itry=1:length(a_try)
    ll_try_samp_beta(itry)=loglik_beta(a_try(itry),obs);
    ll_try_samp_disc(itry)=loglik_beta(a_try(itry),obs,opts_disc);  %log likelihood calc with discrete part, with h and q
end
ind_max_samp_beta=min(find(ll_try_samp_beta==max(ll_try_samp_beta)));
ind_max_samp_disc=min(find(ll_try_samp_disc==max(ll_try_samp_disc)));
%
%optimize with to the samples, leaving out the discrete part
%
[a_best_samp_beta,nll_best_samp_beta,exitflag_beta]=fminbnd(@(x) -loglik_beta(x,[successes tries]),a_limits(1),a_limits(2));
%
%optimize to the finite samples, assuming discrete part known
%
[a_best_samp_disc,nll_best_samp_disc,exitflag_disc]=fminbnd(@(x) -loglik_beta(x,[successes tries],opts_disc),a_limits(1),a_limits(2));
%
%now fit both, using fitted a as starting point
%
ah_init=[a_best_samp_beta;h_init];
[ah_best_samp,nll_best_ah,exitflag_ah,output_ah]=fminsearch(@(x) -loglik_beta(x(1),[successes tries],setfield(opts_disc,'hvec',x(2))),ah_init);
%
%plot
%
subplot(2,1,2);
semilogx(a_try,ll_try_samp_beta,'k');
hold on;
semilogx(a_true,ll_true_samp_beta,'ko');
semilogx(a_try(ind_max_samp_beta),ll_try_samp_beta(ind_max_samp_beta),'kx');
semilogx(a_best_samp_beta,-nll_best_samp_beta,'ks');
%
semilogx(a_try,ll_try_samp_disc,'r');
semilogx(a_true,ll_true_samp_disc,'ro');
semilogx(a_try(ind_max_samp_disc),ll_try_samp_disc(ind_max_samp_disc),'rx');
semilogx(a_best_samp_disc,-nll_best_samp_disc,'rs');
%
semilogx(ah_best_samp(1),-nll_best_ah,'b*');
%
set(gca,'XLim',[min(a_try) max(a_try)]);
yvals=[ll_try_samp_disc,ll_try_samp_beta,ll_true_samp_disc,ll_true_samp_beta];
%set(gca,'YLim',sort(ll_try_samp_disc(ind_max_samp_disc)*[.1 1.1]));
set(gca,'YLim',[min(yvals) max(yvals)*0.9]);
title('finite samples, fit a with and without discrete part');
xlabel('a');
ylabel('ll');
%
legend({'beta only','true a','best a','fit a','beta + disc','true a','best a','fit a',sprintf('fit both, h=%4.2f',ah_best_samp(2))},'FontSize',7);
%
disp(sprintf(' beta only:  true a                                       is %7.3f, with log likelihood %9.3f',a_true,ll_true_samp_beta));
disp(sprintf(' beta only:  best-fit a via optimization (finite samples) is %7.3f, with log likelihood %9.3f',a_best_samp_beta,-nll_best_samp_beta));
disp(sprintf(' assume h,q: true a                                       is %7.3f, with log likelihood %9.3f',a_true,ll_true_samp_disc));
disp(sprintf(' assume h,q: best-fit a via optimization (finite samples) is %7.3f, with log likelihood %9.3f',a_best_samp_disc,-nll_best_samp_disc));
disp(sprintf(' assume   q: best-fit a via optimization (finite samples) is %7.3f, with log likelihood %9.3f, h=%5.3f',ah_best_samp(1),-nll_best_ah,ah_best_samp(2)));

function [ncloser_resample,choice_probs,opts_used]=psg_resample_conform(params,ntrials,partitions,nprobs,ndraws,opts)
%[ncloser_resample,choice_probs,opts_used]=psg_resample_conform(params,ntrials,partitions,nprobs,ndraws,opts)
% simulates a set of similarity judgments taken from underlying probabilities that conform to a specific set of
% rank-choice inequalities, with the probabilities drawn from a Dirichlet distribution with a discrete component
%
% ntrials: size [nt,ncomps], nt=number of triads or tents, ncomps=3 for triads, 6 for tents
% params: parameters of the distribution of probabilities
%   params.a is the Dirichlet exponent
%   params.h is the mass of the discrete component, defaults to 0
% partitions: array of ncomps dimensions, size 3 on each, indicating the
%    excluded regions, typically computed by psg_ineq_logic with arguments 
%    for triads: 'exclude_sym' or 'exclude_umi_trans'
%    for tents: 'exclude_addtree_trans'
% nprobs: number of times to resample the probabilities
% ndraws: number of times to simulate choices for each draw of probabilities
% opts:
%    opts.if_log: 1 to log
%
% ncloser_resample:  simulated judgments, size=[nt ncomps nprobs ndraws]
%      0<=ncloser_resample(it,icomp,iprob,idraw)<=ntrials(it,icomp)
% choice_probs: choice probabilities for each draw, size=[nt ncomps nprobs]
% opts_used: options used
%    nondiscrete_allowed: 1 if an all-continuous distribution can be selected
%    nattempts: number of attempts for each draw, size=[nt nprobs]
%
%  See also:  PSG_UMI_TRIPLIKE_DEMO, PSG_TENTLIKE_DEMO, PSG_INEQ_LOGIC, PSG_CONFORM, PSG_RESAMPLE_CONFORM_TEST.
%
if (nargin<=5)
    opts=struct;
end
opts=filldefault(opts,'if_log',1);
params=filldefault(params,'h',0);
%
nt=size(ntrials,1);
ncomps=size(ntrials,2);
discrete_val=0.5; %discrete part of the distribution
%
ncloser_resample=zeros(nt,ncomps,nprobs,ndraws);
choice_probs=zeros(nt,ncomps,nprobs);
%
if_warn=0;
if ndims(partitions)~=ncomps
    warning('partitions has wrong number of dimensions');
    if_warn=1;
end
if length(partitions(:))~=3^ncomps
    warning('partitions is wrong size');
    if_warn=1;
end
%
%create a template of the number of edges in each location
%
z=[0 1 0];
for k=2:ncomps
    z=[z z+1 z];
end
zbox=reshape(z,[repmat(3,1,ncomps) 1]);
%
%check whether partitions requires values in the discrete part and h=0
%
nondiscrete_allowed=length(intersect(find(partitions(:)==0),find(z(:)==0)));
opts.nondiscrete_allowed=nondiscrete_allowed;
if (params.h==0 & nondiscrete_allowed==0)
    warning('distribution has no discrete part but this is required by partition logic');
    opts_used=opts;
    if_warn=1;
end
if (if_warn>0)
    return
end
%
mults=3.^[0:ncomps-1];
nattempts=zeros(nt,nprobs);
for it=1:nt
    for iprob=1:nprobs
        ifok=0;
        while (ifok==0)
            nattempts(it,iprob)=nattempts(it,iprob)+1;
            %generate a set of probabilities
            probs=zeros(1,ncomps);
            probs=betainv(rand(1,ncomps),params.a,params.a); %assume non-discrete
            discrete=find(rand(1,ncomps)<params.h);
            probs(discrete)=discrete_val;
            %test for conformance to logic
            pd=1+sign(probs-discrete_val); %0,1, or 2 depending on whether probs is <1/2, =1/2, or > 1/2
            if partitions(1+sum(mults.*pd))==0
                ifok=1;
            end
        end
        choice_probs(it,:,iprob)=probs;
        %use probs to simulate ndraws trials
        for icomp=1:ncomps
            ncloser_resample(it,icomp,iprob,:)=binornd(ntrials(it,icomp),probs(icomp),[1 1 1 ndraws]);
        end
    end %iprob
end %it
opts.nattempts=nattempts;
%
opts_used=opts;
return

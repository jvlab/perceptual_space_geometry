function [ncloser_conform,if_flip,opts_used]=psg_conform(ncloser,ntrials,partitions,opts)
% [ncloser_conform,if_flip,opts_used]=psg_conform(ncloser,ntrials,partitions,opts)
% modifies a psg dataset of choices to conform to a specific set of
% rank-choice inequalities
%
% ncloser: size [nt,ncomps] ncomps=3 for triads, 6 for tents
%    each of the nt rows corresponds to a separate triad or tent
%    entries list number of times that the first alternative was chosen as closer
% ntrials: size [nt,ncomps], number of trials
% partitions: array of ncomps dimensions, size 3 on each, indicating the
%    excluded regions, typically computed by psg_ineq_logic with arguments 
%    for triads: 'exclude_sym' or 'exclude_umi_trans'
%    for tents: 'exclude_addtree_trans'
% opts:
%    opts.if_log: 1 to log
%    opts.method: 
%          'flip_one': flip one response if that will move the triad or tent into the allowed region
%          'none': do nothing
%    opts.penalty: 
%          'chi2': penalty for flipping is chi-squared table (default)
%          'nobs': penatly for flipping is number of observations flipped
%
% ncloser_conform:  modified ncloser, as per above
% if_flip: size [nt,ncomps]: 1 if a response is to be flipped, 0 otherwise
% opts_used: options used
%    opts_used.results: some statistics
%
%  See also:  PSG_UMI_TRIPLIKE_DEMO, PSG_TENTLIKE_DEMO, PSG_INEQ_LOGIC.
%
if (nargin<=3)
    opts=struct;
end
opts=filldefault(opts,'method','flip_one');
opts=filldefault(opts,'if_log',1);
opts=filldefault(opts,'penalty','chi2');
%
results=[];
%
nt=size(ncloser,1);
ncomps=size(ncloser,2);
ncloser_conform=ncloser;
if_warn=0;
if ndims(partitions)~=ncomps
    warning('partitions has wrong number of dimensions');
    if_warn=1;
end
if length(partitions(:))~=3^ncomps
    warning('partitions is wrong size');
    if_warn=1;
end
if size(ntrials,1)~=nt | size(ntrials,2)~=ncomps
    warning('ncloser and ntrials incompatible');
    if_warn=1;
end
%
if (if_warn)
    if_flip=[];
    ncloser_conform=[];
else
    ncloser_conform=ncloser;
    if_flip=zeros(nt,ncomps);
    %determine which orthant each set of observations is in
    signs=1+sign(ncloser-ntrials/2); %-1 if < 1/2, 0 if 1/2, 1 if > 1/2
    powers3=3.^[0:ncomps-1];
    multi_index=1+signs*powers3';
    tokens_orig=partitions(multi_index);
    if (opts.if_log)
        disp(sprintf('analyzing %6.0f sets of %2.0f judgments',nt,ncomps));
        disp(sprintf('     with no flips: %6.0f sets in excluded regions',sum(tokens_orig==1)));
    end
    %consider all possible one-judgment flips
    tokens_oneflip=zeros(nt,ncomps);
    results.n_exclude_eachflip=zeros(1,ncomps);
    for icomp=1:ncomps
        signs_flip=signs;
        signs_flip(:,icomp)=2-signs(:,icomp);
        multi_index_flip=1+signs_flip*powers3';
        tokens_oneflip(:,icomp)=partitions(multi_index_flip);
        results.n_exclude_eachflip(1,icomp)=sum(tokens_oneflip(:,icomp));
    end
    results.n_exclude_anyflip=sum(min(tokens_oneflip,[],2));
    canflip=intersect(find(any(tokens_oneflip==0,2)),find(tokens_orig==1));
    nflips=sum(tokens_oneflip==0,2);
    results.n_oneflip=length(canflip);
    canflip_oneway=intersect(canflip,find(nflips==1));
    results.n_oneflip_oneway=length(canflip_oneway);
    if (opts.if_log)
        for icomp=1:ncomps
            disp(sprintf('   flipping obs %2.0f: %6.0f sets in excluded regions',icomp,results.n_exclude_eachflip(1,icomp)));
        end
        disp(sprintf(' best one-obs flip: %6.0f sets in excluded regions',results.n_exclude_anyflip));
        disp(sprintf('                for %6.0f sets, one flip can move out of excluded regions',results.n_oneflip));
        disp(sprintf('                for %6.0f of these sets, there is just one way',results.n_oneflip_oneway));
    end    
    switch opts.penalty
        case 'chi2'
            penalties=2*(2*ncloser-ntrials).^2./(max(ntrials,1));
        case 'nobs'
            penalties=abs(2*ncloser-ntrials);
    end
    ties=setdiff(canflip,canflip_oneway);
    results.n_ties=length(ties);
    switch opts.method
        case 'none' %do nothing
        case 'flip_one'
        %flip all the no-ties
        which_flip=(1-tokens_oneflip)*[1:ncomps]'; %find the single zero element
        for it=canflip_oneway(:)'
            ncloser_conform(it,which_flip(it))=ntrials(it,which_flip(it))-ncloser(it,which_flip(it));
            if_flip(it,which_flip(it))=1;
        end
        %handle the ties
        penalties(tokens_oneflip==1)=Inf; %only need to look at flips that move to non-excluded area
        for it=ties(:)'
            p=find(penalties(it,:)==min(penalties(it,:)));
            if length(p)>1
                p=p(ceil(length(p)*rand(1)));
            end
            ncloser_conform(it,p)=ntrials(it,p)-ncloser(it,p);
            if_flip(it,p)=1;
        end
    end
end
opts_used=opts;
opts_used.results=results;
return

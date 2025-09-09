function [ind_cent,ind_recip,opts_used]=psg_cent_recip(dists,opts)
%[ind_cent,ind_recip,opts_used]=psg_cent_recip(dists,opts);computes centrality and reciprocity indices
% 
% Tversky, A. & Hutchinson, J. W. Nearest neighbor analysis of
% psychological spaces Psychological Review, 93(1), 3–22 (1986). 
%
%  dists: [nstims, nstims], symmetric array of distances
%  opts: 
%    opts.if_log: 1 to log results, defaults to 0eivs
%    opts.tiebreak: 'random' to break randomly (default, as in Tversky & Hutchinson)
%                   'share' to use the average rank, or to share
%    opts.seed: 1 for frozen random numbers (default), 0 for new random numbers each time, <0 for a specific seed
%    random number generator state is saved at beginning and restored at end
%
% ind_cent: cenrality index
%    1 if every vertex is nearest neighbor of exactly one vertex
%    larger if some vertices are nearest-neighbor of several vertices
%    max is (nstims^2-2*nstims+2)/nstims
%  ind_recip: reciprocity index
%    1 if nearest-neighbor relations are reciprocal
%    maximal if one vertex is nearest neighbor of all others
%    max is (nstims-1)/2+1/nstims
%  opts_used: options used, includes intermediate calcuations and maximal values of indices
%
% See also: PSG_CENT_RECIP_DEMO.
%
if (nargin<=1)
    opts=struct;
end
opts=filldefault(opts,'if_log',0);
opts=filldefault(opts,'tiebreak','random');
opts=filldefault(opts,'seed',1);
rand_state=rng;
if opts.seed==1
    rng('default');
elseif  opts.seed<0
    rng('default');
    rand(1,abs(opts.seed));
else %opts.seed=0
    rng('shuffle');
end
nstims=size(dists,1);
dists=dists+diag(Inf(nstims,1)); %make sure that a stim is not the nearest neighbor of itselt
%accumulate nearest-neighbor statistics
nrnb_tbl=zeros(1,nstims); %nearest neighbor of each stimulus
nrnb_cts=zeros(1,nstims); %counts of how many stimuli are the nearest neighbor
nrnb_list=cell(1,nstims); %list of nearest neigbbors (can be more than 1 if a tie)
nrnb_use=cell(1,nstims); %nearest neighbor to use (can be more than 1 if tiebreak='share')
nrnb_ties=zeros(1,nstims); %number of tied nearest neighbors
rank_table=zeros(nstims,nstims);
for istim=1:nstims   
    nrnb_list{istim}=find(dists(:,istim)==min(dists(:,istim)));
    nrnb_ties(istim)=length(nrnb_list{istim});
    nrnb_use{istim}=nrnb_list{istim};
    if nrnb_ties(istim)>1 & strcmp(opts.tiebreak,'random')
        nrnb_use{istim}=nrnb_list{istim}(randi(nrnb_ties(istim)));
    else
        nrnb_use{istim}=nrnb_list{istim};
    end
    nrnb_cts(nrnb_use{istim})=nrnb_cts(nrnb_use{istim})+1/length(nrnb_use{istim});
    %
    %rank the distances, breaking ties as needed
    [us,ui,rank_raw]=unique(dists(istim,:)); %rank_raw is the rank order, but may contain duplicates
    ranks_unique=zeros(1,nstims);
    rank_offset=0;
    for rank_orig=1:max(rank_raw)
        locs=find(rank_raw==rank_orig);
        nlocs=length(locs);
        switch opts.tiebreak
            case 'share'
                rank_new=rank_orig+rank_offset+(nlocs-1)/2; %use the average rank
                ranks_unique(locs)=rank_new;
            case 'random'
                ranks_new=rank_orig+rank_offset+randperm(nlocs)-1; %assign random ranks in the range rank_orig+rank_offset+[0:nlocs-1]
                ranks_unique(locs)=ranks_new;
        end
        rank_offset=rank_offset+nlocs-1;
    end
    rank_table(istim,:)=ranks_unique;
    %
end %istim
%compute indices
ind_cent=mean(nrnb_cts.^2);
rank_table_contrib=zeros(1,nstims);
for istim=1:nstims
    rank_table_contrib(istim)=mean(rank_table(nrnb_use{istim},istim));
end
ind_recip=mean(rank_table_contrib);
%restore random number generator state
rng(rand_state);
%save intermediates
opts.nrnb_list=nrnb_list;
opts.nrnb_ties=nrnb_ties;
opts.nrnb_use=nrnb_use;
opts.nrnb_cts=nrnb_cts;
opts.rank_table=rank_table;
opts.rank_table_contrib=rank_table_contrib;
opts.ind_cent_max=(nstims^2-2*nstims+2)/nstims;
opts.ind_recip_max=(nstims-1)/2+1/nstims;
opts_used=opts;
return
end

function [triad_stats,opts_used]=psg_triad_stats(sessions,opts)
% [triad_stats,opts_used]=psg_triad_stats(sessions,opts) tabulates and
% optionally displays the statistics of a session setup (how many of each stimulus used, etc)
%
% sessions: integer array [ntrials 1+ncompares nsess] of stimulus types
% opts: options field, see psg_defopts
%   opts.if_log: 1 to log
%   opts.if_cumulative: calculate statistics for sessions 1-1, 1-2, 1-3, etc. (defaults to 0)
% 
% triad_stats: structure of statistics computed
% opts_used: options used
%
% 18Nov22: add counts of number of times each triad is used
% 14Jan23: move tally display to psg_stat_tally
% 20Feb23: add opts.if_eachsess
% 
%   See also:  PSG_SESSION_STATS, PSG_SETUP_DEMO, PSG_COND_WRITE, PSG_SESSCONFIG_MAKE, PSG_DEFOPTS, PSG_UMI_STATS, PSG_STAT_TALLY.
if (nargin<2)
    opts=struct;
end
opts=psg_defopts(opts);
opts_used=opts;
ntrials=size(sessions,1);
if opts.if_log
    if isfield(opts,'cond_desc')
        disp(opts.cond_desc);
    end
    disp('triad statistics')
    disp(sprintf('%3.0f trials per session',ntrials));
end
%this also checks consistency and displays a warning if inconsistent
basic_stats=psg_session_stats(sessions,setfield(opts,'if_log',0));
triad_stats=struct;
triad_stats.warnings=[];
triad_stats.warnings=basic_stats.warnings;
%
for isess=0:opts.cond_nsess
    if (isess==0)
        ts=psg_triad_stats_do(sessions,opts);
        triad_stats.combined=ts;
        triad_count_max=length(ts.triad_repeats);
        if ts.triad_checksum~=0
            disp('warning: triad count checksum for all sessions is nonzero');
        end
    else
        if opts.if_eachsess
            ts=psg_triad_stats_do(sessions(:,:,isess),opts);
            triad_stats.eachsess.unique_triads(1,isess)=ts.unique_triads;
            triad_stats.eachsess.unique_comparison_pairs(1,isess)=ts.unique_comparison_pairs;
            triad_stats.eachsess.unique_triplets(1,isess)=ts.unique_triplets;
            triad_stats.eachsess.triad_repeats{1,isess}=ts.triad_repeats;
            triad_stats.eachsess.triads_notshown(1,isess)=ts.triads_notshown;
            if ts.triad_checksum~=0
                disp(sprintf('warning: triad count checksum for session %2.0f is nonzero',isess));
            end
        end
    end
end
if (opts.if_cumulative)
    for isess=1:opts.cond_nsess
        ts=psg_triad_stats_do(sessions(:,:,1:isess),opts);
        triad_stats.cumulative.unique_triads(1,isess)=ts.unique_triads;
        triad_stats.cumulative.unique_comparison_pairs(1,isess)=ts.unique_comparison_pairs;
        triad_stats.cumulative.unique_triplets(1,isess)=ts.unique_triplets;
        triad_stats.cumulative.triad_repeats{1,isess}=ts.triad_repeats;
        triad_stats.cumulative.triads_notshown(1,isess)=ts.triads_notshown;
        if ts.triad_checksum~=0
            disp(sprintf('warning: triad count checksum for sessions 1 - %2.0f is nonzero',isess));
        end
    end
end
if ~isempty(triad_stats.warnings)
    disp(triad_stats.warnings);
end
nstims=opts.cond_nstims;
if (opts.if_log)
    disp(' session  unique     unique    unique');
    disp('          triads   comparison unordered');
    disp('                     pairs    triplets');
    if (opts.if_eachsess)
        for isess=1:opts.cond_nsess
            disp(sprintf('   %4.0f   %6.0f   %6.0f     %6.0f',isess,...
                triad_stats.eachsess.unique_triads(1,isess),...
                triad_stats.eachsess.unique_comparison_pairs(1,isess),...
                triad_stats.eachsess.unique_triplets(1,isess)));
        end
    end
    if (opts.if_cumulative)
        for isess=1:opts.cond_nsess
            disp(sprintf('1 -%4.0f   %6.0f   %6.0f     %6.0f',isess,...
                triad_stats.cumulative.unique_triads(1,isess),...
                triad_stats.cumulative.unique_comparison_pairs(1,isess),...
                triad_stats.cumulative.unique_triplets(1,isess)));
        end
    end
    disp(sprintf('%8s  %6.0f   %6.0f     %6.0f','pooled',...
        triad_stats.combined.unique_triads,...
        triad_stats.combined.unique_comparison_pairs,...
        triad_stats.combined.unique_triplets));
    disp(sprintf('%8s  %6.0f   %6.0f     %6.0f','max poss',...
        nstims*nchoosek(nstims-1,2),nchoosek(nstims,2),nchoosek(nstims,3)));
    %
    ntriads_max=nstims*nchoosek(nstims-1,2);
    s=triad_stats;
    s.combined.triad_count_max=triad_count_max;
    psg_stat_tally('number of triads with a given number of repeats',opts,s,'triads_notshown','triad_repeats','triad_count_max');
end
return
%
function ts=psg_triad_stats_do(sessions,opts)
ts=[];
nsess=size(sessions,3);
ncompares=size(sessions,2)-1;
ntrials=size(sessions,1);
nstims=opts.cond_nstims;
ntriads_max=nstims*nchoosek(nstims-1,2); %max possible number of triads
%
triads_per_trial=nchoosek(ncompares,2);
ntriads_totshown=ntrials*triads_per_trial;
triad_list=zeros(ntriads_totshown,5); %columns are: reference, compare 1, compare 2, trial, session
%
for itrial=1:ntrials
    for isess=1:nsess
        triad_ptrs=(isess-1)*ntrials*triads_per_trial+(itrial-1)*triads_per_trial+[1:triads_per_trial];
        triad_list(triad_ptrs,1)=sessions(itrial,1,isess);
        triad_list(triad_ptrs,[2:3])=nchoosek(sort(sessions(itrial,1+[1:ncompares],isess)),2);
        triad_list(triad_ptrs,[4:5])=repmat([itrial isess],triads_per_trial,1);
    end
end
unique_triads=unique(triad_list(:,[1:3]),'rows');
ts.unique_triads=size(unique_triads,1); %unique triads
ts.unique_comparison_pairs=size(unique(unique_triads(:,[2:3]),'rows'),1); %unique pairs of comparison stimuli
ts.unique_triplets=size(unique(sort(unique_triads,2),'rows'),1); %triplets regardless of order
%tally the number of triads seen once, twice, etc.
compares_lex=zeros(ntriads_totshown*nsess,1);
for triad_ptr=1:ntriads_totshown*nsess
    compares_lex(triad_ptr)=nchoosek2seq(triad_list(triad_ptr,[2:3])); %map each pair of stimulus numbers to a unique integer
end
triads_lex=nstims*(compares_lex-1)+triad_list(:,1); %map each (reference and comparison pair) to a unique integer
counts_sparse=sparse(ones(ntriads_totshown*nsess,1),triads_lex,ones(ntriads_totshown*nsess,1));
counts=full(counts_sparse);
counts_max=max(counts);
ts.triad_repeats=zeros(1,counts_max);
for icount=1:counts_max
    ts.triad_repeats(icount)=sum(counts==icount);
end
ts.triads_notshown=ntriads_max-sum(ts.triad_repeats); %number of triads not used
ts.triad_checksum=ntriads_totshown*nsess-sum([1:counts_max].*ts.triad_repeats);
return

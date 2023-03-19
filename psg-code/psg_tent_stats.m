function [tent_stats,opts_used]=psg_tent_stats(sessions,opts)
% [tent_stats,opts_used]=psg_tent_stats(sessions,opts) tabulates and
% optionally displays the statistics of a session setup (how many of each stimulus used, etc)
% of "tent" configurations, that are relevant to testing a sufficient condition for an
% addree model.  See psg_umi_notes.docx for details.
%
% A "tent" is a set of trials with two components:
%   A triplet (or triangle): 
%     trials with s1 as ref and s2,s3 as comparisons, counted in umi_onesided_count(s1,s2,s3)
%     trials with s2 as ref and s3,s1 as comparisons, counted in umi_onesided_count(s2,s3,s1)
%     trials with s3 as ref and s1,s2 as comparisons, counted in umi_onesided_count(s3,s1,s2)
%     [These are tabulated by psg_umi_stats, in umi_onesided_count and umi_triplet_table]
%   A "tripod", a trial with s4 as ref and all of s1,s2,s3 as comparisons
%
%  The necesary condition is that if D(s4,s1), D(s4,s2), and D(s4,s3) have
%  a particular rank order, then D(s2,s3), D(s3,s1), and D(s1,s2) cannot
%  have the same rank order
%
% sessions: integer array [ntrials 1+ncompares nsess] of stimulus types
% opts: options field, see psg_defopts
%   opts.if_log: 1 to log
%   opts.if_cumulative: calculate statistics for sessions 1-1, 1-2, 1-3, etc. (defaults to 0)
%   opts.if_eachsess: caclulcate statistics for each session (defaults to 1)
%   Note that opts is modified to agree with sessions at cond_nstims, cond_ncompares, and cond_nsess.
% 
% tent_stats: structure of statistics computed, with substructures for individual trials and cumulative trials
%    tent_triplets (from umi_triplets) total number of three-legged tests of the ultrametric inequality
%    tent_triplet_table: (inherited from umi_triplet_table)
%       each row is s1,s2,s3, then umi_onesided_count(s1,s2,s3), umi_onesided_count(s2,s3,s1), umi_onesided_count(s3,s1,s2)
%    tent_table: each row is (s4,s1,s2,s3), then number of trials in tripod and then counts from tent_triplet_table
%    tents: total number of tents
%
% opts_used: options used
%
%   See also:  PSG_SESSION_STATS, PSG_SETUP_DEMO, PSG_COND_WRITE,
%   PSG_SESSCONFIG_MAKE, PSG_DEFOPTS, PSG_TRIAD_STATS, PSG_UMI_STATS, PSG_STAT_TALLY, PSG_QUAD_STATS.
%
opts_prelim=psg_defopts(sessions);
if (nargin<2)
    opts=psg_defopts(sessions);
end
opts=psg_defopts(opts);
opts.cond_nstims=opts_prelim.cond_nstims;
nstims=opts_prelim.cond_nstims;
opts.cond_ncompares=opts_prelim.cond_ncompares;
ncompares=opts_prelim.cond_ncompares;
opts.cond_nsess=opts_prelim.cond_nsess;
nsess=opts_prelim.cond_nstims;
%
ntrials=size(sessions,1);
judgs_per_trial=2*nchoosek(ncompares,2); %each pairwise judgment can be used twice
%
if opts.if_log
    if isfield(opts,'cond_desc')
        disp(opts.cond_desc);
    end
    disp('tent statistics')
    disp(sprintf('%3.0f sessions',size(sessions,3)));
    disp(sprintf('%3.0f trials per session',ntrials));
    disp(sprintf('%3.0f stimuli found',nstims));  
end
%this also checks consistency and displays a warning if inconsistent
basic_stats=psg_session_stats(sessions,setfield(opts,'if_log',0));
tent_stats=struct;
tent_stats.warnings=[];
tent_stats.warnings=basic_stats.warnings;
%
for isess=0:opts.cond_nsess
    if (isess==0)
        ts=psg_tent_stats_do(sessions,opts);
        tent_stats.combined=ts;
    else
        if opts.if_eachsess
            ts=psg_tent_stats_do(sessions(:,:,isess),opts);
            tent_stats.eachsess.tent_triplets(1,isess)=ts.tent_triplets;
            tent_stats.eachsess.tents(1,isess)=ts.tents;
            tent_stats.eachsess.tents_utriplets(1,isess)=ts.tents_utriplets;
            tent_stats.eachsess.tent_notshown(1,isess)=ts.tent_notshown;
            tent_stats.eachsess.tent_least{1,isess}=ts.tent_least;
            tent_stats.eachsess.tent_most{1,isess}=ts.tent_most;
            tent_stats.eachsess.tent_total{1,isess}=ts.tent_total;
        end
    end
end
if (opts.if_cumulative)
    for isess=1:opts.cond_nsess
        ts=psg_tent_stats_do(sessions(:,:,1:isess),opts);
        tent_stats.cumulative.tent_triplets(1,isess)=ts.tent_triplets;
        tent_stats.cumulative.tents(1,isess)=ts.tents;
        tent_stats.cumulative.tents_utriplets(1,isess)=ts.tents_utriplets;
        tent_stats.cumulative.tent_notshown(1,isess)=ts.tent_notshown;
        tent_stats.cumulative.tent_least{1,isess}=ts.tent_least;
        tent_stats.cumulative.tent_most{1,isess}=ts.tent_most;
        tent_stats.cumulative.tent_total{1,isess}=ts.tent_total;
    end
end
if ~isempty(tent_stats.warnings)
    disp(tent_stats.warnings);
end
nstims=opts.cond_nstims;
if (opts.if_log)
    disp(' summary of tent statistics');
    disp('                                 unique');
    disp(' session  number of  number of  triplets');
    disp('          triplets     tents    in tents');
    if (opts.if_eachsess)
        for isess=1:opts.cond_nsess
            disp(sprintf('   %4.0f   %6.0f     %6.0f     %6.0f',...
                isess,tent_stats.eachsess.tent_triplets(1,isess),tent_stats.eachsess.tents(1,isess),...
                tent_stats.eachsess.tents_utriplets(1,isess)));
        end
    end
    disp(sprintf('available %6.0f     %6.0f     %6.0f',...
       nchoosek(nstims,3),nstims*nchoosek(nstims-1,3),nchoosek(nstims,3)));
    if (opts.if_cumulative)
        for isess=1:opts.cond_nsess
            disp(sprintf('1 -%4.0f   %6.0f     %6.0f     %6.0f',isess,...
                tent_stats.cumulative.tent_triplets(1,isess),tent_stats.cumulative.tents(1,isess),...
                tent_stats.cumulative.tents_utriplets(1,isess)));
        end
    end
    disp(sprintf('combined  %6.0f     %6.0f     %6.0f',...
        tent_stats.combined.tent_triplets,tent_stats.combined.tents,tent_stats.combined.tents_utriplets));
    %
    disp(sprintf('available %6.0f     %6.0f     %6.0f',...
       nchoosek(nstims,3),nstims*nchoosek(nstims-1,3),nchoosek(nstims,3)));
    %
    %tallies
    %
    psg_stat_tally('number of tents with at least the given number of repeats per arm',opts,tent_stats,'tent_notshown','tent_least','tent_least_max');
    psg_stat_tally('number of tents with at most  the given number of repeats per arm',opts,tent_stats,'tent_notshown','tent_most', 'tent_most_max');
    psg_stat_tally('number of tents with a total of the given number of repeats'      ,opts,tent_stats,'tent_notshown','tent_total','tent_total_max');
end
opts_used=opts;
return
%
function ts=psg_tent_stats_do(sessions,opts)
ts=[];
%concatenate all sessions
sessions_concat=zeros(0,size(sessions,2));
nsess=size(sessions,3);
sessions_concat=sessions(:,:,1);
if (nsess>1)
    for isess=2:nsess
        sessions_concat=[sessions_concat;sessions(:,:,isess)];
    end
end
ncompares=size(sessions,2)-1;
ntrials_concat=size(sessions_concat,1);
nstims=opts.cond_nstims;
%use umi to get triplets
opts_umi=[];
opts_umi.if_log=0;
opts_umi.if_eachsess=0;
opts_umi.if_cumulative=0;
umi=getfield(psg_umi_stats(sessions,opts_umi),'combined');
ts.tent_triplets=umi.umi_triplets;
ts.tent_triplet_table=umi.umi_triplet_table;
%
trials_cmp=cell(1,nstims); %list of all trials that contain a given stimulus as a comparison
for istim=1:nstims
    trials_cmp{istim}=find(any(sessions_concat(:,1+[1:ncompares])==istim,2));
end
tent_table=zeros(0,8); %c1: reference, c2-4: comparisons,  c5: number of trials with c1 as ref, c6-8: number of trials for the triplet
tents_utriplets=0;
for itriplet=1:ts.tent_triplets
    triplet=ts.tent_triplet_table(itriplet,[1:3]);
    ref_trials=intersect(intersect(trials_cmp{triplet(1)},trials_cmp{triplet(2)}),trials_cmp{triplet(3)});
    if ~isempty(ref_trials)
        refs=sessions_concat(ref_trials,1); %the reference stimuli
        urefs=unique(refs);
        nu=length(urefs);
        tents_found=zeros(nu,8);
        tents_found(:,[2:4 6:8])=repmat(ts.tent_triplet_table(itriplet,:),nu,1); %tents that share the same triplet
        for iu=1:length(urefs)
            tents_found(iu,1)=urefs(iu);
            tents_found(iu,5)=sum(refs==urefs(iu));
        end
        tent_table=[tent_table;tents_found];
        tents_utriplets=tents_utriplets+1;
    end
end
%
ts.tent_table=tent_table;
ts.tents=size(tent_table,1);
ts.tents_utriplets=tents_utriplets;
%
%how many tents have at least a given number of trials for each arm
count_cols=[5:8]; %trial counts are in these columns
counts_least=min(tent_table(:,count_cols),[],2); 
tent_least_max=max(counts_least);
ts.tent_least_max=tent_least_max;
ts.tent_least=zeros(1,tent_least_max);
for icount=1:tent_least_max
    ts.tent_least(icount)=sum(counts_least==icount);
end
ts.tent_notshown=nstims*nchoosek(nstims-1,3)-sum(ts.tent_least);
%
%how many tents have at most a given number of trials for each arm
counts_most=max(tent_table(:,count_cols),[],2); 
tent_most_max=max(counts_most);
ts.tent_most_max=tent_most_max;
ts.tent_most=zeros(1,tent_most_max);
for icount=1:tent_most_max
    ts.tent_most(icount)=sum(counts_most==icount);
end
%
%how many tents have a total given number of trials for each arm
counts_total=sum(tent_table(:,count_cols),2); 
tent_total_max=max(counts_total);
ts.tent_total_max=tent_total_max;
ts.tent_total=zeros(1,tent_total_max);
for icount=1:tent_total_max
    ts.tent_total(icount)=sum(counts_total==icount);
end
%
return

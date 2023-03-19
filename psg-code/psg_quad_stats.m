function [quad_stats,opts_used]=psg_quad_stats(sessions,opts)
% [umi_stats,opts_used]=psg_umi_stats(sessions,opts) tabulates and
% optionally displays the statistics of a session setup (how many of each stimulus used, etc)
% that are relevant to testing quadruples of comparisons, related to the four-point
% condition for addtrees
%
% A quadruplet of stimuli a,b,c,d can be used to test additivity if
%   there are trials that can estimate Q(a,b;c,d) and Q(c,d;b,a), where
%   Q(a,b;c,d)=R(a;c,b)-R(b;d,a) and  Q(c,d;a,b)=R(c;a,d)-R(d;b,c),
%   and R(s1;s2,s3) is the probability of trials in which, with s1 as a reference
%   that s2 is chosen to be closer to s1 than s3.
%   Q(a,b;c,d) and Q(c,d;b,a) requires quad_cycle_table(a,b,d,c) to be nonzero.
%
% Note that Q(a,b;c,d)=-Q(c,d;a,b).
% See psg_umi_notes.docx and psg_umi_notes.pptx.
%
% sessions: integer array [ntrials 1+ncompares nsess] of stimulus types
% opts: options field, see psg_defopts
%   opts.if_log: 1 to log
%   opts.if_cumulative: calculate statistics for sessions 1-1, 1-2, 1-3, etc. (defaults to 0)
%   Note that opts is modified to agree with sessions at cond_nstims, cond_ncompares, and cond_nsess.
% 
%quad_stats: structure of statistics computed, with substructures for individual trials and cumulative trials
%    umi_onesided_count(s1,s2,s3) is the number of trials for which s1 is a reference and s2 and s3 are  
%      comparison stimuli. (So it is symmetric in s2 and s3)
%    quad_cycle_table: [nstims nstims nstims nstims 4], quad_cycle_table(s1,s2,s3,s4,:) is umi_onesided_count(s1,s4,s2),(s2,s1,s3),(s3,s2,s4),(s4,s3,s1)
%    quad_cycle_count: [nstims nstims nstims nstims]: quad_cycle_count(s1,s2,s3,s4) is minimum of quad_cycle_table(s1,s2,s3,s4,:)
%    quad_cycle_any_total: total number of unique entries of quad_cycle_count that are >0 (quad_cycle_count has every cycle 4 times)
%    quad_cycle_[least,most,total] is the number of cycles that have [at least, at most, total] number of trials
%       The maximum number of cycles is 6*nchoosek(nstims,4).
%    quad_loop* are analogous quantities for loops, i.e., cycles independent of direction.  
%       Since quad_cycle_table is independent of the direction around the loop, every cycle is paired with an opposite-direction cycle,
%       and there is one loop for each such pair.
%       The maximum number of loops is 3*nchoosek(nstims,4).
%    quad_trip* are the analogous quantities for triplets of loops.  These are the three inequivalent ways (a,b,c,d); (a,c,d,b); (a,d,b,c) of making loops out of four points.
%       The maximum number of loop triplets is nchoosek(nstims,4).
%
%  20Feb23: add opts.if_eachsess
%
% opts_used: options used
%   See also:  PSG_SESSION_STATS, PSG_SETUP_DEMO, PSG_COND_WRITE, PSG_QUAD_STATS_DEMO,
%   PSG_SESSCONFIG_MAKE, PSG_DEFOPTS, PSG_TRIAD_STATS, PSG_UMI_STATS, PSG_STAT_TALLY.
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
    disp('quadruplet statistics')
    disp(sprintf('%3.0f sessions',size(sessions,3)));
    disp(sprintf('%3.0f trials per session',ntrials));
    disp(sprintf('%3.0f stimuli found',nstims));  
end
%this also checks consistency and displays a warning if inconsistent
basic_stats=psg_session_stats(sessions,setfield(opts,'if_log',0));
quad_stats=struct;
quad_stats.warnings=[];
quad_stats.warnings=basic_stats.warnings;
%
for isess=0:opts.cond_nsess
    if (isess==0)
        ts=psg_quad_stats_do(sessions,opts);
        quad_stats.combined=ts;
    else
        if opts.if_eachsess
            ts=psg_quad_stats_do(sessions(:,:,isess),opts);
            quad_stats.eachsess.quad_cycle_count{1,isess}=ts.quad_cycle_count;
            quad_stats.eachsess.quad_cycle_any_total(1,isess)=ts.quad_cycle_any_total;
            quad_stats.eachsess.quad_cycle_notshown(1,isess)=ts.quad_cycle_notshown;
            quad_stats.eachsess.quad_cycle_least{1,isess}=ts.quad_cycle_least;
            quad_stats.eachsess.quad_cycle_most{1,isess}=ts.quad_cycle_most;
            quad_stats.eachsess.quad_cycle_total{1,isess}=ts.quad_cycle_total;
            %
            quad_stats.eachsess.quad_loop_notshown(1,isess)=ts.quad_loop_notshown;
            quad_stats.eachsess.quad_loop_least{1,isess}=ts.quad_loop_least;
            quad_stats.eachsess.quad_loop_most{1,isess}=ts.quad_loop_most;
            quad_stats.eachsess.quad_loop_total{1,isess}=ts.quad_loop_total;
            %
            quad_stats.eachsess.quad_trip_any_total(1,isess)=ts.quad_trip_any_total;
            quad_stats.eachsess.quad_trip_notshown(1,isess)=ts.quad_trip_notshown;
            quad_stats.eachsess.quad_trip_least{1,isess}=ts.quad_trip_least;
            quad_stats.eachsess.quad_trip_most{1,isess}=ts.quad_trip_most;
            quad_stats.eachsess.quad_trip_total{1,isess}=ts.quad_trip_total;
        end
    end
end
if (opts.if_cumulative)
    for isess=1:opts.cond_nsess
        ts=psg_quad_stats_do(sessions(:,:,1:isess),opts);
        quad_stats.cumulative.quad_cycle_count{1,isess}=ts.quad_cycle_count;
        quad_stats.cumulative.quad_cycle_any_total(1,isess)=ts.quad_cycle_any_total;
        quad_stats.cumulative.quad_cycle_notshown(1,isess)=ts.quad_cycle_notshown;
        quad_stats.cumulative.quad_cycle_least{1,isess}=ts.quad_cycle_least;
        quad_stats.cumulative.quad_cycle_most{1,isess}=ts.quad_cycle_most;
        quad_stats.cumulative.quad_cycle_total{1,isess}=ts.quad_cycle_total;
        %
        quad_stats.cumulative.quad_loop_notshown(1,isess)=ts.quad_loop_notshown;
        quad_stats.cumulative.quad_loop_least{1,isess}=ts.quad_loop_least;
        quad_stats.cumulative.quad_loop_most{1,isess}=ts.quad_loop_most;
        quad_stats.cumulative.quad_loop_total{1,isess}=ts.quad_loop_total;
        %
        quad_stats.cumulative.quad_trip_any_total(1,isess)=ts.quad_trip_any_total;
        quad_stats.cumulative.quad_trip_notshown(1,isess)=ts.quad_trip_notshown;
        quad_stats.cumulative.quad_trip_least{1,isess}=ts.quad_trip_least;
        quad_stats.cumulative.quad_trip_most{1,isess}=ts.quad_trip_most;
        quad_stats.cumulative.quad_trip_total{1,isess}=ts.quad_trip_total;
     end
end
if ~isempty(quad_stats.warnings)
    disp(quad_stats.warnings);
end
nstims=opts.cond_nstims;
if (opts.if_log)
    disp(' summary of tests of quadruples of comparisons');
    disp(' session  number of  number of    number of');
    disp('           quad        quad      triplets of');
    disp('          cycles      loops      quad loops');
     if (opts.if_eachsess)
         for isess=1:opts.cond_nsess
            disp(sprintf('   %4.0f   %6.0f     %6.0f     %6.0f     %6.0f',isess,...
                 quad_stats.eachsess.quad_cycle_any_total(1,isess),...
                quad_stats.eachsess.quad_cycle_any_total(1,isess)/2,...
                quad_stats.eachsess.quad_trip_any_total(1,isess)));
        end
     end
     disp(sprintf('available %6.0f     %6.0f     %6.0f     %6.0f',...
        6*nchoosek(nstims,4),3*nchoosek(nstims,4),nchoosek(nstims,4)));
     if (opts.if_cumulative)
         for isess=1:opts.cond_nsess
             disp(sprintf('1 -%4.0f   %6.0f     %6.0f     %6.0f     %6.0f',isess,...
                quad_stats.cumulative.quad_cycle_any_total(1,isess),...
                quad_stats.cumulative.quad_cycle_any_total(1,isess)/2,...
                quad_stats.cumulative.quad_trip_any_total(1,isess)))
         end
     end
     disp(sprintf('combined  %6.0f     %6.0f     %6.0f     %6.0f',...
         quad_stats.combined.quad_cycle_any_total,...
         quad_stats.combined.quad_cycle_any_total/2,...
         quad_stats.combined.quad_trip_any_total));
     disp(sprintf('available %6.0f     %6.0f     %6.0f     %6.0f',...
        6*nchoosek(nstims,4),3*nchoosek(nstims,4),nchoosek(nstims,4)));
    %
    %tallies
    %
    psg_stat_tally('number of quad cycles with at least the given number of repeats per arm',opts,quad_stats,'quad_cycle_notshown','quad_cycle_least','quad_cycle_least_max');
    psg_stat_tally('number of quad cycles with at most  the given number of repeats per arm',opts,quad_stats,'quad_cycle_notshown','quad_cycle_most', 'quad_cycle_most_max');
    psg_stat_tally('number of quad cycles with a total of the given number of repeats'      ,opts,quad_stats,'quad_cycle_notshown','quad_cycle_total','quad_cycle_total_max');
    %
    psg_stat_tally('number of quad loops with at least the given number of repeats per arm',opts,quad_stats,'quad_loop_notshown','quad_loop_least','quad_loop_least_max');
    psg_stat_tally('number of quad loops with at most  the given number of repeats per arm',opts,quad_stats,'quad_loop_notshown','quad_loop_most', 'quad_loop_most_max');
    psg_stat_tally('number of quad loops with a total of the given number of repeats'      ,opts,quad_stats,'quad_loop_notshown','quad_loop_total','quad_loop_total_max');
    %
    psg_stat_tally('number of quad triplets with at least the given number of repeats per arm',opts,quad_stats,'quad_trip_notshown','quad_trip_least','quad_trip_least_max');
    psg_stat_tally('number of quad triplets with at most  the given number of repeats per arm',opts,quad_stats,'quad_trip_notshown','quad_trip_most', 'quad_trip_most_max');
    psg_stat_tally('number of quad triplets with a total of the given number of repeats'      ,opts,quad_stats,'quad_trip_notshown','quad_trip_total','quad_trip_total_max');

end
opts_used=opts;
return
%
function ts=psg_quad_stats_do(sessions,opts)
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
%
trials_ref=cell(1,nstims);
trials_cmp=cell(1,nstims);
for istim=1:nstims
    trials_ref{istim}=find(sessions_concat(:,1)==istim);
    trials_cmp{istim}=find(any(sessions_concat(:,1+[1:ncompares])==istim,2));
end
% umi_onesided_count(s1,s2,s3) is the number of trials for which s1 is a reference and s2 and s3 are  
%   comparison stimuli. (So it is symmetric in s2 and s3)
umi_onesided_count=zeros(nstims,nstims,nstims);
umi_onesided_triallist=cell(nstims,nstims,nstims);
for icmp=1:nstims-1
    for jcmp=icmp+1:nstims
        ijcmp=intersect(trials_cmp{icmp},trials_cmp{jcmp});
        if ~isempty(ijcmp)
            for iref=1:nstims
                count_list=intersect(trials_ref{iref},ijcmp);
                umi_onesided_triallist{iref,icmp,jcmp}=count_list';
                umi_onesided_count(iref,icmp,jcmp)=length(count_list);
                umi_onesided_count(iref,jcmp,icmp)=length(count_list);               
            end %ifref
        end %ijcmp not empty
    end %jcmp
end %icmp
quad_cycle_table=zeros([nstims nstims nstims nstims 4]);
%    quad_cycle_table: [nstims nstims nstims nstims 4], quad_cycle_table(s1,s2,s3,s4,:) is umi_onesided_count(s1,s4,s2),(s2,s1,s3),(s3,s2,s4),(s4,s3,s1)
for ia=1:nstims
    for ib=1:nstims
        for ic=1:nstims
            u=umi_onesided_count(ib,ia,ic);
            quad_cycle_table(ib,ic,:,ia,1)=u;
            quad_cycle_table(ia,ib,ic,:,2)=u;
            quad_cycle_table(:,ia,ib,ic,3)=u;
            quad_cycle_table(ic,:,ia,ib,4)=u;
        end
    end
end
quad_cycle_count=zeros(nstims,nstims,nstims,nstims); %first index is always lowest, indices 2,3,4 can be in any order
for ia=1:nstims-3
    for ibr=ia+1:nstims-2
        for icr=ibr+1:nstims-1
            for idr=icr+1:nstims
                ibcd=perms([ibr icr,idr]);
                for iperm=1:size(ibcd,1)
                    ib=ibcd(iperm,1);
                    ic=ibcd(iperm,2);
                    id=ibcd(iperm,3);
                    u=min(quad_cycle_table(ia,ib,ic,id,:));
                    quad_cycle_count(ia,ib,ic,id)=u;
                    quad_cycle_count(ib,ic,id,ia)=u;
                    quad_cycle_count(ic,id,ia,ib)=u;
                    quad_cycle_count(id,ia,ib,ic)=u;
                end
            end
        end
    end
end
%
%summarize quad cycles and quad loops
%
ts.umi_onesided_count=umi_onesided_count;
ts.quad_cycle_table=quad_cycle_table;
quad_cycle_table_reshaped=reshape(quad_cycle_table,[nstims^4,4]);
ts.quad_cycle_any_total=sum(min(quad_cycle_table_reshaped,[],2)>0)/4;% each cycle counted four times
ts.quad_cycle_count=quad_cycle_count;
%
%how many quad cycles have at least a given number of trials for each arm
quad_cycle_least=min(quad_cycle_table_reshaped,[],2);
quad_cycle_least_max=max(quad_cycle_least);
ts.quad_cycle_least_max=quad_cycle_least_max;
ts.quad_cycle_least=zeros(1,quad_cycle_least_max);
for icount=1:quad_cycle_least_max
    ts.quad_cycle_least(icount)=sum(quad_cycle_least==icount)/4; %each cycle counted four times
end
ts.quad_cycle_notshown=6*nchoosek(nstims,4)-sum(ts.quad_cycle_least);
%
%how many quad cycles have at most a given number of trials for each arm
quad_cycle_most=max(quad_cycle_table_reshaped,[],2);
quad_cycle_most_max=max(quad_cycle_most);
ts.quad_cycle_most_max=quad_cycle_most_max;
ts.quad_cycle_most=zeros(1,quad_cycle_most_max);
for icount=1:quad_cycle_most_max
    ts.quad_cycle_most(icount)=sum(quad_cycle_most==icount)/4; %each cycle counted four times
end
%
%how many quad cycles have a total given number of trials for each arm
quad_cycle_total=sum(quad_cycle_table_reshaped,2);
quad_cycle_total_max=max(quad_cycle_total);
ts.quad_cycle_total_max=quad_cycle_total_max;
ts.quad_cycle_total=zeros(1,quad_cycle_total_max);
for icount=1:quad_cycle_total_max
    ts.quad_cycle_total(icount)=sum(quad_cycle_total==icount)/4; %each cycle counted four times
end
%
%cycles are paired into loops
ts.quad_loop_least_max=quad_cycle_least_max;
ts.quad_loop_most_max=quad_cycle_most_max;
ts.quad_loop_total_max=quad_cycle_total_max;
ts.quad_loop_notshown=ts.quad_cycle_notshown/2;
ts.quad_loop_least=ts.quad_cycle_least/2;
ts.quad_loop_most=ts.quad_cycle_most/2;
ts.quad_loop_total=ts.quad_cycle_total/2;
%
%look for loop triplets
%
quad_trip_table=zeros([nstims nstims nstims nstims 12]);
for ia=1:nstims-3
    for ib=ia+1:nstims-2
        for ic=ib+1:nstims-1
            for id=ic+1:nstims
                quad_trip_table(ia,ib,ic,id,:)=reshape([...
                    squeeze(quad_cycle_table(ia,ib,ic,id,:));...
                    squeeze(quad_cycle_table(ia,ic,id,ib,:));...
                    squeeze(quad_cycle_table(ia,id,ib,ic,:))],[1 1 1 1 12]);
            end
        end
    end
end
%
%summarize loop triplets
%
ts.quad_trip_table=quad_trip_table;
quad_trip_table_reshaped=reshape(quad_trip_table,[nstims^4,12]);
ts.quad_trip_any_total=sum(min(quad_trip_table_reshaped,[],2)>0);
%
%how many loop triplets have at least a given number of trials for each arm
quad_trip_least=min(quad_trip_table_reshaped,[],2);
quad_trip_least_max=max(quad_trip_least);
ts.quad_trip_least_max=quad_trip_least_max;
ts.quad_trip_least=zeros(1,quad_trip_least_max);
for icount=1:quad_trip_least_max
    ts.quad_trip_least(icount)=sum(quad_trip_least==icount);
end
ts.quad_trip_any_total=sum(ts.quad_trip_least);
ts.quad_trip_notshown=nchoosek(nstims,4)-sum(ts.quad_trip_least);
%
%how many loop triplets have at most a given number of trials for each arm
quad_trip_most=max(quad_trip_table_reshaped,[],2);
quad_trip_most_max=max(quad_trip_most);
ts.quad_trip_most_max=quad_trip_most_max;
ts.quad_trip_most=zeros(1,quad_trip_most_max);
for icount=1:quad_trip_most_max
    ts.quad_trip_most(icount)=sum(quad_trip_most==icount);
end
%
%how many loop triplets have a total given number of trials for each arm
quad_trip_total=sum(quad_trip_table_reshaped,2);
quad_trip_total_max=max(quad_trip_total);
ts.quad_trip_total_max=quad_trip_total_max;
ts.quad_trip_total=zeros(1,quad_trip_total_max);
for icount=1:quad_trip_total_max
    ts.quad_trip_total(icount)=sum(quad_trip_total==icount);
end
%
return

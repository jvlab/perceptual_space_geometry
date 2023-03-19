function [umi_stats,opts_used]=psg_umi_stats(sessions,opts)
% [umi_stats,opts_used]=psg_session_stats(sessions,opts) tabulates and
% optionally displays the statistics of a session setup (how many of each stimulus used, etc)
% that are relevant to ttesting the ultrametric inequality
%
% An ultrametric inequality pair umipair is a pair of stimuli A and B for which
%   there ia a trial with A as the standard and B and C as the reference,
%   and also a trial with B as the standard and A and C as a reference.
% The ultrametric inequality is that d(A,B)<=max(d(A,C),d(B,C)), so that,
%   if both trials are present, at least one of d(A,B)<=d(A,C) and d(A,B)<=d(B,C) must hold.
%
% sessions: integer array [ntrials 1+ncompares nsess] of stimulus types
% opts: options field, see psg_defopts
%   opts.if_log: 1 to log
%   opts.if_cumulative: calculate statistics for sessions 1-1, 1-2, 1-3, etc. (defaults to 0)
%   Note that opts is moditifed to agree with sessions at cond_nstims, cond_ncompares, and cond_nsess.
% 
% umi_stats: structure of statistics computed, with substructures for
%    individual trials and cumulative trials
%    umi_onesided_count(s1,s2,s3) is the number of trials for which s1 is a reference and s2 and s3 are  
%      comparison stimuli. (So it is symmetric in s2 and s3)
%    umi_count(s1,s2,s3) is umi_onesided_count(s1,s2,s3) but zero if umi_onesided_count(s2,s1,s3) is zero.  
%    umi_any(s1,s2) is the total of umi_count(s1,s2,*) and umi_count(s2,s1,*),
%      i.e., number of configs of stimuli that participate in tests of ultrametric inequality
%    umi_nc(s1,s2) is the number of distinct s3's for which umi_count(s1,s2,s3) is nonzero, 
%      i.e., number of trials that participate in tests of ultrametric inequality
%    umi_any_total and umi_nc_total are sums of umi_any and umi_nc, both divided
%      by 2, since they are symmetric and doubly count the trials and stimulus pairs
%    umi_pairs is th total number of umi_any or umi_nc that are >0
% opts_used: options used
% 
%   See also:  PSG_SESSION_STATS, PSG_SETUP_DEMO, PSG_COND_WRITE,
%   PSG_SESSCONFIG_MAKE, PSG_DEFOPTS, PSG_TRIAD_STATS, PSG_UMI_STATS_DEMO.
%
opts_prelim=psg_defopts(sessions);
if (nargin<2)
    opts=psg_defopts(sessions);
end
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
    disp('ultrametric inequality statistics')
    disp(sprintf('%3.0f sessions',size(sessions,3)));
    disp(sprintf('%3.0f trials per session',ntrials));
    disp(sprintf('%3.0f stimuli found',nstims));  
end
%this also checks consistency and displays a warning if inconsistent
basic_stats=psg_session_stats(sessions,setfield(opts,'if_log',0));
umi_stats=struct;
umi_stats.warnings=[];
umi_stats.warnings=basic_stats.warnings;
%
for isess=0:opts.cond_nsess
    if (isess==0)
        ts=psg_umi_stats_do(sessions,opts);
        umi_stats.combined=ts;
    else
        ts=psg_umi_stats_do(sessions(:,:,isess),opts);
        umi_stats.eachsess.umi_any_total(1,isess)=ts.umi_any_total;
        umi_stats.eachsess.umi_nc_total(1,isess)=ts.umi_nc_total;
        umi_stats.eachsess.umi_pairs(1,isess)=ts.umi_pairs;
    end
end
if (opts.if_cumulative)
    for isess=1:opts.cond_nsess
        ts=psg_umi_stats_do(sessions(:,:,1:isess),opts);
        umi_stats.cumulative.umi_any_total(1,isess)=ts.umi_any_total;
        umi_stats.cumulative.umi_nc_total(1,isess)=ts.umi_nc_total;
        umi_stats.cumulative.umi_pairs(1,isess)=ts.umi_pairs;

    end
end
if ~isempty(umi_stats.warnings)
    disp(umi_stats.warnings);
end
nstims=opts.cond_nstims;
if (opts.if_log)
    disp(' summary of tests of ultrametric inequality');
    disp(' session  number of  number of   number of');
    disp('          judgment   stimulus    stimulus');
    disp('           pairs      configs     pairs');
    for isess=1:opts.cond_nsess
        disp(sprintf('   %4.0f   %6.0f     %6.0f     %6.0f',isess,...
            umi_stats.eachsess.umi_any_total(1,isess),...
            umi_stats.eachsess.umi_nc_total(1,isess),...
             umi_stats.eachsess.umi_pairs(1,isess)));
    end
    disp(sprintf('available %6.0f     %6.0f     %6.0f',...
       judgs_per_trial*ntrials,nchoosek(nstims,3)*3,nchoosek(nstims,2))); %T

    if (opts.if_cumulative)
        for isess=1:opts.cond_nsess
            disp(sprintf('1 -%4.0f   %6.0f     %6.0f     %6.0f',isess,...
                umi_stats.cumulative.umi_any_total(1,isess),...
                umi_stats.cumulative.umi_nc_total(1,isess),...
                umi_stats.cumulative.umi_pairs(1,isess)));
        end
    end
    disp(sprintf('combined  %6.0f     %6.0f     %6.0f',...
        umi_stats.combined.umi_any_total,umi_stats.combined.umi_nc_total,umi_stats.combined.umi_pairs));
    %
    disp(sprintf('available %6.0f     %6.0f     %6.0f',...
       judgs_per_trial*ntrials*size(sessions,3),nchoosek(nstims,3)*3,nchoosek(nstims,2))); %T
end
return
%
function ts=psg_umi_stats_do(sessions,opts)
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
% umi_count(s1,s2,s3) is umi_onesided_count(s1,s2,s3) but zero if umi_onesided_count(s2,s1,s3) is zero.  
umi_count=zeros(nstims,nstims,nstims);
% umi_any(s1,s2) is the total of umi_count(s1,s2,*) and umi_count(s2,s1,*)
% umi_nc(s1,s2) is the number of distinct s3's for which umi_count(s1,s2,s3) is nonzero
umi_any=zeros(nstims,nstims);
umi_nc=zeros(nstims,nstims);
for iref=1:nstims-1
    for jref=iref+1:nstims
        counts_ij=umi_onesided_count(iref,jref,:);
        counts_ji=umi_onesided_count(jref,iref,:);
        umi_count(iref,jref,:)=counts_ij.*double(counts_ji>0);
        umi_count(jref,iref,:)=counts_ji.*double(counts_ij>0);
        s=sum(umi_count(iref,jref,:),3)+sum(umi_count(jref,iref,:),3);
        umi_any(iref,jref)=s;
        umi_any(jref,iref)=s;
        umi_nc(iref,jref)=sum(double(umi_count(iref,jref,:)>0),3);
        umi_nc(jref,iref)=sum(double(umi_count(iref,jref,:)>0),3);
    end
end
ts.umi_onesided_count=umi_onesided_count;
ts.umi_onesided_triallist=umi_onesided_triallist;
ts.umi_count=umi_count;
ts.umi_any=umi_any;
ts.umi_any_total=sum(umi_any(:))/2; %double-counted, symmetric
ts.umi_nc=umi_nc;
ts.umi_nc_total=sum(umi_nc(:))/2; %double_counted, symmetric
ts.umi_pairs=sum(umi_any(:)>0)/2; %double_counted, symmetric
return

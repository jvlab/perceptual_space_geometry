function [stats,opts_used]=psg_session_stats(sessions,opts)
% [stats,opts_used]=psg_session_stats(sessions,opts) tabulates and
% optionally displays the statistics of a session setup (how many of each stimulus used, etc)
%
% sessions: integer array [ntrials 1+ncompares nsess] of stimulus types
% opts: options field, see psg_defopts
%   opts.if_log: 1 to log
% 
% stats: structure of statistics computed
% opts_used: options used
% 
%   See also:  PSG_SETUP_DEMO, PSG_COND_WRITE, PSG_SESSCONFIG_MAKE, PSG_TRIAD_STATS, PSG_DEFOPTS.
if (nargin<1)
    opts=struct;
end
opts=psg_defopts(opts);
opts_used=opts;
ntrials=size(sessions,1);
if opts.if_log
    if isfield(opts,'cond_desc')
        disp(opts.cond_desc);
    end
    disp(sprintf('%3.0f trials per session',ntrials));
end
counts_reference=zeros(opts.cond_nstims,opts.cond_nsess);
counts_compares=zeros(opts.cond_nstims,opts.cond_nsess);
stats=struct;
stats.warnings=[];
%
if length(unique(sessions(:)))~=opts.cond_nstims
    wtext=sprintf('number of unique stimuli found in session array (%2.0f) does not match number of stimuli specified in cond configuration (%2.0f)',...
        length(unique(sessions(:))),opts.cond_nstims);
    warning(wtext);
    stats.warnings=strvcat(stats.warnings,wtext);
end
if size(sessions,2)~=(1+opts.cond_ncompares)
    wtext=sprintf('number of stimuli on each trial in session array (%2.0f) does not match number of stimuli on each trial found in cond configuration (%2.0f)',...
        size(sessions,2),opts.cond_ncompares+1);
    warning(wtext);
    stats.warnings=strvcat(stats.warnings,wtext);
end
if size(sessions,3)~=(opts.cond_nsess)
    wtext=sprintf('number of sessions found in session array (%2.0f) does not match number of sessions in cond configuration (%2.0f)',...
        size(sessions,3),opts.cond_nsess);
    warning(wtext);
    stats.warnings=strvcat(stats.warnings,wtext);
end
hbins=[1:opts.cond_nstims];
for isess=1:opts.cond_nsess
    counts_reference(:,isess)=hist(sessions(:,1,isess),hbins);
    counts_compares(:,isess)=hist(reshape(sessions(:,1+[1:opts.cond_ncompares],isess),[1 ntrials*opts.cond_ncompares]),hbins);
end
counts_alluses=counts_reference+counts_compares;
if opts.if_log
    disp('            stim used as  stim used as       all');
    disp('             reference     comparison     instances');
    disp(' session     min   max     min   max     min   max');
    for isess=1:opts.cond_nsess
        disp(sprintf('  %5.0f     %4.0f  %4.0f    %4.0f  %4.0f    %4.0f  %4.0f',isess,...
            min(counts_reference(:,isess)),max(counts_reference(:,isess)),...
            min(counts_compares(:,isess)),max(counts_compares(:,isess)),...
            min(counts_alluses(:,isess)),max(counts_alluses(:,isess))));
    end
    disp(sprintf('%s  %6.0f%6.0f  %6.0f%6.0f  %6.0f%6.0f','all sess',...
        min(sum(counts_reference,2)),max(sum(counts_reference,2)),...
        min(sum(counts_compares,2)),max(sum(counts_compares,2)),...
        min(sum(counts_alluses,2)),max(sum(counts_alluses,2))));
end
stats.counts_reference=counts_reference;
stats.counts_compares=counts_compares;
stats.counts_alluses=counts_alluses;
return


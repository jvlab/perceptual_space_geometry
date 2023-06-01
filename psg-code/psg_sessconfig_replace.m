function [replaced,opts_used,warnings]=psg_sessconfig_replace(sessions,opts)
% [replaced,opts_used,warnings]=psg_sessconfig_replace(sessions,opts)
% modifies arrays that define sessions for a perceptual space geometry experiment
% by replacing one or more of the stimulus tokens with lower-number values
%
% opts: options
%   opts.if_log: 1 to log (defaults to 0)
%   opts.cond_nstims: number of stimuli (typically 25 or 37)
%   opts.cond_nstims_toreplace: number of stimuli to replace
%     the stimuli from cond_nstims-cond_nstims_toreplace+1 to cond_nstims
%     are replaced, randomly by lower-number values
% sessions: dim 1 is trial, dim 2 is stimulus number (1 to 1+ncompares+), dim 3 is session
%   note that replacement is not applied to sessions_sorted, as this will conflict with the sorting
%
% replaced: simlar array to sessions, but with stimulus tokens replaced
% opts_used: options used, modified by replacement
%     opts_used.cond_nstims=opts.cond_nstims-opts.cond_nstims_toreplace
%     opts_used.cond_nstims_toreplace=0
% warnings: warning text, including checking for compatibility of sessions
%     and opts.cond_[nstims|ncompares|nsess]
%
%  See also:  PSG_SESSION_STATS, PSG_DEFOPTS, PSG_SESSCONFIG_MAKE, FACES_MPI_PSG_SETUP.
%
if (nargin<2)
    opts=[];
end
opts=psg_defopts(opts);
opts_used=opts;
nkeep=opts.cond_nstims-opts.cond_nstims_toreplace;
opts_used.cond_nstims=nkeep;
opts_used.cond_nstims_toreplace=0;
%
ncompares=opts.cond_ncompares;
nsess=size(sessions,3);
ncompares=size(sessions,2)-1;
ntrials=size(sessions,1);
replaced=sessions;
warnings=[];
if length(unique(sessions(:)))~=opts.cond_nstims
    wtext=sprintf('number of unique stimuli found in session array (%2.0f) does not match number of stimuli specified in cond configuration (%2.0f)',...
        length(unique(sessions(:))),opts.cond_nstims);
    warning(wtext);
    warnings=strvcat(warnings,wtext);
end
if size(sessions,2)~=(1+opts.cond_ncompares)
    wtext=sprintf('number of stimuli on each trial in session array (%2.0f) does not match number of stimuli on each trial found in cond configuration (%2.0f)',...
        size(sessions,2),opts.cond_ncompares+1);
    warning(wtext);
    warnings=strvcat(warnings,wtext);
end
if size(sessions,3)~=(opts.cond_nsess)
    wtext=sprintf('number of sessions found in session array (%2.0f) does not match number of sessions in cond configuration (%2.0f)',...
        size(sessions,3),opts.cond_nsess);
    warning(wtext);
    warnings=strvcat(warnings,wtext);
end
if (1+ncompares)>nkeep
    wtext='not enough unique stimulus tokens for replacement';
    warnings=strvcat(warnings,wtext);
    warning(wtext);
    return
end
%
if opts.cond_nstims_toreplace>0
    for ireplace=nkeep+1:opts.cond_nstims
        replaced(replaced==ireplace)=0;
    end
    nreplace=sum(replaced==0,2);
    for isess=1:nsess
        for itrial=1:ntrials
            if nreplace(itrial,isess)>0
                nr=nreplace(itrial,isess);
                zloc=find(replaced(itrial,:,isess)==0);
                can_use=setdiff([1:nkeep],replaced(itrial,:,isess)); %do not duplicate a stimulus
                chosen=randperm(length(can_use));
 %               disp(replaced(itrial,:,isess));
 %               disp([itrial isess nr nkeep])
                replaced(itrial,zloc,isess)=can_use(chosen([1:nr]));
            end
        end
    end
end
return

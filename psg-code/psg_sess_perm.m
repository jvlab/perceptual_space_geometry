function [sessions,sessions_sorted,opts_used,desc,warnings,perms_used]=psg_sess_perm(sessions_orig,sessions_sorted_orig,pmode,opts)
% [sessions,sessions_sorted,opts_used,warnings,perms_used]=psg_sess_perm(sessions_orig,sessions_sorted_orig,pmode,opts)
% applies permutatations to stimulus numbers to achieve a more balanced design
%
% This is typically called after psg_sessconfig_make with opts.refseq=2.
%
% sessions_orig: results of psg_sessconfig_make [ntrials nstims nsess]
% sessions_sorted_orig: as in sessions, but trials are ordered by reference stimulus
%     and the positions of the comparison stimuli are sequential
% pmode: a structure describing the permutations
%   pmode.type='primroot': (only current choice): multiplication by successive powers of a primitive root mod nstims
%   pmode.prim: the value of the primitive root
% opts: options
%   opts.if_log: 1 to log (defaults to 0)
%   opts.cond_nstims: number of stimuli (typically 25 or 37)
%   opts.cond_ncompares: number of comparison stimuli in each trial (typically 8)
%   opts.cond_novlp: number of overlaps in two trials with the same reference stimulus (typically 2)
%   opts.cond_nsess: number of session arrays to make
%   opts.refseq: 1 (default) to randomize the sequence of reference stimuli used for context
%                2 for no randomization
%
% sessions: modified session indices, same shape as sessions_orig
% sessions_sorted: sessions_sorted_orig with permutations applied
% opts_used: options used
% desc: descriptor string
% warnings: warning text
% perms_used: size=[nsess nstims]: permutations used (0-indexed)
%
%  See also:  PSG_SESSCONFIG_MAKE, PSG_SESSION_STATS, PSG_DEFOPTS, PSG_SETUP_DEMO, PSG_COND_WRITE, RANDPERM, PRIMROOT.
%
if (nargin<4)
    opts=[];
end
opts=psg_defopts(opts);
opts_used=opts;
if opts.if_log
    disp('doing the permutations')
    if isfield(opts,'cond_desc')
        disp(sprintf('starting with %s',opts.cond_desc));
    else
        opts.cond_desc=[];
    end
end
opts_used=opts;
%this also checks consistency and displays a warning if inconsistent
basic_stats=psg_session_stats(sessions_orig,setfield(opts,'if_log',0));
warnings=basic_stats.warnings;
perms_used=repmat([0:opts.cond_nstims-1],[opts.cond_nsess,1]); %null permutation
switch pmode.type
    case 'primroot'
        desc_add=sprintf('; perms:  %s, mult %2.0f',pmode.type,pmode.prim);
        pwr=1;
        for isess=2:opts.cond_nsess
            pwr=mod(pwr*pmode.prim,opts.cond_nstims);
            perms_used(isess,:)=mod([0:opts.cond_nstims-1]*pwr,opts.cond_nstims);           
        end
    otherwise
        wtext=sprintf('unknown permutation type (%s)',pmode.type);
        warning(wtext);
        warnings=strvcat(warnings,wtext);
        desc_add=[];
end
%apply the permutations
for isess=1:opts.cond_nsess
    perm=perms_used(isess,:);
    sessions(:,:,isess)=1+perm(sessions_orig(:,:,isess)); %%apply the permutation
    sessions_sorted(:,:,isess)=1+perm(sessions_sorted_orig(:,:,isess)); %%apply the permutation
end
desc=cat(2,opts.cond_desc,desc_add);
if (opts.if_log)
    disp(sprintf('yielding %s',desc));
end
return

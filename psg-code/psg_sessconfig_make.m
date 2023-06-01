function [sessions,sessions_sorted,opts_used,desc,warnings]=psg_sessconfig_make(opts)
% [sessions,sessions_sorted,opts_used,warnings]=psg_sessconfig_make(opts)
% creates arrays that define sessions for a perceptual space geometry experiment
%
% opts: options
%   opts.if_log: 1 to log (defaults to 0)
%   opts.cond_nstims: number of stimuli (typically 25 or 37)
%   opts.cond_ncompares: number of comparison stimuli in each trial (typically 8)
%   opts.cond_novlp: number of overlaps in two trials with the same reference stimulus (typically 2)
%   opts.cond_nsess: number of session arrays to make
%   opts.refseq: 1 (default) to randomize the sequence of reference stimuli used for context
%                    (separate randomization for each session and each stimulus)
%                2 for no randomization
%                3 frozen randomization, same for all sessions and stimuli
%   opts.setseq: 1 (default) each session uses independently randomized sets of stimuli (irrelevant if refseq=2)
%                2 for each session uses same sets of stimuli but locations are randomized
%                3 for each session uses same sets of stimuli, no randomization of location
%                4 for each session uses independently randomized sets of stimuli, but no randomization of location
%
% sessions: dim 1 is trial, dim 2 is stimulus number (1 to 1+ncompares), dim 3 is session 
% sessions_sorted: as in sessions, but trials are ordered by reference stimulus
%     and the positions of the comparison stimuli are sorted in the initial
%     stimulus generation procedure -- though this may be changed by
%     psg_sess_perm or psg_sessconfig_replace.
% opts_used: options used
% desc: descriptor string
% warnings: warning text
%
%  17Nov22: add setseq: options for stimulus sets across sessions
%  31Dec22: add refseq 3: frozen randomization
%  01Jun23: fix documentation about sessions_sorted, and add nstims_toreplace to desc if nonzero
%
%  See also:  PSG_SESSION_STATS, PSG_DEFOPTS, PSG_SETUP_DEMO, PSG_COND_WRITE, RANDPERM, LCM, PSG_SESS_PERM, PSG_SESSCONFIG_REPLACE.
%
if (nargin<1)
    opts=[];
end
opts=psg_defopts(opts);
opts_used=opts;
%
nstims=opts.cond_nstims;
ncompares=opts.cond_ncompares;
novlp=opts.cond_novlp;
nsess=opts.cond_nsess;
%
warnings=[];
nstep=ncompares-novlp;
if (mod(nstims-1,nstep))>0
    wtext=sprintf('number of non-overlapping stimuli (%1.0f) is not a factor of nstims-1 (%1.0f)',nstep,nstims-1);
    warning(wtext);
    warnings=strvcat(warnings,wtext);
end
ntrials_per_ref=lcm(nstims-1,nstep)/nstep;
ntrials=ntrials_per_ref*nstims;
if opts.cond_nstims_toreplace==0
    desc=sprintf('nstims=%2.0f, ncompares=%2.0f, novlp=%1.0f, nsess=%2.0f, refseq %s setseq %s',...
        nstims,ncompares,novlp,nsess,opts.refseq_labels{opts.refseq},opts.setseq_labels{opts.setseq});
else
    desc=sprintf('nstims=%2.0f (%2.0f to replace), ncompares=%2.0f, novlp=%1.0f, nsess=%2.0f, refseq %s setseq %s',...
        nstims,opts.cond_nstims_toreplace,ncompares,novlp,nsess,opts.refseq_labels{opts.refseq},opts.setseq_labels{opts.setseq});
end
if opts.if_log
    disp(sprintf('making %2.0f sessions of %4.0f trials with %4.0f stimuli (%4.0f trials per reference stimulus, %4.0f comparison stimuli per trial, %4.0f overlaps)',...
        nsess,ntrials,nstims,ntrials_per_ref,ncompares,novlp));
    disp(sprintf('descriptor: %s',desc));
end
compares=zeros(ntrials_per_ref,ncompares);
% use same sets in each session?
setseq_label=opts.setseq_labels{opts.setseq};
if ~isempty(strfind(setseq_label,'same_sets'))
    nsess_make=1;
else
    nsess_make=nsess;
end
sessions_sorted=zeros(ntrials,ncompares+1,nsess_make);
%
for isess=1:nsess_make
    for istim=1:nstims
        %comp_list=setdiff([1:nstims],istim);
        comp_list=[istim+1:nstims, 1:istim-1]; %this cycles the roles of the stimuli, might make it easier to modify to a balanced design
        switch opts.refseq
            case 1 %random
                comp_order=randperm(nstims-1);
            case 2 %in order
                comp_order=[1:nstims-1];
            case 3 %in frozen random order
                if (istim==1) & (isess==1)
                    comp_order=randperm(nstims-1);
                end
        end
        for itrial=1:ntrials_per_ref
            comp_index=(itrial-1)*nstep+[1:ncompares];
            compares(itrial,:)=comp_list(comp_order(1+mod(comp_index,nstims-1)));
        end %itrial
        trialnos=(istim-1)*ntrials_per_ref+[1:ntrials_per_ref];
        sessions_sorted(trialnos,1,isess)=istim;
        sessions_sorted(trialnos,1+[1:ncompares],isess)=compares;
    end %istim
end %isess
if (nsess_make==1)
    sessions_sorted=repmat(sessions_sorted,[1 1 nsess]);
end
sessions=sessions_sorted;
for isess=1:nsess
    %randomly permute the positions of the comparison stimuli
    if isempty(strfind(setseq_label,'same_positions'))
        for itrial=1:ntrials
            sessions(itrial,1+[1:ncompares],isess)=sessions(itrial,1+randperm(ncompares),isess);
        end
    end
    %randomly permute the trial order
    sessions(:,:,isess)=sessions(randperm(ntrials),:,isess);
end
return

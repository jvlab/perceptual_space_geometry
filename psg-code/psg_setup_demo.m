%psg_setup_demo
% demonstrate setup for perceptual space geometry
%
% file writing will need to be changed to make use of actual stimulus file
% names (and instances of them) -- here, 'type01', etc, is used
%
% to do:
%  * fix colormap of context map: black ref, white for not used, other
%  colors for stimuli tested in multiple contexts
%  * make analysis of context data a function, which includes statistical table of context data
%
%  18Nov22: add option for frozen random numbers
%  18Nov22: add setseq: options for stimulus sets across sessions
%  19Nov22: if multiplicative scheme requested, only redo the sessions if the number of sessions needs to be adjusted
%  19Nov22: allow for a list of actions (plots, stats, write files)
%  20Nov22: split off psg_setup_process
%
%   See also:  PSG_SESSCONFIG_MAKE, PSG_SESSION_STATS, PSG_TRIAD_STATS, PSG_SPOKES_SETUP, PRIMROOT, PSG_SESS_PERM, PSG_SETUP_PROCESS.
%
if ~exist('opts_psg')
    opts_psg=struct;
end
opts_psg=psg_defopts(opts_psg);
%
opts_psg=filldefault(opts_psg,'if_log',1);
if ~exist('psg_alt_setups') %alternative setups
    psg_alt_setups{1}.nstims=37;
    psg_alt_setups{1}.ncompares=8;
    psg_alt_setups{1}.novlp=2;
    psg_alt_setups{1}.nsess=5;
    %
    psg_alt_setups{2}.nstims=25;
    psg_alt_setups{2}.ncompares=8;
    psg_alt_setups{2}.novlp=2;
    psg_alt_setups{2}.nsess=10;
    %
    psg_alt_setups{3}.nstims=19;
    psg_alt_setups{3}.ncompares=8;
    psg_alt_setups{3}.novlp=2;
    psg_alt_setups{3}.nsess=16;
    %
    psg_alt_setups{4}.nstims=29;
    psg_alt_setups{4}.ncompares=10;
    psg_alt_setups{4}.novlp=3;
    psg_alt_setups{4}.nsess=8;
end
disp(sprintf('%1.0f->setup with %3.0f stimuli, %3.0f comparison stimuli per trial, overlap %3.0f; %3.0f sessions (default)',...
    0,opts_psg.cond_nstims,opts_psg.cond_ncompares,opts_psg.cond_novlp,opts_psg.cond_nsess));
for k=1:length(psg_alt_setups)
    disp(sprintf('%1.0f->setup with %3.0f stimuli, %3.0f comparison stimuli per trial, overlap %3.0f',...
        k,psg_alt_setups{k}.nstims,psg_alt_setups{k}.ncompares,psg_alt_setups{k}.novlp));
end
k=getinp('choice (negative to modify)','d',length(psg_alt_setups)*[-1 1],0);
if k~=0
    opts_psg.cond_nstims=psg_alt_setups{abs(k)}.nstims;
    opts_psg.cond_ncompares=psg_alt_setups{abs(k)}.ncompares;
    opts_psg.cond_novlp=psg_alt_setups{abs(k)}.novlp;
    opts_psg.cond_nsess=psg_alt_setups{abs(k)}.nsess;
end
if k<0
    opts_psg.cond_nstims=getinp('nstims','d',[1 1000],opts_psg.cond_nstims);
    opts_psg.cond_ncompares=getinp('ncompares','d',[1 1000],opts_psg.cond_ncompares);
    opts_psg.cond_novlp=getinp('novlp','d',[1 1000],opts_psg.cond_novlp);
end
%
opts_psg.cond_nsess=getinp('number of sessions','d',[1 10000],opts_psg.cond_nsess);
%
for k=1:length(opts_psg.refseq_labels)
    disp(sprintf('%1.0f->method for choosing stimuli in overlap: %s',k,opts_psg.refseq_labels{k}));
end
opts_psg.refseq=getinp('choice','d',[1 length(opts_psg.refseq_labels)],opts_psg.refseq);
%
for k=1:length(opts_psg.setseq_labels)
    disp(sprintf('%1.0f->method for sequencing stimulus sets across sessions: %s',k,opts_psg.setseq_labels{k}));
end
opts_psg.setseq=getinp('choice','d',[1 length(opts_psg.setseq_labels)],opts_psg.setseq);
%
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1]);
%
if (if_frozen~=0)
    rng('default');
    if (if_frozen<0)
        rand(1,abs(if_frozen));
    end
else
    rng('shuffle');
end
%
[sessions,sessions_sorted,opts_used,psg_desc]=psg_sessconfig_make(opts_psg);
opts_psg.cond_desc=psg_desc;
%
disp(sprintf('Analyzing setup with %s',opts_psg.cond_desc));
ntrials=size(sessions,1);
%
%accumulate and display statistics of the configuration
%
stats=psg_session_stats(sessions,setfield(opts_psg,'if_log',1));
%
%do we want to permute the stimulus number via a multiplicative scheme?
%
[primroots,nrp_phi]=primroot(opts_psg.cond_nstims);
if length(primroots)==0
    disp(sprintf('no primitive roots mod %3.0f, so no option for a multiplicative scheme',opts_psg.cond_nstims));
else
    disp(sprintf('primitive roots mod %3.0f:',opts_psg.cond_nstims));
    disp(primroots);
    if_mult=getinp('1 to apply a multiplicative scheme','d',[0 1]);
    if (if_mult)
        disp('available primitive roots')
        disp(primroots);
        prim=NaN;
        while (~ismember(prim,primroots))
            prim=getinp('a primitive root','d',primroots([1 end]),primroots(1));
        end
        disp(sprintf('period is %3.0f',nrp_phi));
        if (nrp_phi~=opts_psg.cond_nsess)
            if_change_nsess=getinp(sprintf('1 to change to have %3.0f sessions instead of %3.0f',nrp_phi,opts_psg.cond_nsess),'d',[0 1]);
        else
            if_change_nsess=0;
        end
        if if_change_nsess
            opts_psg.cond_nsess=nrp_phi;
            [sessions,sessions_sorted,opts_used,psg_desc]=psg_sessconfig_make(opts_psg);
            opts_psg.cond_desc=psg_desc;
            stats=psg_session_stats(sessions,setfield(opts_psg,'if_log',1)); %redo statistics
        end
        %apply permutation based on primitive roots
        pmode=struct;
        pmode.type='primroot';
        pmode.prim=prim;
        [sessions,sessions_sorted,opts_used,psg_desc,warnings,perms_used]=...
            psg_sess_perm(sessions,sessions_sorted,pmode,setfield(opts_psg,'if_log',1));
        opts_psg.cond_desc=psg_desc;
        stats=psg_session_stats(sessions,setfield(opts_psg,'if_log',1));
    end
end
%
% write files, make plots, and calculate statistics
%
psg_setup_process(sessions,sessions_sorted,opts_psg);
%

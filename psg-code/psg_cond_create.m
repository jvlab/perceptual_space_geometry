function [session_cells,perms_used,examps_used,opts_used]=psg_cond_create(sessions,typenames,opts)
% [session_cells,perms_used,examps_used,opts_used]=psg_cond_create(sessions,typenames,opts)
% creates text arrays for condition files from numerical arrays of stimulus types
%
% sessions: integer array [ntrials 1+ncompares nsess] of stimulus types
% typenames: cell array of text strings for each stimulus type
% opts: options field, see psg_defopts
%   opts.example_infix_mode: controls re-use of stimuli
%   opts.example_infix_string: separator between stimulus name and example number
%   opts.example_infix_zpad: number of digits to zero-pad 
%
% session_cells: session_cells{isess}(itrial,istim} is a cell array of cell arrays
%     that can be used for the csv cond file
% perms_used: permutations used, zero-indexed (computed even if not needed)
%   perms_used.allsess{istim}: across all sessions, needed for unique stimuli across all trials in all sessions
%   perms_used.eachsess{istim,isess}: across each session, needed for unique stimuli within a session
% examps_used: examples used, zero-indexed. size is [ntrials 1+ncompares nsess]
% opts_used: options used
%
%   See also:  PSG_SETUP_DEMO, PSG_COND_WRITE, PSG_SESSCONFIG_MAKE, ZPAD, PSG_SESSION_STATS.

if (nargin<3)
    opts=struct;
end
%
opts_used=opts;
stats=psg_session_stats(sessions,setfield(opts,'if_log',0)); %will need statistics to do randomization
%create permutations, even if they won't be used
perms_used.allsess=cell(opts.cond_nstims,1);
perms_used.eachsess=cell(opts.cond_nstims,opts.cond_nsess);
for istim=1:opts.cond_nstims
    for isess=1:opts.cond_nsess
        perms_used.eachsess{istim,isess}=randperm(stats.counts_alluses(istim,isess))-1;
    end
    perms_used.allsess{istim}=randperm(sum(stats.counts_alluses(istim,:),2))-1;
end
ntrials=size(sessions,1);
ncols=size(sessions,2);
nsess=size(sessions,3);
session_cell=cell(ntrials,ncols);
session_cells=cell(nsess,1);
examps_used=zeros(size(sessions));
for isess=1:nsess
    session_cell=cell(ntrials,ncols);
    for icol=1:ncols
        for itrial=1:ntrials
            stimtype=sessions(itrial,icol,isess);
            entry=typenames{stimtype};
            switch opts.example_infix_mode
                case 1 %unique stimuli across all sessions
                    loc=(isess-1)*ncols*ntrials+(icol-1)*ntrials+itrial; %location in 3d array
                    example_seq=sum(sessions(1:loc-1)==stimtype); %which occurrence is this across all sessions?
                    example_perm=perms_used.allsess{stimtype}(1+example_seq);
                case 2 %unique stimuli within sessions
                    loc=(icol-1)*ntrials+itrial; %location in 2d array
                    example_seq=sum(sessions((isess-1)*ncols*ntrials+[1:loc-1])==stimtype) ;%which occurrence is this within this session?
                    example_perm=perms_used.eachsess{stimtype,isess}(1+example_seq);
                otherwise
                    example_perm=0;
            end
            examps_used(itrial,icol,isess)=example_perm;
            if (opts.example_infix_mode)~=4
                entry=cat(2,entry,opts.example_infix_string,zpad(example_perm,opts.example_infix_zpad));
            end
            session_cell{itrial,icol}=entry;
        end %icol
    end %itrial
    session_cells{isess}=session_cell;
end %isess
return

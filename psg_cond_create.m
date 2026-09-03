function [session_cells,perms_used,examps_used,opts_used]=psg_cond_create(sessions,typenames,opts)
% [session_cells,perms_used,examps_used,opts_used]=psg_cond_create(sessions,typenames,opts)
% creates text arrays for condition files from numerical arrays of stimulus types
%
% sessions: integer array [ntrials 1+ncompares nsess] of stimulus types
% typenames: cell array of text strings for each stimulus type
% opts: options field, see psg_defopts
%   opts.example_infix_mode: controls re-use of stimuli
%   opts.example_infix_string: separator between stimulus name and example number, typically underscore (_)
%   opts.example_infix_zpad: number of digits to zero-pad 
%   opts.prefix: a prefix for all file names
%   opts.example_numoffset: optional global offset to avoid repeats across differnt setups
%   opts.example_maxnum: maximum number of examples available for each stimtype, scalar or vector of length equal to opts.cond_nstims
%   opts.example_list: opts_example_list{istim} is the list of examples available;
%       if not supplied, assumed to be 0:example_maxnum-1 or [0:example_maxnum(istim)-1]
%
% session_cells: session_cells{isess}(itrial,istim} is a cell array of cell arrays that can be used for the csv cond file, specifying the stimulus file.
%     Stimulus file is [opts.prefix,typenames{stimtype},opts.example_infix_string,then stimulus instance, zero-padded via opts.example_infix_zpad]
% perms_used: permutations used, zero-indexed (computed even if not needed)
%   perms_used.allsess{istim}: across all sessions, needed for unique stimuli across all trials in all sessions
%   perms_used.eachsess{istim,isess}: across each session, needed for unique stimuli within a session
% examps_used: examples used, zero-indexed, as a pointer into example number. size is [ntrials 1+ncompares nsess]
%   
% opts_used: options used
%
% 07Jul23: opts.prefix added
% 11Sep23: example_numoffset added, affects session_cells but not perms_used or examps_used
% 14Aug26: changes for compatibility with mater_psg_sess_setup: 
%      different numbers of examples available for each stimulus via opts.example_maxnum
%      allow for exclusion of immedeate repeats of stimuli if number of available examples is >=3, via opts.exclude_immed_repeats
%      allow for an explicit listing of example numbers via opts.example_list
%         
%   See also:  PSG_SETUP_DEMO, PSG_COND_WRITE, PSG_SESSCONFIG_MAKE, ZPAD, PSG_SESSION_STATS.
%
if (nargin<3)
    opts=struct;
end
%
opts=filldefault(opts,'prefix','');
opts=filldefault(opts,'exclude_immed_repeats',0);
opts=filldefault(opts,'example_maxnum',Inf);
opts_used=opts;
%
stats=psg_session_stats(sessions,setfield(opts,'if_log',0)); %will need statistics to do randomization
%create permutations, even if they won't be used
perms_used.allsess=cell(opts.cond_nstims,1);
perms_used.eachsess=cell(opts.cond_nstims,opts.cond_nsess);
for istim=1:opts.cond_nstims
    %stats.counts_alluses{istim,isess} is the number of times that stimulus istim is used in session isess
    for isess=1:opts.cond_nsess
        perms_used.eachsess{istim,isess}=randperm(stats.counts_alluses(istim,isess))-1;
    end
    perms_used.allsess{istim}=randperm(sum(stats.counts_alluses(istim,:),2))-1;
end
ntrials=size(sessions,1);
ncols=size(sessions,2);
nsess=size(sessions,3);
session_cells=cell(nsess,1);
examps_used=zeros(size(sessions));
exclude_offset=0;
last_used=NaN(1,opts.cond_nstims);
for isess=1:nsess
    session_cell=cell(ntrials,ncols);
    for itrial=1:ntrials %trial is outer loop beginning 14Aug26, to enable excluding immediate repeats of same example
        for icol=1:ncols
            stimtype=sessions(itrial,icol,isess);
            entry=cat(2,opts.prefix,typenames{stimtype});
            switch opts.example_infix_mode
                case 1 %unique stimuli across all sessions
                    loc=(isess-1)*ncols*ntrials+(icol-1)*ntrials+itrial; %location in 3d array
                    example_seq=sum(sessions(1:loc-1)==stimtype); %which occurrence is this across all sessions?
                    example_perm=perms_used.allsess{stimtype}(1+example_seq)+exclude_offset;
                case 2 %unique stimuli within sessions
                    loc=(icol-1)*ntrials+itrial; %location in 2d array
                    example_seq=sum(sessions((isess-1)*ncols*ntrials+[1:loc-1])==stimtype) ;%which occurrence is this within this session?
                    example_perm=perms_used.eachsess{stimtype,isess}(1+example_seq)+exclude_offset;
                otherwise
                    example_perm=0;
            end
            %take care of possibility that different stimuli have different numbers of available examples
            if length(opts.example_maxnum)==1
                example_maxnum=opts.example_maxnum;
            else
                example_maxnum=opts.example_maxnum(stimtype);
            end
            if example_maxnum<Inf & example_maxnum>0
                example_perm=mod(example_perm,example_maxnum);
            end
            %if requested, exclude immediate repeats
            if example_maxnum>=3 & opts.exclude_immed_repeats==1 & ismember(opts.example_infix_mode,[1 2])
                if example_perm==last_used(stimtype)
                    new_offset=1+floor((example_maxnum-1)*rand(1)); %1, ..., example_maxnum-1, to ensure a different example of this stimulus
                    example_perm=mod(example_perm+new_offset,example_maxnum); %choose a different example
                    exclude_offset=exclude_offset+new_offset;
                end
            end
            examps_used(itrial,icol,isess)=example_perm;
            last_used(stimtype)=example_perm; %to avoid immediate repeats, if requested
            example_index=example_perm+opts.example_numoffset; %
            if isfield(opts,'example_list')
                example_number=opts.example_list{stimtype}(example_index+1);
            else
                example_number=example_index;
            end
            if (opts.example_infix_mode)~=4
                entry=cat(2,entry,opts.example_infix_string,zpad(example_number,opts.example_infix_zpad));
            end
            session_cell{itrial,icol}=entry;
        end %icol
    end %itrial
    session_cells{isess}=session_cell;
end %isess
return

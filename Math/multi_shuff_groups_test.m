% multi_shuff_groups_test: tests multi_shuff_groups
%
%  See also:  MULTI_SHUFF_GROUPS, MULTI_SHUFF_ENUM.
if ~exist('opts_multi')
    opts_multi=struct;
end
opts_multi=filldefault(opts_multi,'if_log',0);
opts_multi=filldefault(opts_multi,'if_sortrows',0);
opts_multi.if_ask=0;
%
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
if (if_frozen~=0) 
    rng('default');
    if (if_frozen<0)
        rand(1,abs(if_frozen));
    end
else
    rng('shuffle');
end
%
if ~exist('setups')
    setups=cell(0);
    %
    setups{1}.gps= [1 1 1 2 2 2];
    setups{1}.tags=[1 1 2 2 2 1];
    %
    setups{2}.gps= [1 1 1 2 2 2];
    setups{2}.tags=[1 1 2 1 1 2];
    % non-consecutive, and missing elements
    setups{3}.gps= [7 7 7 7 3 3 3 3];
    setups{3}.tags=[6 4 1 2 4 2 6 1];
    % non-consecutive, and missing elements
    setups{4}.gps= [6 4 1 2 4 2 6 1];
    setups{4}.tags=[7 7 7 7 3 3 3 3];
    % lots of reduction
    setups{5}.gps= [1 1 1 2 2 2 3 3 3];
    setups{5}.tags=[8 3 4 3 4 8 4 3 8];
    %two group lengths
    setups{6}.gps= [1 1 2 2 3 3 3 4 4 4];
    setups{6}.tags=[];
    %two group lengths, bigger
    setups{6}.gps= [1 1 2 2 3 3 3 7 4 4 4 7];
    setups{6}.tags=[];
    %five small groups
    setups{7}.gps= [1 1 2 2 3 3 4 4 5 5];
    setups{7}.tags=[1 1 1 1 2 2 1 1 2 2];
end
nsetups=length(setups);
nvariants=6;
results=cell(nsetups,nvariants);
for isetup=1:nsetups
    disp(sprintf(' doing setup %2.0f',isetup));
    disp(' groups:')
    disp(setups{isetup}.gps);
    disp(' tags:')
    disp(setups{isetup}.tags);
    for ivar=1:nvariants
        opts_multi_use=opts_multi;
        switch ivar
            case 1 %exhaustive, no tags, no reduction
                label='exhaustive, no tags, no reduction';
                opts_multi_use.if_exhaust=1;
                opts_multi_use.tags=[];
                opts_multi_use.if_reduce=0;
            case 2 %exhaustive, no tags, reduction
                label='exhaustive, no tags, reduction';
                opts_multi_use.if_exhaust=1;
                opts_multi_use.tags=[];
                opts_multi_use.if_reduce=1;
            case 3 %exhaustive, tags, no reduction
                label='exhaustive, tags, no reduction';
                opts_multi_use.if_exhaust=1;
                opts_multi_use.tags=setups{isetup}.tags;
                opts_multi_use.if_reduce=0;
            case 4 %exhaustive, tags, reduction
                label='exhaustive, tags, reduction';
                opts_multi_use.if_exhaust=1;
                opts_multi_use.tags=setups{isetup}.tags;
                opts_multi_use.if_reduce=1;
            case 5 %random, no tags
                label='random, no tags';
                opts_multi_use.if_exhaust=0;
                opts_multi_use.tags=[];
            case 6 %random, tags
                label='random, tags';
                opts_multi_use.if_exhaust=0;
                opts_multi_use.tags=setups{isetup}.tags;
        end %ivar
        %
        tic;
        [shuffs,gp_info,opts_used]=multi_shuff_groups(setups{isetup}.gps,opts_multi_use);
        elapsed=toc;
        %correct number made?
        shuffs_made=size(shuffs,1);
        ok_string_count='OK';
        if shuffs_made~=opts_used.nshuffs
            ok_string_count='BAD';
        end
        %are the permutations unique (if expected to be unique)?
        ok_string_unique='OK';
        unique_shuffs_made=size(unique(shuffs,'rows'),1);
        if shuffs_made~=unique_shuffs_made
            if opts_multi_use.if_exhaust==1
                ok_string_unique='BAD';
            else
                ok_string_unique='DUP';
            end
        end
        disp(sprintf('% 40s: expected: %10.0f, found: %10.0f (%3s), unique: %10.0f (%3s); time: %12.5f sec',...
            label,opts_used.nshuffs,shuffs_made,ok_string_count,unique_shuffs_made,ok_string_unique,elapsed));
        results{isetup,ivar}.label=label;
        results{isetup,ivar}.setup=setups{isetup};
        results{isetup,ivar}.shuffs=shuffs;
        results{isetup,ivar}.gp_info=gp_info;
        results{isetup,ivar}.opts_used=opts_used;
    end
end

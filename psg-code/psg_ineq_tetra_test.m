% psg_ineq_tetra_test: test whether the inequality structure has the requisite symmetry
%
%   See also:  PSG_INEQ_LOGIC_TETRA, PSG_INEQ_LOGIC, PSG_INEQ_LOOKUP.
%
if ~exist('opts') opts=struct;end
opts=filldefault(opts,'if_log',0);
triad_table=psg_ineq_triads('tetra');
nc=12;
nv=4; %four stimuli
%
all_types=fieldnames(psg_ineq_logic());
tetra_types=cell(0);
for it=1:length(all_types)
    if contains(all_types{it},'tetra')
        tetra_types{end+1}=all_types{it};
    end
end
for it=1:length(tetra_types)
    disp(sprintf('analyzing %s',tetra_types{it}));
    partitions=psg_ineq_logic(nc,tetra_types{it});
    perm_list=sortrows(perms(1:nv)); %a permutation
    for ip=1:size(perm_list,1)
        perm_use=perm_list(ip,:);
        triad_table_perm=perm_use(triad_table);
        %
        %create partitions with the permuted indices, but then check that it is unchanged
        %
        [triad_ptrs,if_flips]=psg_ineq_lookup(triad_table,triad_table_perm); %which dimensions will this be?
        if any(triad_ptrs==0)
            disp(sprintf(' permutation %1.0f %1.0f %1.0f %1.0f: some triads not found',perm_use))
        else
            %
            %note that here, we first inverse-permute and then do the flips
            %(to create the tetra tables, we first flip and then permute)
            %
            partitions_toflip=permute(partitions,triad_ptrs);
            partitions_new=partitions_toflip;
            %flip as needed
            for id=1:nc
                if if_flips(triad_ptrs(id))
                    partitions_new=flip(partitions_new,id); %which dimensions need to flip
                end
            end
            ndiff=sum(double(partitions_new(:)~=partitions(:)));
            disp(sprintf(' permutation %1.0f %1.0f %1.0f %1.0f: %7.0f mismatches',perm_use,ndiff))
        end
    end %ip
end

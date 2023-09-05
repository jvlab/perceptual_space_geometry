% psg_ineq_tetra_test: test whether the inequality structure has the requisite symmetry
%
%   See also:  PSG_INEQ_LOGIC_TETRA, PSG_INEQ_LOGIC, PSG_INEQ_LOOKUP.
%
sources=struct;
%the inequality structures computed via psg_ineq_logic_tetra, and where they come from 
sources.exclude_trans_tetra='exclude_trans_tent';
sources.exclude_addtree_tetra='exclude_addtree';
sources.exclude_addtree_trans_tetra='exclude_addtree_trans';
sources.exclude_trans_tetra_cyc4='exclude_trans_cyc4';
%
if ~exist('opts') opts=struct;end
opts=filldefault(opts,'if_log',0);
triad_table=psg_ineq_triads('tetra');
nt=4; %four points
%
%
all_types=fieldnames(psg_ineq_logic());
tetra_types=cell(0);
for it=1:length(all_types)
    if contains(all_types{it},'tetra')
        tetra_types{end+1}=all_types{it};
    end
end
tetra_types=fieldnames(sources);
nc=12; %number of ineqs (dimensions)
nv=4; %number of tents
for it=1:length(tetra_types)
    disp(sprintf('analyzing %s',tetra_types{it}));
    partitions=psg_ineq_logic(nc,tetra_types{it});
    perm_list=sortrows(perms(1:nv)); %a permutation
    mismatch_range=[0 0];
    for ip=1:size(perm_list,1)
        perm_use=perm_list(ip,:);
        triad_table_perm=perm_use(triad_table);
        %
        %create partitions with the permuted indices, but then check that relabeling leaves partitions unchanged
        %
        [triad_ptrs,if_flips]=psg_ineq_lookup(triad_table,triad_table_perm); %which dimensions will this be?
        if any(triad_ptrs==0)
            disp(sprintf(' permutation %1.0f %1.0f %1.0f %1.0f: some triads not found',perm_use))
        else
            %
            partitions_new=partitions;
            %flip as needed
            for id=1:nc
                if if_flips(id)
                    partitions_new=flip(partitions_new,id); %which dimensions need to flip
                end
            end
            partitions_new=ipermute(partitions_new,triad_ptrs); %note inverse permutation
            ndiff=sum(double(partitions_new(:)~=partitions(:)));
            if (ndiff>0)
                disp(sprintf(' permutation %1.0f %1.0f %1.0f %1.0f: %7.0f mismatches',perm_use,ndiff))
            end
            mismatch_range=[min(mismatch_range(1),ndiff),max(mismatch_range(2),ndiff)];
        end
    end %ip
    disp(sprintf('mismatches across permutations: %7.0f to %7.0f',mismatch_range))
    %check source from psg_ineq_logic_tetra
    if isfield(sources,tetra_types{it})
       [partitions_check,tetra_table,psub,ou_tetra]=psg_ineq_logic_tetra(sources.(tetra_types{it}));
        nbad=sum(double(partitions(:)~=partitions_check(:)));
        disp(sprintf('verifying %30s against re-creation from %25s: %8.0f mismatches',tetra_types{it},sources.(tetra_types{it}),nbad));
        nv=length(psub);
        for jv=2:nv
            for iv=1:jv-1
                ndiffs=sum(double(psub{iv}(:)~=psub{jv}(:)));
                for k=1:nt %inverse permutation
                    newtent_iv(k)=find(ou_tetra.tents(iv,:)==k);
                    newtent_jv(k)=find(ou_tetra.tents(jv,:)==k);
                end
                disp(sprintf(' differences for component %3.0f (%1.0f %1.0f %1.0f -> %1.0f) and %3.0f (%1.0f %1.0f %1.0f-> %1.0f): %8.0f',...
                    iv,newtent_iv,jv,newtent_jv,ndiffs));
            end
        end
    else
        disp(sprintf('source of %s not found',tetra_types{it}));
    end
    disp(' ');
end

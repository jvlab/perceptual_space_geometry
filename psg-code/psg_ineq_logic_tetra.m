function [partitions,triad_table,psub]=psg_ineq_logic_tetra(ineq_type,opts)
% [partitions,triad_table,psub]=psg_ineq_logic_tetra(ineq_type,opts) create a 3^12-element array 
% for the excluded inequalities for a tetrahedron, from its component tripods
% addtree mask of a tetrahedron
%
% ineq_type: inequality type for the tripod
%   exclude_trans_tent: conditions that symmetry and transitivity in a tent
%   exclude_addtree_trans: conditions that exclude addtree inequality or transitivity
% opts: if_log=1 to log intermediate calculations
%
% partitions: an array of size repmat(3,1,12) of 0's and 1's, specifiying
%    the forbidden inequality combinations in a tetrahedron (1: forbidden) 
% triad_table: array of size [12 3] indicating the triads for the 12 dimensions of partitions
% psub: the masks corresponding to each of the four tripods that overlap in the tetrahedron
%
% See also:  PSG_INEQ_LOGIC, PSG_INEQ_TRIADS, PSG_INEQ_LOOKUP.
%
if (nargin<2)
    opts=[];
end
opts=filldefault(opts,'if_log',0);
triad_table=psg_ineq_triads('tetra');
tent_table=psg_ineq_triads('tent');
nc=12;
nc_tripod=6;
nv=4; %four vertices
psub=cell(1,nv);
p_tent=psg_ineq_logic(nc_tripod,ineq_type);  
p_tent_aug=repmat(p_tent,[ones(1,nc_tripod),3*ones(1,nc-nc_tripod)]); %replicate on unused dimensions
for iv=1:nv
    psub{iv}=zeros(repmat(3,1,nc));
    %determine the exclusions for a tent with vertex at iv 
    tent_table_mod=mod(tent_table+iv-1,nv)+1;
    [triad_ptrs,if_flips]=psg_ineq_lookup(triad_table,tent_table_mod); %which dimensions will this be?
    %triad_ptrs points to where each dimension of p_tent belongs in the tetrahedral partition
    triad_ptrs_aug=[triad_ptrs;setdiff([1:nc],triad_ptrs)']; %unused dimensions at end
    if_flips_aug=[if_flips;zeros(nc-nc_tripod,1)];
    if (opts.if_log)
        disp(sprintf('calculations for vertex %1.0f',iv));
        disp('triads (dimensions)');
        disp(triad_ptrs')
        disp('if_flips');
        disp(if_flips')
        disp('triad_ptrs_aug')
        disp(triad_ptrs_aug');
    end
    p_tent_aug_flipped=p_tent_aug;
    for id=1:nc_tripod
        if if_flips(id)
            p_tent_aug_flipped=flip(p_tent_aug_flipped,id); %which dimensions need to flip
        end
    end
    psub{iv}=permute(p_tent_aug_flipped,triad_ptrs_aug); %
    if (opts.if_log)
        for id=1:nc
            margsum=sum(psub{iv},setdiff([1:nc],id));
            disp(sprintf(' vertex %3.0f dim %3.0f: marginal sums: %8.0f %8.0f %8.0f, from tripod dim %3.0f with flip %1.0f',...
                iv,id,margsum(:)',triad_ptrs_aug(id),if_flips_aug(id)));
        end
    end
end
partitions=min(psub{1}+psub{2}+psub{3}+psub{4},1);
return

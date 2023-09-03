function [partitions,tetra_table,psub,opts_used]=psg_ineq_logic_tetra(ineq_type,opts)
% [partitions,tetra_table,psub,opts_used]=psg_ineq_logic_tetra(ineq_type,opts) create a 3^12-element array 
% for the excluded inequalities for a tetrahedron, from its component tripods
% addtree mask of a tetrahedron
%
% ineq_type: inequality type for the tripod
%   exclude_trans_tent: conditions that symmetry and transitivity in a tent
%   exclude_addtree_trans: conditions that exclude addtree inequality or transitivity
% opts: 
%   if_log=1 to log intermediate calculations
%   tents: array of size [:,4] specifying, in each row, the vertices to be used as tents.
%   defaults to [1 2 3 4;2 3 4 1;3 4 2 1;4 1 2 3];
% partitions: an array of size repmat(3,1,12) of 0's and 1's, specifiying
%    the forbidden inequality combinations in a tetrahedron (1: forbidden) 
% tetra_table: array of size [12 3] indicating the triads for the 12 dimensions of partitions
% psub: the masks corresponding to each of the four tripods that overlap in the tetrahedron
%
% See also:  PSG_INEQ_LOGIC, PSG_INEQ_TRIADS, PSG_INEQ_LOOKUP, PSG_INEQ_TETRA_TEST.
%
if (nargin<2)
    opts=[];
end
opts=filldefault(opts,'if_log',0);
opts=filldefault(opts,'tents',[1 2 3 4;2 3 4 1;3 4 2 1;4 1 2 3]);
opts_used=opts;
%
tetra_table=psg_ineq_triads('tetra');
tent_table=psg_ineq_triads('tent');
nt=4;
nc=12;
nc_tripod=6;
nv=size(opts.tents,1); %typically four vertices
psub=cell(1,nv);
p_tent=psg_ineq_logic(nc_tripod,ineq_type);  
p_tent_aug=repmat(p_tent,[ones(1,nc_tripod),3*ones(1,nc-nc_tripod)]); %replicate on unused dimensions
%
triad_ptrs=zeros(size(tent_table,1),nv);
triad_ptrs_aug=zeros(size(tetra_table,1),nv);
if_flips=zeros(size(tent_table,1),nv);
if_flips_aug=zeros(size(tetra_table,1),nv);
for iv=1:nv     %determine the exclusions for the tent with vertices [1 2 3 4] relabeled as per opts.tents{iv,:)
    psub{iv}=zeros(repmat(3,1,nc));
    tents_use=opts.tents(iv,:);
    %
    % tents_use(tetra_table) has the 12 comparisons of a tetrahedron, but
    % with tokens relabeled per tents_use.  We then find the rows that are
    % in a standard tent (i.e., 4 at the vertex, [123] around the base)
    % and then, via permute, put them in thr proper rows of psub{iv}.
    % So psub{iv} has the constraints corresponding to a [123]->4 tent, relabled.
    %
    [triad_ptrs_aug(:,iv),if_flips_aug(:,iv)]=psg_ineq_lookup(tent_table,tents_use(tetra_table)); %
    triad_ptrs_aug(find(triad_ptrs_aug(:,iv)==0),iv)=[nc_tripod+1:nc]';
    %
    if (opts.if_log)
        disp(sprintf('calculations for tent %1.0f',iv));
        disp(sprintf('tent: %1.0f %1.0f %1.0f %1.0f',opts.tents(iv,:)));
        disp(tents_use(tetra_table));
        disp('if_flips_aug');
        disp(if_flips_aug(:,iv)')
        disp('triad_ptrs_aug')
        disp(triad_ptrs_aug(:,iv)');
    end
    p_tent_aug_flipped=p_tent_aug;
    psub{iv}=permute(p_tent_aug_flipped,triad_ptrs_aug(:,iv)); 
    for id=1:nc
        if if_flips_aug(id,iv)
            psub{iv}=flip(psub{iv},id); %which dimensions need to flip
        end
    end
%
    if (opts.if_log)
        for id=1:nc
            margsum=sum(psub{iv},setdiff([1:nc],id));
            minvals=min(psub{iv},[],id);
            maxvals=max(psub{iv},[],id);
            if_const=double(all(minvals(:)==maxvals(:)));
            disp(sprintf(' tent %3.0f dim %3.0f: marginal sums: %8.0f %8.0f %8.0f, from tripod dim %3.0f with flip %1.0f, if_const=%2.0f',...
                iv,id,margsum(:)',triad_ptrs_aug(id),if_flips_aug(id),if_const));
        end
    end
end
if (opts.if_log)
    for jv=2:nv
        for iv=1:jv-1
            ndiffs=sum(double(psub{iv}(:)~=psub{jv}(:)));
            for k=1:nt %inverse permutation
                newtent_iv(k)=find(opts.tents(iv,:)==k);
                newtent_jv(k)=find(opts.tents(jv,:)==k);
            end
            disp(sprintf(' differences for tent %3.0f (%1.0f %1.0f %1.0f -> %1.0f) and %3.0f (%1.0f %1.0f %1.0f-> %1.0f): %8.0f',...
                iv,newtent_iv,jv,newtent_jv,ndiffs));
        end
    end
end
partitions=zeros(repmat(3,1,nc));
for iv=1:nv
    partitions=partitions+psub{iv};
end
partitions=min(partitions,1);
return

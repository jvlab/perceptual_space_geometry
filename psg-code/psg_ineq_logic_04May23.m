function [partitions,opts_used]=psg_ineq_logic(nc,ineq_type,opts)
% [partitions,opts_used]=psg_ineq_logic(nc,ineq_type,opts) sets up the
% partitions of a "probability n-cube" that satisfy a set of inequalities
% and thereby exclude a model for a perceptual geometry
%
% in all cases, transitivity of comparing numbers is assumed, i.e,. if x>y and y>z then x>z
%
% if called with no arguments, partitions returns a structur of the allowed
%   values of ineq_type, each of which has a field nc indicating the allowed
%   range of nc
% nc: number of rank-choice probability estimates, typically 3 or 6
%   If nc=3, the three rank-choice probabilities are in the order
%   P(d(a,b)<d(a,c)), P(d(b,c)<d(b,a)), P(d(c,a)<d(c,b)), so 
%   values that are < 1/2, > 1/2, and = 1/2 means that
%   d(a,b)>d(a,c), d(b,c) < d(b,a), d(c,a)=d(c,b) which is INconsistent with the UMI,
%   since this is an isosceles triangle with the unequal side (a,b) the longest.
%      Consistency:
%         partitions(2,3,1)=partitions(3,1,2)=partitions(1,2,3)=0
%      For (2,3,1): d(a,b)=d(a,c); d(b,c)< d(b,a); d(c,a)>d(c,b), i.e.,
%      isosceles with d(a,b) longest.
%  If nc=6, the six rank-choice probabilities are in the order
%  P(d(z,b)<d(z,c)), P(d(z,c)<d(z,a)), P(d(z,a)< d(z,b), P(d(a,b)<d(a,c)), P(d(b,c)<d(b,a)), P(d(c,a)<d(c,b))
% ineq_type: type of inequality
%   'all': covers all partitions (just for checking)
%   'exclude_sym'           (nc=3): conditions that exclude symmetry (transitivity assumed)
%   'exclude_trans'         (nc=3): conditions that exclude transitivity (symmetry of distance assumed) -- same as exclude_sym
%   'exclude_umi'           (nc=3): conditions that exclude the ultrametric inequality but may
%        allow for non-transitivity (symmetry of distance assumed)
%   'exclude_umi_trans'     (nc=3): conditions that exclude the ultrametric inequality or transitivity
%      (computed recursively from exclude_umi, exclude_trans)
%   'exclude_trans_tent     (nc=6): conditions that exclude transitivity in a tent
%      (computed recursively from exclude_trans)
%   'exclude_addtree'       (nc=6): conditions that exclude addtree inequality (symmetry of distance assumed)
%   'exclude_addtree_trans' (nc=6): conditions that exclude addtree inequality or transitivity
%      (computed recursively from exclude_addtree, exclude_trans_tent)
% opts: options
%   opts.if_log: 1 to log, -1 to not log and not calculate edge_counts or partitions_nz
%
% partitions: array of nc dimensions (one for each rank-choice),
%  3 values on each dimension, corresponding to rank choice< 1/2, = 1/2, > 1/2
%  (size is 3^nc)
%  An entry of 0 means that this combination is NOT excluded
%  An entry of 1 means that this combination is part of the exclusion
% opts_used: options used
% opts_used.edge_counts: array of size nc+1, edge_counts(k+1) is the number of nonzero entries in
%  partitions that have k coordinates = 2 (i.e., on the rank choice=1/2 edge)
% opts_used.partitions_nz: an array of nc dimensions, 2 values on each dimension, 
%  excluding the "=1/2" values of partitions
%
% note that trans_tent, addtree, and addtree_trans do not have tetrahedral symmetry -- transformations for swapping a and z:
% orig P(d(z,b)<d(z,c)), P(d(z,c)<d(z,a)), P(d(z,a)< d(z,b), P(d(a,b)<d(a,c)), P(d(b,c)<d(b,a)), P(d(c,a)<d(c,b))
% swap P(d(a,b)<d(a,c)), P(d(a,c)<d(a,z)), P(d(a,z)< d(a,b), P(d(z,b)<d(z,c)), P(d(b,c)<d(b,z)), P(d(c,z)<d(c,b))
% d1 and d4 interchange; d2,d3,d5,d6 are unique
%
% 11Mar23: add special call with no arguments
% 12Mar23: check for legitimate arguments based on avail_types, add descriptive fields to avail_types.(ineq_type)
%
% See also:  PSG_UMI_TRIPLIKE, PSG_UMI_TRIPLIKE_DEMO, PSG_CHECK_PROBS, PSG_TENTLIKE_DEMO, PSG_INEQ_APPLY, PSG_PROBS_CHECK,
%   PSG_INEQ_EDGECOUNT, PSG_INEQ_LOGIC_DEMO.
%
flip_invar={'all','none','exclude_sym','exclude_trans','exclude_trans_tent'};  %logic types that are unchanged if < replaced by >
cycle_invar={'all','none','exclude_sym','exclude_trans','exclude_umi','exclude_umi_trans'}; %logic types that are unchanged if dimensions are cycled
avail_types=struct();
avail_types.all.nc=[];
avail_types.none.nc=[];
avail_types.exclude_sym.nc=3;
avail_types.exclude_trans=avail_types.exclude_sym;
avail_types.exclude_umi.nc=3;
avail_types.exclude_umi_trans.nc=3;
avail_types.exclude_trans_tent.nc=6;
avail_types.exclude_addtree.nc=6;
avail_types.exclude_addtree_trans.nc=6;
ineq_types=fieldnames(avail_types);
for iq=1:length(ineq_types)
    fn=ineq_types{iq};
    avail_types.(fn).flip_invar=double(~isempty(strmatch(fn,flip_invar,'exact')));  
    avail_types.(fn).cycle_invar=double(~isempty(strmatch(fn,cycle_invar,'exact')));
end
if (nargin==0)
    partitions=avail_types;
    return
end
if (nargin<=2)
    opts=struct;
end
%
if ~isfield(avail_types,ineq_type)
    warning(sprintf('inequality type %s unrecognized',ineq_type));
    return
else
    nc_allowed=avail_types.(ineq_type).nc;
    if ~isempty(nc_allowed)
        if ~ismember(nc,nc_allowed)
             warning(sprintf('wrong number of rank-choice probabilities (%1.0f) for %s',nc,ineq_type));
             return
        end
    end
end
%    
threes=[repmat(3,1,nc) 1];
opts=filldefault(opts,'if_log',0);
partitions=zeros(threes);
if_warn=0;
lt=1; %token for rank choice probabilty less than 1/2 (first distance judged as >)
eq=2; %token for rank choice probabilty equal to 1/2 (first distance judged as =)
gt=3; %token forrank  choice probability greater than 1/2 (first distance judged as <)
le=[lt:eq]; %rank choice prob <=1/2, first distance >= second distance
ge=[eq:gt]; %rank choice prob >=1/2, first distance <= second distance
switch ineq_type
    case 'all'
        partitions(:)=1;
    case 'none'
    case {'exclude_sym','exclude_trans'} 
        %symmetry is excluded (assuming transitivity) if there is a consistent cycle of inequalities
        %and at least one inequality is strict
        %
        %transitivity is excluded (assuming symmetry) under the same
        %conditions, which can also be phrased
        % two equalities and one strict inequality
        % one equality and two strict inequalities of same sign
        % no equalities; three strict inequalities of the same sign
        for id=1:3
            is=double([1:nc]~=id); %one inequality must be strict
            partitions(lt:lt+is(1),lt:lt+is(2),lt:lt+is(3))=1;
            partitions(gt-is(1):gt,gt-is(2):gt,gt-is(3):gt)=1;
        end
    case 'exclude_umi'
        % ultrametric inequality excluded if there is one equality and two strict inequalities
        % implying an isosceles triangle with the odd leg longest
        % partitions(1,3,2)=1;
        % partitions(3,2,1)=1;
        % partitions(2,1,3)=1;                       
        %ultrametric inequality is excluded
        % if there are two equalities and one strict inequality
        % partitions([1 3],2,2)=1;
        % partitions(2,[1 3],2)=1;
        % partitions(2,2,[1 3])=1; 
        for id=1:3
             partitions(id,1+mod(id+1,3),1+mod(id,3))=1;
             im=double([1:nc]==id);
             partitions([eq-im(1) eq+im(1)],[eq-im(2) eq+im(2)],[eq-im(3) eq+im(3)])=1;
        end
        %ultrametric inequality is excluded if there are no equalities, with inconsistent sign
        %    partitions([1 3],[1 3],[1 3])=1;
        partitions([lt gt],[lt gt],[lt gt])=1;
        %
        %partitions(1,1,1)=0; %consistent sign (but this is inconsistent with transitivity)
        %partitions(3,3,3)=0; %consistent sign (but this is inconsistent with transitivity)           
        partitions(1)=0;
        partitions(end)=0;
    case 'exclude_umi_trans'
        p_umi=psg_ineq_logic(nc,'exclude_umi');
        p_trans=psg_ineq_logic(nc,'exclude_trans');
        partitions=p_umi | p_trans;
    case 'exclude_trans_tent'
        p_trans=psg_ineq_logic(3,'exclude_trans');
        %exclude transitivity if either the tripod part (first 3 coords) or the triangle part (last 3 coords)
        %are inconsistent with transitivity
        partitions=repmat(p_trans,[1 1 1 3 3 3]) | repmat(reshape(p_trans,[1 1 1 3 3 3]),[3 3 3 1 1 1]); 
     case 'exclude_addtree'
        %for each rotation order of dimension, cannot have one pairsum strictly greater than the other two
        %  P(d(z,b)<d(z,c)), P(d(z,c)<d(z,a)), P(d(z,a)< d(z,b), P(d(a,b)<d(a,c)), P(d(b,c)<d(b,a)), P(d(c,a)<d(c,b))
        %
        z6=zeros(repmat(3,1,nc));
        %d(z,c)>=d(z,a) and d(a,b)>=d(b,c), at least one inequality strict
        zc_ab_gt_za_bc=z6;
        zc_ab_gt_za_bc(:,le,:,:,ge,:)=1;
        zc_ab_gt_za_bc(:,eq,:,:,eq,:)=0;
        %d(z,c)>=d(z,b) and d(a,b)>=d(c,a), at least one inequality strict
        zc_ab_gt_zb_ca=z6;
        zc_ab_gt_zb_ca(ge,:,:,le,:,:)=1;
        zc_ab_gt_zb_ca(eq,:,:,eq,:,:)=0;
        %
        zc_ab_largest=and(zc_ab_gt_za_bc,zc_ab_gt_zb_ca);
        %
        %now apply this after cycling by [a,b,c]
        partitions=zc_ab_largest | permute(zc_ab_largest,[2 3 1 5 6 4]) | permute(zc_ab_largest,[3 1 2 6 4 5]);           
    case 'exclude_addtree_trans'
        p_addtree=psg_ineq_logic(nc,'exclude_addtree');
        p_trans_tent=psg_ineq_logic(nc,'exclude_trans_tent');
        partitions=p_addtree | p_trans_tent;
end
if opts.if_log>=0
    %
    %create a template of the number of edges in each location to use as a mask
    %
    z=[0 1 0];
    for k=2:nc
        z=[z z+1 z];
    end
    opts.edge_counts=psg_ineq_edgecount(partitions,opts.if_log);
    pz=partitions(find(z==0));
    opts.partitions_nz=reshape(pz,[repmat(2,1,nc) 1]);
end
opts_used=opts;
return

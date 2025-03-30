function [mixent,mixmats,shuff_gp_origs]=multi_shuff_mixent(shuffs,gps)
% [mixent,mixmats,shuff_gp_origs]=[shuffs,gp_info,opts_used]=multi_shuff_mixent(shuffs,gps)
% determines group mixing and mixing entropy of shuffles
% 
% shuffs: the shuffles; each row is a perm of [1:n], and shuffles(is,:) goes into gps(:)
%   nshuffs=size(shuffs,1)
% gps: a row vector of length n, that indicates the assignment of each
%   element to [1:ngps].  Not all of [1:ngps] must be present.
%   size(gps,2) must match size(shuffs,2)
%
% mixent:column vector, length nshuffs. mixents(ishuff) is entropy (log 2) of mixing
%  matrix mixmats(:,:,ishuff).  zero if no mixing.
% mixmats(igp,jgp,ishuff) is the number of elements in initial group igp that
%   come from group jgp in shuffle(ishuff,:)
% shuff_gp_origs: cell array{1,ngps}, shufff_gp_orig{igp}(ishuff,:) is original group membership after shuffling
% 
%  See also: MULTI_SHUFF_GROUPS, HLID_GEOM_TRANSFORM_STATS, LOGNZ.
%
nshuffs=size(shuffs,1);
ngps=max(gps);
%
mixent=zeros(nshuffs,1);
mixmats=zeros(ngps,ngps,nshuffs);
shuff_gp_origs=cell(1,ngps);
nsets_gp=zeros(1,ngps);
for igp=1:ngps
    nsets_gp(igp)=sum(gps==igp);
    shuff_gp_origs{igp}=zeros(1,nsets_gp(igp)); %original group membership in each group
end
nsets=length(gps);
%determine group membership, using logic of hlid_geom_transform_stats
for ishuff=1:nshuffs
    for igp=1:ngps
        gp_select=find(gps(shuffs(ishuff,:))==igp); %shuffled datasets for group igp
        shuff_gp_origs{igp}(ishuff,:)=gps(gp_select);
        for jgp=1:ngps
            mixmats(igp,jgp,ishuff)=sum(shuff_gp_origs{igp}(ishuff,:)==jgp);
        end
    end
end
%for mixing entropy, row sums and column sums are equal to group membership
marg_ent=-sum(lognz(nsets_gp/nsets).*(nsets_gp/nsets));
mixmats_res=reshape(mixmats,ngps*ngps,nshuffs);
table_ents=-sum(lognz(mixmats_res/nsets).*(mixmats_res/nsets),1);
mixent=(table_ents'-marg_ent)/log(2);
return
end
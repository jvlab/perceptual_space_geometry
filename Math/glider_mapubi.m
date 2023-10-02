function [blockcounts,ubi,optsused]=glider_mapubi(map,vecs,g,opts)
% [blockcounts,ubi,optsused]=glider_mapubi(map,vecs,g,opts) counts the unique configurations
%  within arbitrary gliders, and creates a map of indices into them
%  Entries in map must match exactly for two blocks to be considered the same
%
%  map: an image, composed of integers [0:g-1]
%  vecs: an array (nvecs x 2) of displacements that define the glider
%
%  Note that only the relative displacement of vecs (i.e., vecs-vecs(1,:))
%  matter.  vecs should also be row-sorted, since the first pixel specified in 
%  vecs is the least-significant bit in the ordering of the unique blocks
%
%  e.g., for g=3 and two rows in vecs, and assuming vecs is row-sorted
%  unique block 1 has gray level 0 in vecs(1,:), gray level 0 in vecs(2,:)
%  unique block 2 has gray level 1 in vecs(1,:), gray level 0 in vecs(2,:)
%  unique block 3 has gray level 2 in vecs(1,:), gray level 0 in vecs(2,:)
%  unique block 4 has gray level 0 in vecs(1,:), gray level 1 in vecs(2,:), ...
%  unique block 9 has gray level 2 in vecs(1,:), gray level 2 in vecs(2,:)
%
%  g: number of gray levels, defaults to 2
%  opts:
%     opts.mapubi_bc: boundary conditions 0->rectangle (default) or 1->periodic
%  Note that when opts.mapubi_bc=1,gliders are displaced by 1 pixel because
%  of the way that rowmods and colmods are calculated
%
%  blockcounts: array of block counts, one row for each unique block
%     blockcounts is a row vector, of length g^size(vecs,1)
%  ubi:   indices of each pixel into [1:g^nvecs] indicating which kind of configuration
%         is at each location
%   if opts.mapubi_bc=0, then size(ubi)=size(map)-[extremes of vecs]+1
%   if opts.mapubi_bc=1, then size(ubi)=size(map)
%  optsused: options used
%
%  Basic logic is similar to MAPUBI, however
%    * map assumed to be [0:g-1]
%    * general shape allowed for a glider
%    * ubi has one element for all possible blocks, not just those that are present
%    * blockcounts is a single row, corresponding to the last column of an array
%         in mapubi that also specifies the block contents
%
%  See also:  MAPUBI, NARY2INT, GLIDER_MAPUBI_RED, UBI_CODEL_MAP, MAPUBI, UBI_CODEL_MAP_TEST, MLIS_BTCSTATS, BTCSTATS.
%
if (nargin<=2) g=2; end;
if (nargin<=3) opts=[]; end
opts=filldefault(opts,'mapubi_bc',0);
optsused=opts;
%
%
minvec=min(vecs,[],1);
maxvec=max(vecs,[],1);
block=maxvec-minvec+1;
%
if (opts.mapubi_bc==1)
    rowmods=1+mod([1:(size(map,1)+block(1)-1)],size(map,1));
    colmods=1+mod([1:(size(map,2)+block(2)-1)],size(map,2));
    map=[map(rowmods,colmods)];
end
%
nvecs=size(vecs,1);
nrc=size(map)-block+1;
ri=[1:nrc(1)];
ci=[1:nrc(2)];
% make an exhaustive list of blocks and their counts
blocks=zeros([nrc nvecs]);
for ivec=1:nvecs
    blocks(:,:,ivec)=map(ri+vecs(ivec,1)-minvec(1),ci+vecs(ivec,2)-minvec(2));
end %ir
ubi=1+nary2int(blocks,g,3);
blockcounts=histc(ubi(:),[1:g^nvecs]);
blockcounts=blockcounts(:)'; %force a row
return

function [eivals,coords,opts_used]=psg_isomap(dists,opts)
% [eivals,coords,opts_used]=psg_isomap(dists,opts) finds eigenvalues and coordinates of MDS solution 
% of isomap (and standard) embeddings
%
% A Global Geometric Framework for Nonlinear Dimensionality Reduction
% Tenenbaum, de Silva, and Langford, Science 290, p. 2319 (2000)
%
%  dists: [nstims, nstims]: symmetric array of distances
%  opts: 
%    opts.min_nn: minimum number of neatest neighbors, defaults to 1, 0 does standard multidimensional scaling
%    opts.if_log: 1 for min_nbr_conn to log results, defaults to 0
%
%  eivals: a vector of eigenvalues of isomap embedding.  
%    Negative values correspond to hyperbolic coordinates. 
%    Given in descending order of absolute value
%  coords: coordinates.  sqrt(abs(eivals)) has already been applied
%    Note that if eivals(k) is negative, the coordinate needs to be interpreted as a hyperbolic coord
%  opts_used: options used, includes intermdiate cacluations if opts.min_nn>0:
%   nbr_mindist: minimum distance that gives each node nn neighbors (if min_nn>0)
%   nbr_mtx: graph distance on nearest-neighbor graph (0: no connection)
%   nbr_graph: graph structure for nearest-neighbor graph
%   isomap_dists: distances on nearest-neighbor graph
%
% Requires graph toolkit.
%
% See also: PSG_ISOMAP_DEMO, DOMDS, MIN_NBR_CONN, DISTANCES, COOTODSQ
if (nargin<=1)
    opts=struct;
end
opts=filldefault(opts,'min_nn',1);
opts=filldefault(opts,'if_log',0);
opts_nn=opts;
if opts.min_nn>0
    [nbr_mindist,nbr_mtx,nbr_graph,opts_used]=min_nbr_conn(dists,opts);
    isomap_dists=distances(nbr_graph);
    opts.nbr_mindist=nbr_mindist;
    opts.nbr_mtx=nbr_mtx;
    opts.nbr_graph=nbr_graph;
    [eivals_raw,eivecs_raw]=domds(isomap_dists,1); % [eival,eivec]=domds(distmtx,p) does a multidimensional scaling
else
    [eivals_raw,eivecs_raw]=domds(dists,1);   
end
opts_used=opts;
%sort by descending abs value
[eivals_sort,ix]=sort(abs(eivals_raw),'descend');
eivals_raw=eivals_raw(ix);
eivecs_raw=eivecs_raw(:,ix);
eivals=real(eivals_raw)/2; %per do_mds, the eigenvalues need to be divided by 2, prior to square root, to obtain coords that recapitulte the distances.
coords=eivecs_raw.*repmat(sqrt(abs(eivals')),[size(eivecs_raw,1),1]); %scale by eivals
return

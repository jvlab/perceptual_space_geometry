function [rmtx,opts_used]=incidence2rigidity(imtx,d,opts)
% [rmtx,opts_used]=incidence2rigidity(imtx,d,opts) computes a rigidity matrix from an incidence matrix using random edge lengths
%
% See rigidity_notes.docx,
% Wikipedia https://en.wikipedia.org/wiki/Rigidity_matroid
% Wolfram Mathworld, https://mathworld.wolfram.com/RigidGraph.html
%
% The corresponding structure is infinitesimally rigid iff rank(rmtx)=2*number of vertices-3 (Mathworld), quoting Grasseger (2023),
% but this ref appears to apply only to dimension 2
%
% imtx: incidence matrix for a non-directed graph, each row is a vertex, each column is an edge; two elements equal to 1 in each row, rest are zero
% d: dimension in which graph is embedded
% opts: options (not yet used)
%
% rmtx: rigidity matrix
% opts_used: options used
%
%  See also: RIGIDITY_COMPLETE_DEMO, RIGIDITY_CONSTRAINTS.
%
if (nargin<=2)
    opts=struct();
end
edges=size(imtx,2);
vtcs=size(imtx,1);
rmtx=zeros(edges,vtcs*d);
%generate random coordinates
c=randn(vtcs,d);
for ie=1:edges
    v=find(imtx(:,ie)>0);
    z=c(v(1),:)-c(v(2),:);
    rmtx(ie,d*(v(1)-1)+[1:d])=z;
    rmtx(ie,d*(v(2)-1)+[1:d])=-z;
end
opts_used=opts;
return

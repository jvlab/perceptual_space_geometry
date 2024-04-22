function d=ord_char_simul_dist(type,pts,S)
% d=ord_char_simul_dist(type,pts,coords,ranks,type)
% 
% computes any of several distances in a toy model of a perceptual space
%
% type: one of 'ultra','addtree','addtree_weighted','euclidean'
% pts: indices of two points, must be 1:size(coords,2)
% S: a structure that describes the domain, with 
%   S.coords: coordinates of the points
%   S.edges: points that are paired
%   S.edgelengths: pairwise distances for each edge
%   S.ranks: ranks of the points (for a hierarchy)
%   S.graph: corresponding graph structure, must be a tree 
%
% d: the distance
%
% if called with no arguments, returns the available options for distances
%
% uses shortestpath in graph toolbox
% 
% See also: ORD_CHAR_SIMUL_DEMO.
%
if (nargin==0)
    d={'ultra','addtree','addtree_weighted','euclidean','line','ring_as_graph'};
    return
else
    d=[];
    if pts(1)==pts(2)
        d=0;
    else
        switch type
            case 'ultra'
                 path=shortestpath(S.graph,pts(1),pts(2));
                 d=max(S.ranks(path));
            case 'addtree'
                path=shortestpath(S.graph,pts(1),pts(2));
                d=length(path);
            case 'addtree_weighted'
                [path,d]=shortestpath(S.graph,pts(1),pts(2));
            case 'euclidean'
                d=sqrt(sum((S.coords(pts(1),:)-S.coords(pts(2),:)).^2));
            case 'line'
                d=abs(diff(pts));
            case 'ring_as_graph'
                d=abs(diff(pts));
                d=min(d,S.npts-d);
        end
    end
end
return

function [d,adj_model,transform,opts_used]=psg_geo_procrustes(ref,adj,opts)
% [d,adj_model,transform,opts_used]=psg_geo_procrustes(ref,adj,opts) finds a procrustes model
% with standardized input and output variables
%
% ref: reference coordinates, size=[npts,ref_dim]
% adj: coordinates of dataset to be adjusted, size=[npts,adj_dim]
% opts: options, can be omitted or empty
%     opts.if_scale: 1 to allow scaling (default: 0)
%
% d: goodness of fit, mean squared deviation of adj_model from res, divided by mean squared dev of ref
% adj_model: coordinates of adj, after transformation
% transform:
%   transform.b: scalar, will be 1 if opts.if_scale=0
%   transform.T: matrix of size [adj_dim max(ref_dim,adj_dim)]
%   transform.c: offset, row of size [1 max(ref_dim, adj_dim)], will be zeros if ref and adj are centered
%   adj_model=transform.b*adj*transform.T+repmat(transform.c,npts,1)
% opts_used: options used
%
%   See also: PROCRUSTES, PSG_GEO_AFFINE, PSG_GEOMODELS_TEST, PSG_GEOMODELS_DEFINE, PSG_GEO_GENERAL,
%     PSG_GET_TRANSFORM.
%
if (nargin<=2) opts=struct; end
opts=filldefault(opts,'if_scale',0);
adj_dim=size(adj,2);
ref_dim=size(ref,2);
opts_used=opts;
%
npts=size(adj,1);
if ref_dim<adj_dim
    ref=[ref,zeros(npts,adj_dim-ref_dim)];
end
%
if opts.if_scale==0
    Scaling=false;
else
    Scaling=true;
end
%
% Procrustes
%
[d,adj_model,transform]=procrustes(ref,adj,'Scaling',Scaling);
transform.c=transform.c(1,:); %clean up redundant rows
return
end

function [d,transform]=psg_geo_pwaffine_obj(params,ref,adj,u,a,opts)
% [d,transform]=psg_geo_pwaffine_obj(params,ref,adj,u,a,opts) computes the objective
% function for a piecewise affine model, with unconstrained parameters 
% in the tangent plane (see psg_geo_pwaffine)
%
% params: vector of length [1 dim_x]
% ref: reference dataset, [npts dim_xy], dim_xy >= dim_x
% adj: dataset to adjust, [npts dim_x]
% u: coordinate frame for parameterization, size [dim_x dim_x]
%    See psg_geo_pwaffine_va for relationship of u to the cutplanes.
% a: reference value for cutpoint
% opts: options for psg_geo_pwaffine_va
%
% d: residuals
% transform: transform of best piecewise model (has acut, vcut)
%
%   See also:  PSG_GEO_PWAFFINE, PSG_GEO_PWAFFINE_VA.
%
if (nargin<=5)
    opts=struct;
end
v=params*u+u(1,:); %vector in the tangent plane
vnorm=sqrt(sum(v.^2));
vcut=v/vnorm;
acut=a*vnorm;
[d,transform]=psg_geo_pwaffine_va(ref,adj,vcut,acut,opts);
end

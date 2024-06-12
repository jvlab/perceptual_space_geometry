function coords_new=psg_geomodels_apply(model_class,coords,transform)
% coords_new=psg_geomodels_apply(model_class,coords,transform) applies any of severalsurf
% geometric models to a coordinate set
%
% model_class: one of the model classes in psg_geomodels_define, typically one of 
% 'mean','procrustes','affine','projective','pwaffine'
% coords: coordinates: size is [npts,ndims]
% transform: a structure containing the transform parameters, typically
%    transform.b: scale
%    transform.T: [ndims,ndims]: a post-multiplying matrix
%    transorms.c: offset, [1 ndims]
%    also may include a cutpoint, or a projective parameter
%
% coords_new: coordinates after transformation, size is [npts,ndims]
%
%   See also:  PSG_GEOMODELS_RUN, PROCRUSTES, PSG_GEO_GENERAL, PSG_GEOMODELS_DEFINE,
%     PERSP_APPLY, PSG_PWAFFINE_APPLY.
%
switch model_class
    case {'mean','procrustes','affine'}
       npts=size(coords,1);
       coords_new=transform.b*coords*transform.T+repmat(transform.c,npts,1);
    case 'projective'
        %note change of variable names between psg_geo_projective and persp_apply
        coords_new=persp_apply(transform.T,transform.c,transform.p,coords);
    case 'pwaffine'
        coords_new=psg_pwaffine_apply(transform,coords);
    otherwise
        warning(sprintf('unknown model class %s',model_class));
        coords_new=nan(size(coords));
end
return

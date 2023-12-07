function model_types_def=psg_geomodels_define()
%model_types_def=psg_geomodels_define() sets up the definitions of geometric model types
%
% model_types_def: a structure that defines the models and their hierarchical relationships
%
%   See also: PSG_GEOMODELS_TEST, PSG_GEO_GENERAL.
%
model_types_def=struct;
model_types_def.model_types={'mean','procrustes_noscale','procrustes_scale','affine_nooffset','affine_offset','projective'};
%
model_types_def.mean.opts=struct;
model_types_def.mean.nested={};
%
model_types_def.procrustes_noscale.opts.if_scale=0;
model_types_def.procrustes_noscale.nested={'mean'};
%
model_types_def.procrustes_scale.opts.if_scale=1;
model_types_def.procrustes_scale.nested={'mean','procrustes_noscale'};
%
model_types_def.affine_nooffset.opts.if_offset=0;
model_types_def.affine_nooffset.nested={'mean','procrustes_noscale','procrustes_scale'};
%
model_types_def.affine_offset.opts.if_offset=1;
model_types_def.affine_offset.nested={'mean','procrustes_noscale','procrustes_scale'};
%
model_types_def.projective.opts.method='fmin';
model_types_def.projective.nested={'mean','procrustes_noscale','procrustes_scale','affine_offset'};
%
mnames=model_types_def.model_types;
for im=1:length(mnames)
    mname=mnames{im};
    mclass=mname(1:-1+min(find(cat(2,mname,'_')=='_')));
    model_types_def.(mname).class=mclass;
end
return
end

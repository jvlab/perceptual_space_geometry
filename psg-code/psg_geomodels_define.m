function model_types_def=psg_geomodels_define(if_select)
%model_types_def=psg_geomodels_define(if_select) sets up the definitions of geometric model types
%
% if_select: if present and nonzero, can select model types
% model_types_def: a structure that defines the models and their hierarchical relationships
%
% 05Dec23: add piecewise affine
% 21Dec23: add dof 
%    dof is an array whose entry in row (i+1) and column (j+1) is the
%    coefficient of (nx)^i*(ny)^j in a polynomial that gives the number of degrees of
%    freedom in a model, where nx is the number of dimensions in input space for the transformation
%    and ny is the number of dimensions in the output space, and ny>=nx.
%
%   See also: PSG_GEOMODELS_TEST, PSG_GEO_GENERAL, PSG_GEOMODELS_ILLUS, PSG_GEO_PWAFFINE_TEST,
%     PSG_GEOMODELS_NDOF.
%
% 28May24: add if_seloect; add min_inputdims, minimum number of dimensions in input space required to fit

if nargin==0
    if_select=0;
end
model_types_def=struct;
model_types_def.model_types={'mean','procrustes_noscale','procrustes_scale','affine_nooffset',...
    'affine_offset','projective','pwaffine','pwaffine_2'};
%
model_types_def.mean.opts=struct;
model_types_def.mean.nested={};
model_types_def.mean.dof=[0];
%
model_types_def.procrustes_noscale.opts.if_scale=0;
model_types_def.procrustes_noscale.nested={'mean'};
model_types_def.procrustes_noscale.dof=[0 0;-1 2;-1 0]/2;  %[(ny-(nx+1)/2)*nx]
%
model_types_def.procrustes_scale.opts.if_scale=1;
model_types_def.procrustes_scale.nested={'mean','procrustes_noscale'};
model_types_def.procrustes_scale.dof=model_types_def.procrustes_noscale.dof+[1 0;0 0; 0 0]; %add 1
%
model_types_def.affine_nooffset.opts.if_offset=0;
model_types_def.affine_nooffset.nested={'mean','procrustes_noscale','procrustes_scale'};
model_types_def.affine_nooffset.dof=[0 0;0 1];%  ny*nx (a general matrix)
%
model_types_def.affine_offset.opts.if_offset=1;
model_types_def.affine_offset.nested={'mean','procrustes_noscale','procrustes_scale'};
model_types_def.affine_offset.dof=[0 1;0 1];%  ny*nx+ny
%
model_types_def.projective.opts.method='fmin';
model_types_def.projective.opts.if_display=1;
model_types_def.projective.nested={'mean','procrustes_noscale','procrustes_scale','affine_offset'};
model_types_def.projective.dof=[0 0;1 1];% ny*nx+nx % a general projective matrix but excluding offset ny)
%
model_types_def.pwaffine.opts.method='fmin';
model_types_def.pwaffine.opts.if_display=1;
model_types_def.pwaffine.opts.n_cuts_model=1;
model_types_def.pwaffine.nested={'mean','procrustes_noscale','procrustes_scale','affine_offset'};
model_types_def.pwaffine.dof=[0 1;0 1];% ny*nx+ny
%
%
model_types_def.pwaffine_2.opts.method='fmin';
model_types_def.pwaffine_2.opts.if_display=1;
model_types_def.pwaffine_2.opts.n_cuts_model=2;
model_types_def.pwaffine_2.nested={'mean','procrustes_noscale','procrustes_scale','affine_offset','pwaffine'};
model_types_def.pwaffine_2.dof=[0 2;0 1];% ny*nx+2*ny
model_types_def.pwaffine_2.min_inputdims=2;
%
mnames=model_types_def.model_types;
for im=1:length(mnames)
    mname=mnames{im};
    mclass=mname(1:-1+min(find(cat(2,mname,'_')=='_')));
    model_types_def.(mname).class=mclass;
    model_types_def.(mname)=filldefault(model_types_def.(mname),'min_inputdims',1);
end
mnames_orig=mnames;
if (if_select)
    for im=1:length(mnames)
        disp(sprintf('%1.0f->%s',im,mnames{im}))
    end
    remove_list=getinp('list of models to remove','d',[0 length(mnames)],0);
    remove_list=setdiff(remove_list,0);
    for irem=1:length(remove_list)
        im=remove_list(irem);
        remove_name=mnames_orig{im};
        model_types_def=rmfield(model_types_def,remove_name);
        %remove model from list of model types
        ir=strmatch(remove_name,model_types_def.model_types,'exact');
        model_types_def.model_types={model_types_def.model_types{1,setdiff(1:length(model_types_def.model_types),ir)}};
        %remove from nesteds
        mnames=model_types_def.model_types;
        for im=1:length(mnames)
            mname=mnames{im};
            ir=strmatch(remove_name,model_types_def.(mname).nested,'exact');
            if ir>0
               model_types_def.(mname).nested={model_types_def.(mname).nested{1,setdiff(1:length(model_types_def.(mname).nested),ir)}};
            end
        end
    end
end
return
end

function ndof=psg_geomodels_ndof(mdef,nx,ny,model_type)
% ndof=psg_geomodels_ndof(mdef,nx,ny,model_type) returns the number of degrees of
% freedom for a geometric model specified by psg_geomodels_define.
%
% mdef: the model definition structure, from psg_geomodels_define, may be
%    empty
% nx: dimension of the input dataset (the set to be adjusted)
% ny: dimension of the target dataset (the reference)
%   typically want ny>=nx
% model_type: string, one of the fields of mdef
%   if model_type is omitted or empty, then ndof is a structure, with one entry
%   for each available model
%
% ndof: number of degrees of freedom
%   This does not include translation of the output, except for affine_offset.
%   For other model types, to include translation of the output, add ny to ndof
%
%    See also: PSG_GEOMODELS_DEFINE.
%
if isempty(mdef)
    mdef=psg_geomodels_define;
end
if nargin<=3
    model_type=[];
else
    if ~isfield(mdef,model_type)
       warning(sprintf('model type %s not recognized',model_type));
       model_type=[];
    end
end
if isempty(model_type)
    model_types=fieldnames(mdef);
    ndof=struct;
    for im=1:length(model_types)
        if isfield(mdef.(model_types{im}),'class')
            ndof.(model_types{im})=psg_geomodels_ndof_do(mdef.(model_types{im}),nx,ny);
        end
    end
else
    ndof=psg_geomodels_ndof_do(mdef.(model_type),nx,ny);
end
return

function np=psg_geomodels_ndof_do(md,nx,ny)
if ~isfield(md,'dof')
    np=0;
    return
else
    dof=md.dof;
    vx=nx.^[0:size(dof,1)-1]';
    vy=ny.^[0:size(dof,2)-1];
    vxy=vx*vy; %monomials
    np=sum(sum(dof.*vxy));
end
return



%psg_geo_pwaffine_test: test various strategies for piecewise-affine fit
%
%  See also:  PSG_PWAFFINE_APPLY, PSG_GEOMODELS_DEFINE, PSG_GEO_LAYOUTS_SETUP, PSG_GEO_TRANSFORMS_SETUP.
%
if ~exist('dim_max') dim_max=3; end %max dimension for layouts
%
%model fitting params
model_types_def=psg_geomodels_define();
model_types=model_types_def.model_types;
nmodels=length(model_types);
if ~exist('layouts')
    layouts=psg_geo_layouts_setup(dim_max);
end
nlayouts=length(layouts);
%
%set up transforms
%
if ~exist('transforms')
    transforms=psg_geo_transforms_setup(dim_max);
end
ntransforms=length(transforms);
%
%set up structures needed for fitting data
%
ds=cell(nlayouts,1);
for il=1:nlayouts
    for id=1:dim_max
        ds{il}{id}=layouts{il}.coords(:,1:id);
    end
    disp(sprintf(' layout %2.0f (%3.0f points):  %s',il,layouts{il}.npts,layouts{il}.label));
end
%
if ~exist('layout_list') layout_list=[1 4:6]; end
if ~exist('transform_list') transform_list=6; end %piecewise linear
layout_list=getinp('layouts to use','d',[1 nlayouts],layout_list);
for it=1:ntransforms
    disp(sprintf(' %2.0f->%s',it,transforms{it}.label));
end
transform_list=getinp('transforms to simulate','d',[1 ntransforms],transform_list);
%
%compute and plot transforms and fits
transform_fit=cell(ntransforms,nlayouts);
coords_fit=cell(ntransforms,nlayouts);
d_fit=zeros(ntransforms,nlayouts);
for ilptr=1:length(layout_list)
    il=layout_list(ilptr);
    disp(' ');
    disp(sprintf('layout %s',layouts{il}.label));
    %
    for itptr=1:length(transform_list)
        it=transform_list(itptr);
        coords_orig=ds{il}{dim_max};
        switch transforms{it}.class
            case {'mean','procrustes','affine'}
                coords_transformed=transforms{it}.params.b*coords_orig*transforms{it}.params.T+repmat(transforms{it}.params.c,layouts{il}.npts,1);
            case 'projective'
                coords_transformed=persp_apply(transforms{it}.params.T,transforms{it}.params.c,transforms{it}.params.p,coords_orig);
            case 'pwaffine'
                coords_transformed=psg_pwaffine_apply(transforms{it}.params,coords_orig);
            otherwise
                coords_transformed=orig;
        end %model class
        model_type='pwaffine';
        opts_model=model_types_def.(model_type).opts;
        model_class=model_types_def.(model_type).class;
        opts_model=model_types_def.(model_type).opts;
        [d_fit(it,il),coords_fit{it,il},transform_fit{it,il},opts_model_used]=psg_geo_general(coords_transformed,coords_orig,model_class,opts_model);
        disp(sprintf(' transforming with %30s, fitting with model type %30s: d=%8.5f',transforms{it}.label,model_type,d_fit(it,il)));
        %
    end %itptr
end %il

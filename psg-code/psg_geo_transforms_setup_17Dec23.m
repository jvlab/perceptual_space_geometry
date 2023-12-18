function transforms=psg_geo_transforms_setup(dim_max)
% transforms=psg_geo_transforms_setup(dim_max) creates some illustrative perceptual space transformations
% dim_max: dimension of coordinates
% transforms: a cell array of transformation specifications
%    transforms{il}.label: text label
%    transforms{il}.model_type: model type, psg_geomodels_define
%    transforms{il}.class: model type, psg_geomodels_define
%    transforms{il}.params: parameters of transformation: b (scale), T (linear transform), c (offset)
%
%See also: PSG_GEOMODELS_ILLUS, PSG_GEO_LAYOUTS_SETUP, PSG_GEOMODELS_DEFINE, PSG_GEO_PWAFFINE_TEST.
%
transforms=cell(0);
rotang12=pi/3;
rotang13=pi/10;
rot12=[cos(rotang12) sin(rotang12);-sin(rotang12) cos(rotang12)];
rot13=[cos(rotang13) sin(rotang13);-sin(rotang13) cos(rotang13)];
%
rot=eye(dim_max);
rot([1 2],[1 2])=rot12;
if dim_max>=3
    rot3=eye(dim_max);
    rot3([1 3],[1 3])=rot13;
    rot=rot*rot3;
end
aff=eye(dim_max);
proj=zeros(dim_max,1);
vcut=zeros(1,dim_max);
tdif=zeros(1,dim_max);
%   quantities specified up to 3d, extend beyond with zeros
affratios=[0.8 1.5 1.2 0.7];
projvec=[0.1 0.07 -0.3]';
vcutvec=[.2 .8 .1];
vcutvec=vcutvec./sqrt(sum(vcutvec.^2));
tdifvec=[.5 0 1];
for id=1:min(3,dim_max)
    aff(id,id)=affratios(id);
    proj(id)=projvec(id);
    vcut(id)=vcutvec(id);
    tdif(id)=tdifvec(id);
end
%
transforms{1}.label='original';
transforms{1}.model_type='procrustes_noscale';
transforms{1}.class='procrustes';
transforms{1}.params.T=eye(dim_max);
transforms{1}.params.b=1;
transforms{1}.params.c=zeros(1,dim_max);
%
transforms{2}=transforms{1};
transforms{2}.label='procrustes_noscale';
transforms{2}.params.T=rot;
%
transforms{3}=transforms{2};
transforms{3}.label='procrustes_scale';
transforms{3}.model_type='procrustes_scale';
transforms{3}.params.b=1.25;
%
transforms{4}=transforms{2};
transforms{4}.label='affine_nooffset';
transforms{4}.class='affine';
transforms{4}.model_type='affine_nooffset';
transforms{4}.params.T=aff*rot;
%
transforms{5}=transforms{4};
transforms{5}.label='projective';
transforms{5}.class='projective';
transforms{5}.model_type='projective';
transforms{5}.params.p=proj;
%
transforms{6}=transforms{4};
transforms{6}.label='piecewise affine';
transforms{6}.class='pwaffine';
transforms{6}.model_type='pwaffine';
transforms{6}.params.vcut=vcut;
transforms{6}.params.acut=-0.5;
transforms{6}.params.T=repmat(transforms{6}.params.T,[1 1 2]);
transforms{6}.params.T(:,:,1)=transforms{6}.params.T(:,:,2)+vcut'*tdif; %will do nothing if orthog to vcut
transforms{6}.params.c=zeros(2,dim_max);
%
ntransforms=length(transforms);
for it=1:ntransforms
    transforms{it}=filldefault(transforms{it},'model_type',transforms{it}.label);
end
return

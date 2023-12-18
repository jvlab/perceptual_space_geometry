function transforms=psg_geo_transforms_setup(dim_max)
% transforms=psg_geo_transforms_setup(dim_max) creates some illustrative perceptual space transformations
% dim_max: dimension of coordinates
% transforms: a cell array of transformation specifications
%    transforms{il}.label: text label
%    transforms{il}.model_type: model type, psg_geomodels_define
%    transforms{il}.class: model type, psg_geomodels_define
%    transforms{il}.params: parameters of transformation: b (scale), T (linear transform), c (offset)
%
% 18Dec23:  fix discontinuity with piecewise affine, and add second example
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
vcut2=zeros(1,dim_max);
tdif2=zeros(1,dim_max);
acuts=[-0.50 1.32];
%   quantities specified up to 3d, extend beyond with zeros
affratios=[0.8   1.50  1.20 0.70];
projvec=  [0.10  0.07 -0.30]';
vcutvec=  [0.20  0.80  0.10];
tdifvec=  [0.50  0.00  1.00];
%
vcutvec2= [0.70 -0.40 -0.20];
tdifvec2= [0.40 -0.25  0.65];
%
vcutvec=vcutvec./sqrt(sum(vcutvec.^2));
vcutvec2=vcutvec2./sqrt(sum(vcutvec2.^2));
for id=1:min(3,dim_max)
    aff(id,id)=affratios(id);
    proj(id)=projvec(id);
    vcut(id)=vcutvec(id);
    tdif(id)=tdifvec(id);
    vcut2(id)=vcutvec2(id);
    tdif2(id)=tdifvec2(id);
end
Tbase=aff*rot;
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
transforms{4}.params.T=Tbase;
%
transforms{5}=transforms{4};
transforms{5}.label='projective';
transforms{5}.class='projective';
transforms{5}.model_type='projective';
transforms{5}.params.p=proj;
%
%create piecewise affine that agrees on cutplane
%
transforms{6}.label='piecewise affine';
transforms{6}.class='pwaffine';
transforms{6}.model_type='pwaffine';
transforms{6}.params.vcut=vcut;
transforms{6}.params.acut=acuts(1);
T=cat(3,Tbase+vcut'*tdif,Tbase);
transforms{6}.params.T=T;
transforms{6}.params.c=psg_geo_transforms_getc(dim_max,T,vcut,acuts(1));
%
transforms{7}.label='piecewise affine B';
transforms{7}.class='pwaffine';
transforms{7}.model_type='pwaffine';
transforms{7}.params.vcut=vcut2;
transforms{7}.params.acut=acuts(2);
T=cat(3,Tbase+vcut2'*tdif2,Tbase);
transforms{7}.params.T=T;
transforms{7}.params.c=psg_geo_transforms_getc(dim_max,T,vcut2,acuts(2));
%
ntransforms=length(transforms);
for it=1:ntransforms
    transforms{it}=filldefault(transforms{it},'model_type',transforms{it}.label);
    transforms{it}.params=filldefault(transforms{it}.params,'b',1);
    transforms{it}.params=filldefault(transforms{it}.params,'c',zeros(size(transforms{it}.params.T,3),dim_max));
end
return

function c=psg_geo_transforms_getc(dim_max,T,vcut,acut)
% c=psg_geo_transforms_getc(dim_max,T,vcut,acut) is a utiltiy to adjust offsets so that 
% piecewise linear transfomrs match on the cut hyperplane;
vtest=[1:dim_max]==ones(1,dim_max); %a vector [1 0 0 0 ..0]
vorth=vtest-vtest*vcut'*vcut; %orthogonal to vcut
vtry=vorth+vcut*acut; %on the cut plane
c=-[vtry*T(:,:,1);vtry*T(:,:,2)];
return

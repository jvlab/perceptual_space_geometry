function layouts=psg_geo_layouts_setup(dim_max)
% layouts=psg_geo_layouts_setup(dim_max) creates some illustrative layouts for perceptual spaces
% dim_max: dimension of coordinates
% layouts: a cell array of stimulus setups
%    layouts{il}.npts: number of points
%    layouts{il}.label: text label
%    layouts{il}.coords: stimuls coordinates, siz is [npts dim_max]
%
%See also: PSG_GEOMODELS_ILLUS, PSG_GEO_TRANSFORMS_SETUP, PSG_GEO_PWAFFINE_TEST.
%
layouts=cell(0);
nxy=[5 5];
layouts{1}.label=sprintf('%1.0f x %1.0f grid',nxy);
[xg,yg]=meshgrid([1:nxy(1)],[1:nxy(2)]);
layouts{1}.coords=[xg(:),yg(:)];
%
nxy=[5 7];
layouts{2}.label=sprintf('%1.0f x %1.0f grid',nxy);
[xg,yg]=meshgrid([1:nxy(1)],[1:nxy(2)]);
layouts{2}.coords=[xg(:),yg(:)];
%
nrays=8;
neachray=3;
layouts{3}.label=sprintf('%1.0f rays, %1.0f points each',nrays,neachray);
angs=2*pi*floor([0:nrays*neachray-1]/neachray)/nrays;
rads=1+mod([0:nrays*neachray-1],neachray);
layouts{3}.coords=[[0; cos(angs(:)).*rads(:)],[0; sin(angs(:)).*rads(:)]];
layouts{3}.keeprays=1;
%
nrays=4;
neachray=6;
layouts{4}.label=sprintf('%1.0f rays, %1.0f points each',nrays,neachray);
angs=2*pi*floor([0:nrays*neachray-1]/neachray)/nrays;
rads=1+mod([0:nrays*neachray-1],neachray);
layouts{4}.coords=[[0; cos(angs(:)).*rads(:)],[0; sin(angs(:)).*rads(:)]];
%
ncirc=24;
rads=3.;
layouts{5}.label=sprintf('%2.0f point circle, radius %2.0f',ncirc,rads);
angs=2.*pi*[0:ncirc-1]/ncirc;
layouts{5}.coords=[[0; cos(angs(:)).*rads],[0; sin(angs(:)).*rads]];
%
nxy=[11 11];
layouts{6}.label=sprintf('%1.0f x %1.0f grid (0.5 side length)',nxy);
[xg,yg]=meshgrid([1:nxy(1)],[1:nxy(2)]);
layouts{6}.coords=[xg(:),yg(:)]/2;
%add extra dimensions, subtract mean, add needed fields
for il=1:length(layouts)
    layouts{il}.npts=size(layouts{il}.coords,1);
    if size(layouts{il}.coords,2)<dim_max
        layouts{il}.coords=[layouts{il}.coords,zeros(layouts{il}.npts,dim_max-size(layouts{il}.coords,2))];
    end
    layouts{il}.coords=layouts{il}.coords-repmat(mean(layouts{il}.coords,1),layouts{il}.npts,1);
    layouts{il}=filldefault(layouts{il},'keeprays',0);
end
return
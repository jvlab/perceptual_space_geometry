% irgb_gzranking_getimages: one-time script to get the Giesel-Zaidi images
%
% unzip GizelZaidiImages.zip into TempIages, and run this once to create irgb_gzimages.mat, which
% contains all 261 images as full-color, cut down to the 150 x 150 stimulus size.
%
%   See also:  IRGB_GZRANKING_READ, IRGB_GZRANKING_ANALYZE.
%
if ~exist('filename_rank_def') filename_rank_def='C:/Users/jdvicto/Documents/jv/ENCL/Zaidi/RankingData.mat'; end
if ~exist('filename_img_def') filename_img_def='C:/Users/jdvicto/Documents/jv/ENCL/Zaidi/ImageNames.mat'; end
if ~exist('path_images_def') path_images_def='C:/Users/jdvicto/Documents/jv/ENCL/Zaidi/TempImages/'; end
if ~exist('imagedata_fn_def') imagedata_fn_def='./irgb_gzranking_imagedata.mat';end
%
filename_rank=getinp('Giesel-Zaidi ranking data file','s',[],filename_rank_def);
filename_img=getinp('Giesel-Zaidi image name file','s',[],filename_img_def);
path_img=getinp('path to images','s',[],path_images_def);
%
s=irgb_gzranking_read(filename_rank,filename_img);
%
image_data=struct;
nimgs=length(s.image_names_all_short);
image_data.names=s.image_names_all_short;
for k=1:nimgs
    fn=image_data.names{k};
    fn=fn(1+max(find(fn=='/')):end);
    image_data.rgbvals(:,:,:,k)=imread(cat(2,path_img,fn));
end
if (getinp('1 to save','d',[0 1]));
    imagedata_fn=getinp('file name','s',[],imagedata_fn_def);
    save(imagedata_fn,'image_data');
    disp('saved.');
end

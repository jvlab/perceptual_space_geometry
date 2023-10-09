% irgb_gzranking_getimgs: one-time script to get the Giesel-Zaidi images
%
% unzip GizelZaidiImages.zip OR Ranking4Qasim, subdirectoryies with "cut",
% into TempIages, and run this once to create irgb_gzimages.mat, which
% contains all 261 images as full-color, cut down to the 150 x 150 stimulus size.
%
%   See also:  IRGB_GZRANKING_READ, IRGB_GZRANKING_ANALYZE.
%
if ~exist('filename_rank_def') filename_rank_def='C:/Users/jdvicto/Documents/jv/ENCL/Zaidi/RankingData.mat'; end
if ~exist('filename_img_def') filename_img_def='C:/Users/jdvicto/Documents/jv/ENCL/Zaidi/ImageNames.mat'; end
if ~exist('path_images_def') path_tempimages_def='C:/Users/jdvicto/Documents/ENCL/Zaidi/TempImages/'; end
    %
filename_rank=getinp('Giesel-Zaidi ranking data file','s',[],filename_rank_def);
filename_img=getinp('Giesel-Zaidi image name file','s',[],filename_img_def);
path_img=getinp('path to images','s',[],path_images_def);
%
s=irgb_gzranking_read(filename_rank,filename_img);
%

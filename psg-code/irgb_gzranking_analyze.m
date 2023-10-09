% irgb_gzranking_analyze: analyzes the Giesel-Zaidi ranking data
%
%   See also:  IRGB_GZRANKING_READ.
%
if ~exist('filename_rank_def') filename_rank_def='C:/Users/jdvicto/Documents/jv/ENCL/Zaidi/RankingData.mat'; end
if ~exist('filename_img_def') filename_img_def='C:/Users/jdvicto/Documents/jv/ENCL/Zaidi/ImageNames.mat'; end
%
filename_rank=getinp('Giesel-Zaidi ranking data file','s',[],filename_rank_def);
filename_img=getinp('Giesel-Zaidi imgage name file','s',[],filename_img_def);
%
s=irgb_gzranking_read(filename_rank,filename_img);
%
figure;
set(gcf,'Position',[100 100 1200 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name','G-Z ranking data');
%
n_props=length(s.props);
n_subjs=length(s.subjs);
%
crange=[min(s.rankings(:)),max(s.rankings(:))];
    for iprop=1:n_props
    subplot(1,n_props,iprop);
    imagesc(s.rankings(:,:,iprop),crange);
    colorbar;
    title(s.props{iprop});
    set(gca,'XTickLabel',s.subjs);
    ylabel('rated image index');
    xlabel('subject ID');
end

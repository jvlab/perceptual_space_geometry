% irgb_gzranking_analyze: analyzes the Giesel-Zaidi ranking data
%
%   See also:  IRGB_GZRANKING_READ.
%
if ~exist('filename_def') filename_def='C:/Users/jdvicto/Documents/jv/ENCL/Zaidi/RankingData.mat'; end
%
filename=getinp('Giesel-Zaidi ranking data file','s',[],filename_def);
s=irgb_gzranking_read(filename_def);
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

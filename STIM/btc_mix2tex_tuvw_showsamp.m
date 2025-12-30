% btc_mix2tex_tuvw_showsamp: special-purpose script to show a sample texture made with t,u,v,w specified
%
%   See also: BTDC_MIX2TEX_DEMO.
%
if ~exist('sampsize' sampsize=32; end
figure;
imagesc(map([1:sampsize],[1:sampsize]),[0 1]);
hold on;
colormap gray;axis equal;axis tight;
set(gca,'XTick',[]);
set(gca,'YTick',[]);
title(sprintf('t %5.2f u %5.2f v %5.2f w %5.2f',mean(tvec(:,6:9))))

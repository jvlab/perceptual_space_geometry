%psg_isomap_demo_animals_heatmap
%heatmaps of participation ratios
figure;
set(gcf,'Position',[100 100 1200 800]);
dim_list=[2:7];
%
subplot(2,3,1);
load('psg_isomap_demo_animals_texture.mat')
imagesc(median(part_ratio(dim_list,1+[0 nn_list],:),3)',[1 7]);
set(gca,'YTickLabel',strvcat('std',num2str(nn_list(:))));
set(gca,'XTick',[1:length(dim_list)]);
set(gca,'XTickLabel',dim_list);
xlabel('dim');
ylabel('nn');
title('median participation ratio, texture')
colorbar
%
subplot(2,3,2);
load('psg_isomap_demo_animals_intermediate_texture.mat')
imagesc(median(part_ratio(dim_list,1+[0 nn_list],:),3)',[1 7]);
set(gca,'YTickLabel',strvcat('std',num2str(nn_list(:))));
set(gca,'XTick',[1:length(dim_list)]);
set(gca,'XTickLabel',dim_list);
xlabel('dim');
ylabel('nn');
title('median participation ratio, intermed texture')
colorbar
%
subplot(2,3,4);
load('psg_isomap_demo_animals_intermediate_object.mat')
imagesc(median(part_ratio(dim_list,1+[0 nn_list],:),3)',[1 7]);
set(gca,'YTickLabel',strvcat('std',num2str(nn_list(:))));
set(gca,'XTick',[1:length(dim_list)]);
set(gca,'XTickLabel',dim_list);
xlabel('dim');
ylabel('nn');
title('median participation ratio, intermed obj')
colorbar
%
subplot(2,3,5);
load('psg_isomap_demo_animals_image.mat')
imagesc(median(part_ratio(dim_list,1+[0 nn_list],:),3)',[1 7]);
set(gca,'YTickLabel',strvcat('std',num2str(nn_list(:))));
set(gca,'XTick',[1:length(dim_list)]);
set(gca,'XTickLabel',dim_list);
xlabel('dim');
ylabel('nn');
title('median participation ratio, image')
colorbar
%
subplot(2,3,6);
load('psg_isomap_demo_animals_word.mat')
imagesc(median(part_ratio(dim_list,1+[0 nn_list],:),3)',[1 7]);
set(gca,'YTickLabel',strvcat('std',num2str(nn_list(:))));
set(gca,'XTick',[1:length(dim_list)]);
set(gca,'XTickLabel',dim_list);
xlabel('dim');
ylabel('nn');
title('median participation ratio, word')
colorbar
%
%heatmaps of Euclidean ratios
figure;
set(gcf,'Position',[100 100 1200 800]);
dim_list=[2:7];
%
subplot(2,3,1);
load('psg_isomap_demo_animals_texture.mat')
imagesc(median(pwr_ratio(dim_list,1+[0 nn_list],:),3)',[0.7 1]);
set(gca,'YTickLabel',strvcat('std',num2str(nn_list(:))));
set(gca,'XTick',[1:length(dim_list)]);
set(gca,'XTickLabel',dim_list);
xlabel('dim');
ylabel('nn');
title('median Euclidean pwr ratio, texture')
colorbar
%
subplot(2,3,2);
load('psg_isomap_demo_animals_intermediate_texture.mat')
imagesc(median(pwr_ratio(dim_list,1+[0 nn_list],:),3)',[0.7 1]);
set(gca,'YTickLabel',strvcat('std',num2str(nn_list(:))));
set(gca,'XTick',[1:length(dim_list)]);
set(gca,'XTickLabel',dim_list);
xlabel('dim');
ylabel('nn');
title('median Euclidean pwr ratio, intermed texture')
colorbar
%
subplot(2,3,4);
load('psg_isomap_demo_animals_intermediate_object.mat')
imagesc(median(pwr_ratio(dim_list,1+[0 nn_list],:),3)',[0.7 1]);
set(gca,'YTickLabel',strvcat('std',num2str(nn_list(:))));
set(gca,'XTick',[1:length(dim_list)]);
set(gca,'XTickLabel',dim_list);
xlabel('dim');
ylabel('nn');
title('median Euclidean pwr ratio, intermed obj')
colorbar
%
subplot(2,3,5);
load('psg_isomap_demo_animals_image.mat')
imagesc(median(pwr_ratio(dim_list,1+[0 nn_list],:),3)',[0.7 1]);
set(gca,'YTickLabel',strvcat('std',num2str(nn_list(:))));
set(gca,'XTick',[1:length(dim_list)]);
set(gca,'XTickLabel',dim_list);
xlabel('dim');
ylabel('nn');
title('median Euclidean pwr ratio, image')
colorbar
%
subplot(2,3,6);
load('psg_isomap_demo_animals_word.mat')
imagesc(median(pwr_ratio(dim_list,1+[0 nn_list],:),3)',[0.7 1]);
set(gca,'YTickLabel',strvcat('std',num2str(nn_list(:))));
set(gca,'XTick',[1:length(dim_list)]);
set(gca,'XTickLabel',dim_list);
xlabel('dim');
ylabel('nn');
title('median Euclidean pwr ratio, word')
colorbar


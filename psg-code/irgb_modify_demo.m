%irgb_modify_demo: demonstrate irgb_modify
%
%  See also: IRGB_MODIFY.
%
if ~exist('file_list')
    file_list={...
        'Cloth_Cloth_52.png',...
        'Cloth_Heavy_6.png',...
        'Cloth_Rough_11.png',...
        'Cloth_Stiff_2.png',...
        'Cloth_Synthetic_2.png',...
        'Cloth_Synthetic_5.png',...
        'Cloth_Velvet_11.png',...
        'Cloth_Waterproof_4.png',... 
        'Cloth_Wool_3.png'};
end
if ~exist('file_dir')
    file_dir='./GieselZaidiImages/';
end
nfiles=length(file_list)
img_origs=cell(1,nfiles);
imgs=cell(1,nfiles);
stats=cell(1,nfiles);
if ~exist('opts')
    opts=struct;
end
opts.range=[0 65536];
opts.if_show=1;
for ifile=1:nfiles
    img_origs{ifile}=double(imread(cat(2,file_dir,file_list{ifile})));
    label=strrep(file_list{ifile},'.png','');
    [imgs{ifile},stats{ifile},opts_used]=irgb_modify(img_origs{ifile},setfield(opts,'label',label));
     set(gcf,'Position',[100 200 950 600]); %make it fit on a landscape page
    disp(sprintf('processing image %3.0f: %s',ifile,file_list{ifile}));
end
disp('use figs_to_ps_landscape2 to save')


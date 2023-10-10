% irgb_gzranking_analyze: analyzes the Giesel-Zaidi ranking data
%
%   See also:  IRGB_GZRANKING_READ, IRGB_GZRANKING_GETIMAGES, BTC_DEFINE, IRGB_BTCSTATS, IRGB_MODIFY, IRGB_MODIFY_DICT.
%
if ~exist('filename_rank_def') filename_rank_def='C:/Users/jdvicto/Documents/jv/ENCL/Zaidi/RankingData.mat'; end
if ~exist('filename_img_def') filename_img_def='C:/Users/jdvicto/Documents/jv/ENCL/Zaidi/ImageNames.mat'; end
if ~exist('imagedata_fn_def') imagedata_fn_def='./irgb_gzranking_imagedata.mat';end
if ~exist('workspace_fn_def') workspace_fn_def='./irgb_gzranking_analyze_09Oct23.mat';end
if ~getinp('create a workspace from scratch','d',[0 1],0)
    workspace_fn=getinp('workspace file name','s',[],workspace_fn_def);
    load(workspace_fn);
else
    %
    filename_rank=getinp('Giesel-Zaidi ranking data file','s',[],filename_rank_def);
    filename_img=getinp('Giesel-Zaidi image name file','s',[],filename_img_def);
    imagedata_fn=getinp('Giesel-Zaidi image data file','s',[],imagedata_fn_def);
    %
    s=irgb_gzranking_read(filename_rank,filename_img); %read ranking data
    load(imagedata_fn); %read image data
    n_allimgs=length(image_data.names);
    disp(sprintf('image data loaded, %4.0f images',n_allimgs));
    %
    %for image stat calculation
    %
    btc_dict=btc_define;
    btc_n=length(btc_dict.codel); %number of coords, typically 10
    if ~exist('downsamples_def') downsamples_def=[1 2 4 8 16 32]; end
    if ~exist('opts_detrend') opts_detrend=[]; end
    if ~exist('opts_mlis') opts_mlis=[]; end
    %
    %for image manipulation
    %
    if ~exist('opts_modify') opts_modify=struct;end
    opts_modify=filldefault(opts_modify,'nrandpase',4);
    opts_modify=filldefault(opts_modify,'range',[0 65536]); %for png
    opts_modify.if_show=0;
    manip_dict=irgb_modify_dict;
    manip_list=fieldnames(manip_dict);
    if ~exist('manips') manips={'orig_bw','bw_whiten','bw_randph'}; end %manipulations to do
    %
    if ~exist(opts_detrend) opts_detrend=struct; end
    %
    downsamples=getinp('downsamplings for calculating image statistics','d',[1 64],downsamples_def);
    %
    if_frozen_psg=getinp('1 for frozen random numbers, 0 for new random numbers each time for session configuration, <0 for a specific seed','d',[-10000 1]);
    %
    if (if_frozen_psg~=0)
        rng('default');
        if (if_frozen_psg<0)
            rand(1,abs(if_frozen_psg));
        end
    else
        rng('shuffle');
    end
    %
    imgstats=struct;
    hw=waitbar(0,sprintf('calculating image stats for %3.0f images',n_allimgs));
    for k=1:n_allimgs
        [imgs,stats,opts_modify_used]=irgb_modify(image_data.rgbvals(:,:,:,k),opts_modify);
        for im=1:length(manips)
            manip=manips{im};
            modify_name=manip_dict.(manip).modify_name;
            if k==1
                imgstats.(manip)=zeros(length(downsamples),btc_n,n_allimgs);
            end
            img_modified_stack=imgs.(modify_name);
            nstack=size(img_modified_stack,3);
            btcstats_stack=zeros(length(downsamples),btc_n,nstack);
            for istack=1:nstack
                [btcstats_stack(:,:,istack),opts_mlis_used,opts_detrend_used]=...
                    irgb_btcstats(img_modified_stack(:,:,istack),downsamples,opts_mlis,opts_detrend);
            end
            imgstats.(manip)(:,:,k)=mean(btcstats_stack,3);
        end %next manipulation
        waitbar(k/n_allimgs,hw);
    end
    close(hw);
    %
    workspace_fn=getinp('workspace file name to save (e.g., ./irgb_gzranking_analyze_*.mat','s',[]);
    clear imgs image_data img_modified_stack
    save (workspace_fn);
end
%
%show ratings
%
figure;
set(gcf,'Position',[200 200 1200 800]);
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
%
disp('imgstats');
disp(imgstats);



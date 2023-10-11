% irgb_gzranking_analyze: analyzes the Giesel-Zaidi ranking data
%
% to do: correlations with pca of image stats
%
%   See also:  IRGB_GZRANKING_READ, IRGB_GZRANKING_GETIMAGES, BTC_DEFINE, IRGB_BTCSTATS, IRGB_MODIFY, IRGB_MODIFY_DICT.
%
if ~exist('filename_rank_def') filename_rank_def='C:/Users/jdvicto/Documents/jv/ENCL/Zaidi/RankingData.mat'; end
if ~exist('filename_img_def') filename_img_def='C:/Users/jdvicto/Documents/jv/ENCL/Zaidi/ImageNames.mat'; end
if ~exist('imagedata_fn_def') imagedata_fn_def='./irgb_gzranking_imagedata.mat';end
if ~exist('workspace_fn_def') workspace_fn_def='./irgb_gzranking_analyze_09Oct23.mat';end
if ~getinp('1 to create a workspace from scratch','d',[0 1],0)
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
%
%create an "average subject" 
%(can improve on this with subject-specific weightings, etc.)
%
s.subjs{end+1}='avg';
s.rankings(:,end+1,:)=mean(s.rankings,2);
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
disp('s');
disp(s);
disp('imgstats');
disp(imgstats);
%    
%show correlations of image stats with ratings
%
n_manips=length(manips);
n_ranked=length(s.ranked_list);
n_scales=length(downsamples);
corrs=zeros(n_props,btc_n,n_scales,n_subjs,n_manips);
corrs_pvals=zeros(n_props,btc_n,n_scales,n_subjs,n_manips);
for imanip=1:n_manips
    disp(' ')
    disp(sprintf('correlations between image statistics of images (%s) and ratings',manips{imanip}));
    imgstats_allscales=permute(imgstats.(manips{imanip})(:,:,s.ranked_list),[3 2 1]); %dims are image, btc_stat, scale)
    for isubj=1:n_subjs
        ratings=reshape(s.rankings(:,isubj,:),[n_ranked,n_props]);
        for iscale=1:n_scales
            disp(sprintf('subject %5s, scale %5.0f',s.subjs{isubj},downsamples(iscale)));
            [corrs(:,:,iscale,isubj,imanip),corrs_pvals(:,:,iscale,isubj,imanip)]=corr(ratings,imgstats_allscales(:,:,iscale));
            disp('corr');
            disp(corrs(:,:,iscale,isubj,imanip));
            disp('pval');
            disp(corrs_pvals(:,:,iscale,isubj,imanip));
        end
    end %isubj
end %imanip
%
%plot: each page a subject
%on each page: row is property, col is manip, and within that: scale x btc
%
corr_range=max(abs(corrs(:)));
if ~exist('p_list') p_list=[0.05 0.01,.001]; end
if ~exist('p_symb') p_symb={'.','x','*'}; end
%  
for isubj=1:n_subjs
    figure;
    set(gcf,'Position',[100 100 1100 800]);
    tstring=cat(2,'correlations of ratings and img stats: subj:',s.subjs{isubj});
    set(gcf,'Name',tstring);
    set(gcf,'NumberTitle','off');
    for imanip=1:n_manips
        for iprop=1:n_props
            subplot(n_props,n_manips,imanip+(iprop-1)*n_manips);
            imagesc(reshape(corrs(iprop,:,:,isubj,imanip),[btc_n,n_scales])',[-1 1]*corr_range);
            title(cat(2,s.props{iprop},' ',manips{imanip}),'Interpreter','none');
            set(gca,'YTick',[1:n_scales]);
            set(gca,'YTickLabel',downsamples);
            set(gca,'XTick',[1:btc_n]);
            set(gca,'XTickLabel',btc_dict.codel');
            %show p-values, fdr-corrected in each
            hold on;
            pvals=reshape(corrs_pvals(iprop,:,:,isubj,imanip),[btc_n,n_scales])';
            for ipv=1:length(p_list)
                pcrit_fdr=fdr(pvals(:),p_list(ipv));
                for iscale=1:n_scales
                    for ibtc=1:btc_n
                        if pvals(iscale,ibtc)<pcrit_fdr
                            plot(ibtc,iscale,cat(2,'k',p_symb{ipv}));
                        end
                    end
                end
            end
            %
            colorbar;
        end
    end
    axes('Position',[0.01,0.01,0.01,0.01]); %for text
    text(0,0,tstring,'Interpreter','none','FontSize',10);
    axis off;
end

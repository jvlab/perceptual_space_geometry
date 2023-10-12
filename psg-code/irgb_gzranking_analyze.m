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
%show rankings
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
    ylabel('ranked image index');
    xlabel('subject ID');
end
disp('s');
disp(s);
disp('imgstats');
disp(imgstats);
%
n_manips=length(manips);
n_ranked=length(s.ranked_list);
n_scales=length(downsamples);
n_pcs=getinp('number of principal components','d',[2 10],5);
if_showcorrs=getinp('1 to show table of correlations','d',[0 1],0);
if ~exist('p_list') p_list=[0.05 0.01,.001]; end %for highlighting on plots
if ~exist('p_symb') p_symb={'.','x','*'}; end
%    
%compute correlations of image stats with rankings
%
corrs=zeros(n_props,btc_n,n_scales,n_subjs,n_manips);
corrs_pvals=zeros(n_props,btc_n,n_scales,n_subjs,n_manips);
for imanip=1:n_manips
    if (if_showcorrs)
        disp(' ')
        disp(sprintf('correlations between image statistics of images (%s) and rankings',manips{imanip}));
    end
    imgstats_allscales=permute(imgstats.(manips{imanip})(:,:,s.ranked_list),[3 2 1]); %dims are image, btc_stat, scale)
    for isubj=1:n_subjs
        rankings=reshape(s.rankings(:,isubj,:),[n_ranked,n_props]);
        for iscale=1:n_scales
            [corrs(:,:,iscale,isubj,imanip),corrs_pvals(:,:,iscale,isubj,imanip)]=corr(rankings,imgstats_allscales(:,:,iscale));
            if (if_showcorrs)
                disp(sprintf('subject %5s, scale %5.0f',s.subjs{isubj},downsamples(iscale)));
                disp('corr');
                disp(corrs(:,:,iscale,isubj,imanip));
                disp('pval');
                disp(corrs_pvals(:,:,iscale,isubj,imanip));
            end
        end
    end %isubj
end %imanip
%
%plot correlations of rankings and img stats: each page a subject
%on each page: row is property, col is manip, and within that: scale x btc
%
corr_range=max(abs(corrs(:)));
%  
for isubj=1:n_subjs
    figure;
    set(gcf,'Position',[100 100 1100 800]);
    tstring=cat(2,'correlations of rankings and img stats: subj:',s.subjs{isubj});
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
%    
%compute correlations of image stats with princpal components of rankings
%
disp(sprintf('computing pcs of image statistics'));
%
imgstats_pcs=struct;
for imanip=1:n_manips
    % x = u*diag(sdiag).v' for [u,s,v]=svd(x)
    imgstats_pcs.(manips{imanip}).u=zeros(n_allimgs,n_pcs,n_scales);
    imgstats_pcs.(manips{imanip}).sd=zeros(n_pcs,n_scales);
    imgstats_pcs.(manips{imanip}).v=zeros(btc_n,n_pcs,n_scales);
    for iscale=1:n_scales
        x=permute(imgstats.(manips{imanip})(iscale,:,:),[3 2 1]);
        [u,sd,v]=svd(x);
        vred=v(:,1:n_pcs);
        if iscale>1 %align the pc's across scales by finding eiv of previous scale with largest dot-product
            vdots=vred'*vred_prev;
            for ipc=1:n_pcs
                vd=vdots(ipc,:);
                vmatch=min(find(abs(vd)==max(abs(vd))));
                if vd(vmatch)<0
                    vred(:,vmatch)=-vred(:,vmatch);
                    u(:,vmatch)=-u(:,vmatch);
                end
            end %ipc
        end
        imgstats_pcs.(manips{imanip}).u(:,:,iscale)=u(:,1:n_pcs);
        imgstats_pcs.(manips{imanip}).sd(:,iscale)=diag(sd(1:n_pcs,1:n_pcs));
        imgstats_pcs.(manips{imanip}).v(:,:,iscale)=vred;
        vred_prev=vred;
    end %iscale
end %imanip
%
%plot principal components of rankings
%
for imanip=1:n_manips
    figure;
    set(gcf,'Position',[100 100 1100 800]);
    tstring=cat(2,'pcs of image stats: ',manips{imanip});
    set(gcf,'Name',tstring);
    set(gcf,'NumberTitle','off');
    ncols=max(3,n_scales);
    nrows=3;
    pc_maxval=max(abs(imgstats_pcs.(manips{imanip}).v(:)));
    for iscale=1:n_scales
        subplot(nrows,ncols,iscale);
        imagesc(imgstats_pcs.(manips{imanip}).v(:,:,iscale),[-1 1]*pc_maxval);
        set(gca,'XTick',[1:n_pcs]);
        set(gca,'XTickLabel',[1:n_pcs]);
        set(gca,'YTick',[1:btc_n]);
        set(gca,'YTickLabel',btc_dict.codel');
        xlabel('pc');
        title(sprintf('scale %2.0f',downsamples(iscale)));
        %
        subplot(nrows,ncols,iscale+ncols);
        plot([1:n_pcs],imgstats_pcs.(manips{imanip}).sd(:,iscale),'k.-');
        set(gca,'YLim',[0 max(abs(imgstats_pcs.(manips{imanip}).sd(:)))]);
        set(gca,'XLim',[0 n_pcs+1]);
        set(gca,'XTick',[1:n_pcs]);
        xlabel('pc');
        ylabel('wt');
    end %iscale
    subplot(nrows,ncols,1+2*ncols);
    hcb=colorbar;
    set(hcb,'Ticks',[0 .5 1])
    set(hcb,'TickLabels',[-1 0 1]*pc_maxval)
    axis off;
    %
    axes('Position',[0.01,0.01,0.01,0.01]); %for text
    text(0,0,tstring,'Interpreter','none','FontSize',10);
    axis off;
end %imanip
corrs_pcs=zeros(n_props,n_pcs,n_scales,n_subjs,n_manips);
corrs_pcs_pvals=zeros(n_props,n_pcs,n_scales,n_subjs,n_manips);
%
%compute correlations between pcs and rankings
%the pcs are in imgstats_pcs.(manips{imanip}).u, dim 1: image, dim2: pc, dim 3: scale
%
for imanip=1:n_manips
    pcs_allscales=imgstats_pcs.(manips{imanip}).u(s.ranked_list,:,:); %dims are image, pc, scale)
    for isubj=1:n_subjs
        rankings=reshape(s.rankings(:,isubj,:),[n_ranked,n_props]);
        for iscale=1:n_scales
            [corrs_pcs(:,:,iscale,isubj,imanip),corrs_pcs_pvals(:,:,iscale,isubj,imanip)]=corr(rankings,pcs_allscales(:,:,iscale));
        end
    end %isubj
end %imanip
%
%plot: each page a subject
%on each page: row is property, col is manip, and within that: scale x pc
%
corr_pcs_range=max(abs(corrs_pcs(:)));
%  
for isubj=1:n_subjs
    figure;
    set(gcf,'Position',[100 100 1100 800]);
    tstring=cat(2,'correlations of rankings and pcs of img stats: subj:',s.subjs{isubj});
    set(gcf,'Name',tstring);
    set(gcf,'NumberTitle','off');
    for imanip=1:n_manips
        for iprop=1:n_props
            subplot(n_props,n_manips,imanip+(iprop-1)*n_manips);
            imagesc(reshape(corrs_pcs(iprop,:,:,isubj,imanip),[n_pcs,n_scales])',[-1 1]*corr_pcs_range);
            title(cat(2,s.props{iprop},' ',manips{imanip}),'Interpreter','none');
            set(gca,'YTick',[1:n_scales]);
            set(gca,'YTickLabel',downsamples);
            set(gca,'XTick',[1:n_pcs]);
            set(gca,'XTickLabel',[1:n_pcs]);
            xlabel('pc');
            %show p-values, fdr-corrected in each
            hold on;
            pvals=reshape(corrs_pcs_pvals(iprop,:,:,isubj,imanip),[n_pcs,n_scales])';
            for ipv=1:length(p_list)
                pcrit_fdr=fdr(pvals(:),p_list(ipv));
                for iscale=1:n_scales
                    for ipc=1:n_pcs
                        if pvals(iscale,ipc)<pcrit_fdr
                            plot(ipc,iscale,cat(2,'k',p_symb{ipv}));
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
%
% compute and display cross-subject correlations
%
disp(' ');
disp('subjs:');
disp(s.subjs);
corrs_intersubj_rankings=zeros(n_subjs,n_subjs,n_props);
corrs_intersubj_rankings_pvals=zeros(n_subjs,n_subjs,n_props);
for iprop=1:n_props
    disp(sprintf('intersubject correlations for %s              p-values',s.props{iprop}));
    [corrs_intersubj_rankings(:,:,iprop),corrs_intersubj_rankings_pvals(:,:,iprop)]=corr(s.rankings(:,:,iprop));
    disp([corrs_intersubj_rankings(:,:,iprop) corrs_intersubj_rankings_pvals(:,:,iprop)]);
end %iprop
%
% compute and display cross-property correlations
%
disp(' ');
disp('props:');
disp(s.props);
corrs_interprop_rankings=zeros(n_props,n_props,n_subjs);
corrs_interprop_rankings_pvals=zeros(n_props,n_props,n_subjs);
for isubj=1:n_subjs
    disp(sprintf('interproperty correlations for %s             p-values',s.subjs{isubj}));
    [corrs_interprop_rankings(:,:,isubj),corrs_interprop_rankings_pvals(:,:,iprop)]=corr(permute(s.rankings(:,isubj,:),[1 3 2]));
    disp([corrs_interprop_rankings(:,:,isubj) corrs_interprop_rankings_pvals(:,:,isubj)]);
end %iprop

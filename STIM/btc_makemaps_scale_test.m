%btc_makemaps_scale_test: test btc_makemaps_scale
%
if ~exist('specs')
    specs{1}.g=.4;
    specs{2}.b=.6;
    specs{3}.c=.6;
    specs{4}.d=.6;
    specs{5}.e=.6;
    specs{6}.t=.9;
    specs{7}.u=.9;
    specs{8}.v=.9;
    specs{9}.a=.8;
    specs{10}.b=.5;
    specs{10}.c=-.5;
end
if ~exist('npxls') npxls=[32 48]; end
if ~exist('nmaps') nmaps=3; end
if ~exist('block_stagger') block_stagger=[1 1;1 2; 2 1; 2 2;1 4]; end
nspecs=length(specs);
nblock_stagger=length(block_stagger);
%
dict=btc_define;
%
opts=[];
opts.nmaps=nmaps;
opts.show=[];
opts.area=npxls;
%
optsused=cell(nspecs,nblock_stagger);
errs=cell(nspecs,nblock_stagger);
imgs=cell(nspecs,nblock_stagger);
metro=cell(nspecs,nblock_stagger);
methods=cell(nspecs,1);
%
for ispec=1:nspecs
    tstring=[];
    lets=fieldnames(specs{ispec});
    for ibtc=1:length(lets)
        let=lets{ibtc};
        tstring=cat(2,tstring,sprintf('%s=%5.3f ',let,specs{ispec}.(let)));
    end
    tstring=deblank(tstring);
    augcoords=btc_augcoords(specs{ispec},dict);
    methods{ispec}=augcoords.method{1};
    figure;
    set(gcf,'Position',[50 100 1200 750]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name','btc_makemaps_scale_test');
    set(gcf,'Name',tstring);
    for iblock_stagger=1:nblock_stagger
        opts_scale=opts;
        opts_scale.blocked=block_stagger(iblock_stagger,1);
        opts_scale.staggered=block_stagger(iblock_stagger,2);
        [imgs{ispec,iblock_stagger},optsused{ispec,iblock_stagger},errs{ispec,iblock_stagger},metro{ispec,iblock_stagger}]=...
            btc_makemaps_scale(methods{ispec},opts_scale,dict);
        ng=size(methods{ispec}.p2x2,1); %number of gray levels
        for imap=1:nmaps
            subplot(nmaps,nblock_stagger,iblock_stagger+(imap-1)*nblock_stagger);
            im_max=(ng-1)*opts_scale.staggered.^2;
            imagesc(imgs{ispec,iblock_stagger}(:,:,imap),[0 im_max]);
            axis equal;
            axis tight;
            colormap gray;
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            if (imap==1)
                title(sprintf('block %2.0f stagger %2.0f',block_stagger(iblock_stagger,:)));
            end
        end
    end
    axes('Position',[0.02,0.08,0.01,0.01]); %for text
    text(0,0,tstring,'Interpreter','none');
axis off;

end
    

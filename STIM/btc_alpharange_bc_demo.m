%btc_alpharange_bc_demo:b demonstrate btc_alpharange, as a function of second-order params b and c
% creates heatmaps of range of available alpha, and also sample maps
%
%   See also:  BTC_DEFINE, BTC_ALPHARANGE, BTC_AUGCOORDS, GENMRFM,
%   BTC_ALPHARANGE_BG_DEMO.
%
dict=btc_define;
if ~exist('mapsize') mapsize=64;end
nsamps=getinp('number of samples in [0 1]','d',[2 1000],50);
n=2*nsamps+1;
cvals=[-nsamps:nsamps]/nsamps;
bvals=[-nsamps:nsamps]/nsamps;
%
alpha_vals=NaN(n,n,3); %dim3 for maxent, min and max
labels={'maxent','min','max','range'};
%
aug_opts=struct;
aug_opts.ifstd=1;
aug_opts.nocheck=1;
for ic=1:n
    for ib=1:n
        c=cvals(ic);
        b=bvals(ib);
        if 1==1 %always feasible
            spec=struct;
            spec.b=b;
            spec.c=c;
            augcoords=btc_augcoords(spec,dict,aug_opts);
            method=augcoords.method{1};
            alpha_vals(ib,ic,1)=method.corrs.alpha;
            p2x2=method.p2x2;
            [amin,amax,p2x2_extremes]=btc_alpharange(p2x2);
            alpha_vals(ib,ic,2)=amin;
            alpha_vals(ib,ic,3)=amax;
        end
    end
end
figure;
set(gcf,'Position',[100 100 1000 800]);
for isub=1:4
    subplot(2,2,isub)
    if isub==4
        imagesc(cvals,bvals,alpha_vals(:,:,3)-alpha_vals(:,:,2),[0 2]);
    else
        imagesc(cvals,bvals,alpha_vals(:,:,isub),[-1 1]);
    end
    set(gca,'YDir','normal');
    xlabel('c');
    ylabel('b');
    axis tight;
    axis equal;
    colorbar;
    title(labels{isub});
end
%
%sample maps
%
while(1==1)
    bc=getinp('values of b and c','f',[-1 1],[0 0]);
    spec.b=bc(1);
    spec.c=bc(2);
    augcoords=btc_augcoords(spec,dict,aug_opts);
    method=augcoords.method{1};
    alpha_maxent=method.corrs.alpha;
    p2x2=method.p2x2;
    [alpha_min,alpha_max,p2x2_extremes]=btc_alpharange(p2x2);
    disp(sprintf(' alpha values (maxent, min, max): %7.4f %7.4f %7.4f',alpha_maxent,alpha_min,alpha_max))
    %
    opts.show=1;
    opts.pblocks=p2x2_extremes(:,:,:,:,1);
    img_min=genmrfm(opts,mapsize);
    title(sprintf('b %7.3f, c %7.3f, min alpha %7.3f',bc,alpha_min));
    opts.pblocks=p2x2_extremes(:,:,:,:,2);
    img_max=genmrfm(opts,mapsize);
    title(sprintf('b %7.3f, c %7.3f, max alpha %7.3f',bc,alpha_max));
end

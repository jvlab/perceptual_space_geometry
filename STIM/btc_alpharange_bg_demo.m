%btc_alpharange_bg_demo:b demonstrate btc_alpharange, as a function of first- and second-order params b and g
% creates heatmaps of range of available alpha, and also sample maps
%
%   See also:  BTC_DEFINE, BTC_ALPHARANGE, BTC_AUGCOORDS, GENMRFM.
%
dict=btc_define;
if ~exist('mapsize') mapsize=64;end
nsamps=getinp('number of samples in [0 1]','d',[2 1000],50);
n=2*nsamps+1;
gvals=[-nsamps:nsamps]/nsamps;
bvals=[-nsamps:nsamps]/nsamps;
%
alpha_vals=NaN(n,n,3); %dim3 for maxent, min and max
labels={'maxent','min','max'};
%
aug_opts=struct;
aug_opts.ifstd=1;
aut_opts.nocheck=1;
for ig=1:n
    for ib=1:n
        g=gvals(ig);
        b=bvals(ib);
        if b>=2*abs(g)-1 %(b,g) feasible
            spec=struct;
            spec.b=b;
            spec.g=g;
            augcoords=btc_augcoords(spec,dict,aug_opts);
            method=augcoords.method{1};
            alpha_vals(ib,ig,1)=method.corrs.alpha;
            p2x2=method.p2x2;
            [amin,amax,p2x2_extremes]=btc_alpharange(p2x2);
            alpha_vals(ib,ig,2)=amin;
            alpha_vals(ib,ig,3)=amax;
        end
    end
end
figure;
set(gcf,'Position',[100 100 1000 800]);
for isub=1:3
    subplot(2,2,isub)
    imagesc(gvals,bvals,alpha_vals(:,:,isub),[-1 1]);
    set(gca,'YDir','normal');
    xlabel('g');
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
    bg=getinp('values of b and g','f',[-1 1],[0 0]);
    spec.b=bg(1);
    spec.g=bg(2);
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
    title(sprintf('b %7.3f, g %7.3f, min alpha %7.3f',bg,alpha_min));
    opts.pblocks=p2x2_extremes(:,:,:,:,2);
    img_max=genmrfm(opts,mapsize);
    title(sprintf('b %7.3f, g %7.3f, max alpha %7.3f',bg,alpha_max));
end

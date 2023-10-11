function  [statvec,opts_mlis_used,opts_detrend_used]=mlis_btcstats(img,downsamples,opts_mlis,opts_detrend)
% [statvec,opts_mlis_used,opts_detrend_used]=mlis_btcstats(img,downsamples,opts_mlis,opts_detrend)
% computes binary texture statistics from a rectangular image, after detrending, downsampling, and binarizing
% 
% The starting phase for downsampling is stepped through the grid of possible values [0:d-1,0:d-1]
%   for a downsampling at scale d.
% Does not whiten.
%
% img: a rectangular image, black=0, white=1
% downsamples: a vector of integers, typically in [1 2 4 8 16 32]
% opts_mlis
%    mils_opts.binarize is 'median','mean', or a numeric value for quantile; defaults to 'median'
% opts_detrend
%    opts_detrend.method is 'constant','linear','bilinear','quadratic', defaults to 'quadratic'
%
% statvec: array of size [length(downsamples), nbtc]
%   Each row is the binary texture stats for each level of downsampling
%   (NaN's if any dimension of img is <2*downsample)
%
%  See slao:  BTCSTATS_ROI_DETREND, BTCSTATS_ROI_TIFF_RUN, IRGB_PSG_SESS_SETUP, IRGB_PSG_IMGS_SETUP,
%      BTC_DEFINE, MLIS_BTCSTATS.
if (nargin<=2)
    opts_mlis=struct;
end
if (nargin<=3)
    opts_detrend=struct;
end
opts_mlis=filldefault(opts_mlis,'binarize','median');
opts_detrend=filldefault(opts_detrend,'method','quadratic');
opts_mlis_used=opts_mlis;
opts_detrend_used=opts_detrend;
btc_dict=btc_define;
btc_n=length(btc_dict.codel); %number of coords, typically 10
statvec=zeros(length(downsamples),btc_n);
%
switch opts_detrend.method
    case 'constant'
        terms=[0 0];
    case 'linear'
        terms=[0 0;1 0;0 1];
    case 'bilinear'
        terms=[0 0;1 0;0 1;1 1];
    case 'quadratic'
        terms=[0 0;1 0;0 1;2 0;1 1;0 2];
    otherwise
        warning(sprintf('Detrend method %s not recognized.  No detrending.',opts_detrend.method));
        terms=[];       
end
if ~isempty(terms)
    img_detrended=btcstats_roi_detrend(img,ones(size(img)),terms);
else
    img_detrended=img;
end
%
for idown=1:length(downsamples)
    downsample=downsamples(idown);
    maxoffs=min(size(img)-2*downsample,downsample-1);
    if all(maxoffs>=0)
        opts_mlis.btc_down=downsample;
        %weighted sum of statistics -- since each phase shifts have different numbers of pixels after downsampling
        statvec_accum=zeros(1,btc_n);
        npxls_accum=0;
        for xoff=0:maxoffs(1)
            xmax=xoff+downsample*floor((size(img,1)-xoff)/downsample);
            xrange=[xoff+1:xmax];
            for yoff=0:maxoffs(2)
                ymax=yoff+downsample*floor((size(img,2)-yoff)/downsample);
                yrange=[yoff+1:ymax];
                npxls=length(xrange)*length(yrange)/(downsample.^2);
%                [xoff yoff length(xrange) length(yrange) npxls]
                statvec_accum=statvec_accum+npxls*mlis_btcstats(img_detrended(xrange,yrange),opts_mlis,btc_dict);
                npxls_accum=npxls_accum+npxls;
            end %yoff
        end %xoff
        statvec(idown,:)=statvec_accum/npxls_accum;
    else %downsample too big
        statvec(idown,:)=NaN;
    end
end %downsample
return

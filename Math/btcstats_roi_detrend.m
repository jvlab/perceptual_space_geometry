function [img_detrended,coefs,regressors]=btcstats_roi_detrend(img,roi,terms)
% [img_detrended,coefs,regressors]=btcstats_roi_detrend(img,roi,terms) does a polynomial detrending of an image
% within a mask
%
% img: 2-d array of intensity values
% roi: 2-d binary array indicating region to detrend, 1 to include
% terms: [nterms 2] array of polynomial terms to include
%   rows are as follows: [0 0]: constant term; [1 0] and [0 1], linear terms for first and second dimension, etc.
%   non-constant terms are computed w.r.t. mean values of x and y coords of roi
%
% img_detrended: image with the trends removed, set to NaN where roi=0
% coefs: column of size nterms, the regression coefficients
% regrssors: regressors, size=[size(img) nterms], each slice has NaN's outside of roi
%
% See also:  BTCSTATS_ROI_TIFF_RUN, REGRESS.
%
img_detrended=NaN(size(img));
img_detrended(roi==1)=img(roi==1); %put NaN's where roi=0
nterms=size(terms,1);
coefs=zeros(1,nterms);
regressors=NaN([prod(size(img)) nterms]);
%
nroi=sum(roi(:)==1);
%
xycoords=cell(2,1);
xycoords_roi=zeros(nroi,2);
xycoords{1}=repmat(1:size(img,1),[1 size(img,2)]);
xycoords{2}=repmat(1:size(img,2)',[size(img,1) 1]);
for ixy=1:2
    xycoords_roi(:,ixy)=xycoords{ixy}(roi==1);
end
xycoords_roi=xycoords_roi-mean(xycoords_roi,1);
%create regressors
regressor_columns=ones(nroi,nterms);
for iterm=1:nterms
    for ixy=1:2
        if terms(iterm,ixy)>0
            regressor_columns(:,iterm)=regressor_columns(:,iterm).*xycoords_roi(:,ixy).^terms(iterm,ixy);
        end %pwr>0?
    end %ixy
    regressors(roi==1,iterm)=regressor_columns(:,iterm);
end %iterm
regressors=reshape(regressors,[size(img),nterms]);
coefs=regress(img(roi==1),regressor_columns);
img_detrended(roi==1)=img(roi==1)-regressor_columns*coefs;
return

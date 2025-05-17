function [edirs_pred,allmag,allmag_pred,allmag_std]=btc_soid_predict(edirs,edirs_thr,qform_fitted,dict,verbose)
% [edirs_pred,allmag,allmag_pred,allmag_std]=btc_soid_predict(edirs,edirs_thr,qform_fitted,dict,verbose) calculates
% predited in-plane threshold magnitudes from a quadratic form
%
% edirs: an edirs structure (see btc_soid_test) for which edirs.xy.uvecs are the test directions
% edirs_thr: an edirs structure containing, for each plane xy, a structure thresh_mags of magnitudes
% qform_fitted: the quadratic form, length(dict.codel) x length(dict.codel) (i.e,. 10 x 10)
% dict: the standard dictionary
% verbose=1 to log
%
% edirs_pred:  an edirs structure, with the predicted threshold magnitudes
% allmag: a column vector of all measured threshold magnitudes, stripped from edirs_thr
% allmag_pred: a column vector of all predicted threshold magnitudes, stripped from edirs_pred
% allmag_std: a column vector of measured theshold standard devs,
%    determined from error bars of edirs_thr.thresh_mags_eblo, edirs_thr.thresh_mags.ebhi
%
%   April 2014:  added allmag_std as output, if error bars are available in input
%
%   See also:  BTC_SOID_FIT, BTC_SOID_FIND, BTC_AUGCOORDS, BTC_DEFINE.
%
%extract each plane
planes=char(fieldnames(edirs));
nplanes=size(planes,1);
edirs_pred=[];
allmag=[];
allmag_pred=[];
allmag_std=[];
cl=0.95; % error bar confidence limit
eb_to_std=1/norminv(1-(1-cl)/2);
if (nargin<=4) verbose=0; end
for iplane=1:nplanes
    if (verbose==1)
        disp(sprintf('determining predicted thresholds in plane %s',planes(iplane,:)));
    end
    edir=getfield(edirs,planes(iplane,:));
    edir_thr=getfield(edirs_thr,planes(iplane,:));
    allmag=[allmag;edir_thr.thresh_mags];
    if isfield(edir_thr,'thresh_mags_eblo') & isfield(edir_thr,'thresh_mags_ebhi')
        allmag_std=[allmag_std;eb_to_std*(edir_thr.thresh_mags_ebhi-edir_thr.thresh_mags_eblo)/2];
    end
    edir_pred=[];
    edir_pred.ndirs=edir.ndirs;
    %for each direction, find the threshold
    vecs_inplane=edir.uvecs;
    spec=[];
    for icond=1:edir.ndirs
       for ix=1:2
            %the "3-ix" is so that the SECOND coordinate of vecs_inplane 
            %is the FIRST btc coordinate (in planes)
            spec=setfield(spec,planes(iplane,ix),vecs_inplane(icond,3-ix));
       end
       [specx,avec,results,ou]=btc_soid_find(spec,1,qform_fitted,dict);
       thr_inplane=[getfield(specx,planes(iplane,2)),getfield(specx,planes(iplane,1))];
       edir_pred.thresh_mags(icond,1)=sqrt(sum(thr_inplane.^2));
   end
   edirs_pred=setfield(edirs_pred,planes(iplane,:),edir_pred);
   allmag_pred=[allmag_pred;edir_pred.thresh_mags];
end
%use BTC_SOID_FIT to find the magnitude
% [qfit,results,ou_fit,ou_dict]=btc_soid_fit(edirs_setup,edirs_data,opts_fit,dict)
return

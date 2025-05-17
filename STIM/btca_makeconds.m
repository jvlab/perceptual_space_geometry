function [conds,opts_used]=btca_makeconds(parsed,opts)
% [conds,opts_used]=btca_makeconds(parsed,opts) creates an array with the BigConds format for augmented btc experiments
%
%  parsed.npred is the number of textures to predict (arbitrary directions)
%  parsed.nmeas is the number of axes to measure (axis directions)
%   these are from ptbxcvt_btclabs with aflag=1
%     for in-sample measurements: (axes in btc space)
%  opts.meas_mults: list of mutlipliers for axis directions (defaults to [0.2 0.4 0.6 0.8 1.0]
%  opts.meas_repts: number of repetitions for the measured directions, defaults to 1
%     for out-of-sample measurements: (e.g., natural textues)
%  opts.pred_mults: list of mutlipliers for axis directions (defaults to [1]
%  opts.pred_repts: number of repetitions for the predicted directions, defaults to 3
%  opts.ifplot=1 to plot, defaults to 0
%  opts.tf_range: target flag range: 0 for only background structured, 1 for only target structured, [0 1] for both (default)
%     (added 07Jul17)
%
%  modification 28 Jan 13 to allow for opts.meas_mults and opts.pred_mults
%  to have several rows; they are used cyclically
%
% BIG KLUDGE:  Columns 3 and 4 indicate which direction to use, first direction is along X-axis,
% then counterclockwise, occupying first pi radians ONLY
% so that negative values for multipliers can be used
%     the first parsed.nmeas directions are the measurement directions
%     the next parsed.npred directions are the predictiondirections
% note that this does not specify what the directions are; this is specified
%     in ptbxcvt_btcgetnex: need.texture_meas_dirs,need.texture_meas_maxmults
%
% conds has a row for each condition
%    Column 1: 1 or 0 to indicate whether this condition is to be included in current block, usually 1
%    Column 2: Overall multiplier for parameter strengths
%    Column 3: Value of parameter 1 (e.g., alpha or theta) [to be multiplied by Column 2]
%    Column 4: Value of parameter 2 (e.g., gamma) [to be multiplied by Column 2]
%    Column 5: 0 if texture is background (with coinflip texture target), 1 if roles are reversed
%    Column 6: Target position: 1->top, 2-> right, 3->bottom, 4->left
%    Column 7: subject's response
%    Column 8: index into table of trial descriptors (i.e., which row of Conds)
%    Columns 7 and 8 are ignored
%
%   See also:  PTBXCVT_MAKERAY, PTBXCVT_DEMO, BTC_BIGCONDS_SHOW, MAKEBIGCONDP_BTC, BTCA_DEFLAYOUT.
% 
if nargin<=1 opts=[]; end
layout=btca_deflayout;
opts=filldefault(opts,'meas_mults',layout.meas_mults);
opts=filldefault(opts,'meas_repts',layout.meas_repts);
opts=filldefault(opts,'pred_mults',layout.pred_mults);
opts=filldefault(opts,'pred_repts',layout.pred_repts);
opts=filldefault(opts,'ifplot',0);
opts=filldefault(opts,'tf_range',[0 1]);
opts_used=opts;
%
conds=[];
%
ndirs=parsed.nmeas+parsed.npred;
for idir=1:ndirs
    ang=(idir-1)*pi/ndirs;
    if (idir<=parsed.nmeas)
%       conds=[conds; ptbxcvt_makeray(cos(ang),sin(ang),opts.meas_mults,opts.meas_repts)];
%       next four lines allow for cyclic use of several rows of meas_mults
        irow=idir;
        nrows=size(opts.meas_mults,1);
        mults=opts.meas_mults(mod(irow-1,nrows)+1,:);
        conds=[conds; ptbxcvt_makeray(cos(ang),sin(ang),mults,opts.meas_repts,opts.tf_range)];
    else
%       conds=[conds; ptbxcvt_makeray(cos(ang),sin(ang),opts.pred_mults,opts.pred_repts)];
%       next four lines allow for cyclic use of several rows of meas_mults
        irow=idir-parsed.nmeas;
        nrows=size(opts.pred_mults,1);
        mults=opts.pred_mults(mod(irow-1,nrows)+1,:);
        conds=[conds; ptbxcvt_makeray(cos(ang),sin(ang),mults,opts.pred_repts,opts.tf_range)];
    end
end
%
if (opts.ifplot)
    ptbxcvt_plot(conds+repmat([0 0 0 0 0 0 1 1],size(conds,1),1),setfield([],'pool',1)); %plot it
    title('conds');
end
return

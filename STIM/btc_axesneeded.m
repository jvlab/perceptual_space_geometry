function [axes_sym,axes_full]=btc_axesneeded(codes,dict)
% [axes_sym,axes_full]=btc_axesneeded(codes,dict) determines which coordinate
% axes are needed to fully characterize the metric along the axes given in codes.
%
% codes: a string (size=[1 n]) of symbols drawn from dict.codel
%   if codes is empty, then it defaults to dict.codel
% dict: the dictionary of 2x2 coordinates, from btc_define (or, btc_define([] if not passed)
%
% axes_sym: array of size[: 1] of the necessary axes, assuming replication by symmetry
%     special case: 'te' is substituted for 'td' for the symmetric set
% axes_full: array of size [n 1] of the necessary axes, without symmetry
%
%   See also:  BTC_EXPTNAME, BTC_DEFINE, BTC_ROTCODE, BTC_HVI, BTC_PAIRSNEEDED.
%
if (nargin<=1)
    dict=btc_define([]);
end
if isempty(codes)
    codes=dict.codel;
end
%
axes_sym=[];
axes_full=[];
%
nin=size(codes,2);
if (nin<1) return; end
%
% this is overkill, bue we keep the same logic as btc_pairsneeded
allaxes=codes';
axes_rots=[];
for icomb=1:size(allaxes,1)
    pcand=allaxes(icomb);
    axes_full=[axes_full;pcand];
    axes_thisrot=[];
    pcand_hvi=btc_hvi(pcand,dict);
    for irot=0:3;
        axes_thisrot=[axes_thisrot;btc_rotcode(pcand,irot,dict)];
        axes_thisrot=[axes_thisrot;btc_rotcode(pcand_hvi,irot,dict)];
    end
    %is this combination a rotation of what is already included in axes_sym?
    if isempty(strmatch(pcand,axes_rots,'exact'));
        axes_sym=[axes_sym;pcand];
        axes_rots=[axes_rots;axes_thisrot];
    end
end
return

function [planes_sym,planes_full,planes_sym_unsub,planes_mapping]=btc_pairsneeded(codes,dict)
% [planes_sym,planes_full,planes_sym_unsub,planes_mapping]=btc_pairsneeded(codes,dict) determines which pairs
% of coordinates (i.e., which planes) are needed to fully characterize
% the metric along the axes given in codes.
%
% codes: a string (size=[1 n]) of symbols drawn from dict.codel
%   if codes is empty, then it defaults to dict.codel
% dict: the dictionary of 2x2 coordinates, from btc_define (or, btc_define([] if not passed)
%
% planes_sym: array of size[: 2] of the necessary planes, assuming replication by symmetry
%     special case: 'te' is substituted for 'ud' for the symmetric set
% planes_full: array of size[n*(n-1)/2 2] of the necessary planes, without symmetry
% planes_sym_unsub: array without special case subsitutions
% planes_mapping indicates mapping of reduced subset (planes_sym, planes_sym_unsub) to the full set
%   (above added Sept 2013)
%
%   See also:  BTC_EXPTNAME, BTC_DEFINE, BTC_ROTCODE, BTC_HVI, BTC_AXESNEEDED.
%
if (nargin<=1)
    dict=btc_define([]);
end
if isempty(codes)
    codes=dict.codel;
end
%
planes_sym=[];
planes_full=[];
planes_sym_unsub=[];
planes_mapping=[];
%
nin=size(codes,2);
if (nin<=1) return; end
%
allpairs=nchoosek([1:length(codes)],2);
pairs_rots=[];
for icomb=1:size(allpairs,1)
    pcand=btc_exptname(codes(allpairs(icomb,:)),dict);
    planes_full=[planes_full;pcand];
    pairs_thisrot=[];
    pcand_hvi=btc_hvi(pcand,dict);
    for irot=0:3;
        pairs_thisrot=[pairs_thisrot;btc_exptname(btc_rotcode(pcand,irot,dict))];
        pairs_thisrot=[pairs_thisrot;btc_exptname(btc_rotcode(pcand_hvi,irot,dict))];
    end
    %is this combination a rotation of what is already included in planes_sym?
    if isempty(strmatch(pcand,pairs_rots,'exact'));
        uplanes=unique(pairs_thisrot,'rows');
        planes_sym_unsub=[planes_sym_unsub;pcand];
        planes_mapping.sym_unsub.(pcand)=uplanes;
        if (pcand=='ud')
            pcand='te';
        end
        planes_sym=[planes_sym;pcand];
        planes_mapping.sym.(pcand)=uplanes;
        %
        pairs_rots=[pairs_rots;pairs_thisrot];
    end
end
return

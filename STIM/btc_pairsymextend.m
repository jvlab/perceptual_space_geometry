function [pairs,qones,qone]=btc_pairsymextend(pairname,codes,dict)
% [pairs,qones,qone]=btc_pairsymextend(pairname,codes,dict) determines the coordinate pairs
% that are equivalent to a given coordinate pair, by rotation and hv flip
%
% pairname: a pair of letter codes for btc coordinates (also works if pairname is a singleton)
% codes: a string (size=[1 n]) of symbols drawn from dict.codel
%   if codes is empty, then it defaults to dict.codel
% dict: the dictionary of 2x2 coordinates, from btc_define (or, btc_define([] if not passed)
%
% pairs: a strvcat of the pairs that are equivalent by rotation or hv flip
% qones: a symmetric matrix of size length(codes), containing 1's in locs corresponding to all elements of pairs 
% qone:  a symmetric matrix of size length(codes), containing 1's in locs corresponding to pair
%
% size(pairs,1)=multiplicity of the plane (under rotational and h-v interchange symmetries)
%
%   See also:  BTC_EXPTNAME, BTC_DEFINE, BTC_ROTCODE, BTC_HVI, BTC_PAIRSNEEDED.
%
if (nargin<=2)
    dict=btc_define([]);
end
if (nargin<=1)
    codes=[];
end
if isempty(codes)
    codes=dict.codel;
end
%
qones=zeros(length(codes));
qone=zeros(length(codes));
%
for ic=1:2
    nc(ic)=find(pairname(min(ic,length(pairname)))==codes);
end
qone(nc(1),nc(2))=1;
qone=double((qone+qone')>0);
%
pair_hvi=btc_exptname(btc_hvi(pairname,dict));
pairs=[pairname;pair_hvi];
for irot=1:3;
    pairs=[pairs;btc_exptname(btc_rotcode(pairname,irot,dict))];
    pairs=[pairs;btc_exptname(btc_rotcode(pair_hvi,irot,dict))];
end
pairs=unique(pairs,'rows');
%
%make sure that each element of pairs is in codes
%
pairs_cand=pairs;
pairs=[];
for icand=1:size(pairs_cand,1)
    if all(ismember(pairs_cand(icand,:),codes))
        pairs=strvcat(pairs,pairs_cand(icand,:));
    end
end
%put ones into qones
for isym=1:size(pairs,1)
    for ic=1:2
        nc(ic)=find(pairs(isym,min(ic,length(pairname)))==codes);
    end
    qones(nc(1),nc(2))=1;
end
qones=double((qones+qones')>0);
return

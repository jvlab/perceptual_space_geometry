function counts=btc_map2counts(img)
% counts=btc_map2counts(img) calculates all of the block probabilities in a
% [0 1] image
%
% img: an array of [0 1]
%
% counts ordering conforms to that of genmrfm, getcorrs_p2x2, etc.
% counts can be converted into probabilities by p2x2=counts/sum(counts(:))
%
%    See also:  GETCORRS_P2X2, BTC_CORRS2VEC, MAPUBI.
%
counts_ubi=mapubi(img,[2 2]);
%deal with possiblity that some blocks are not present
indices=counts_ubi(:,4)+2*counts_ubi(:,3)+4*counts_ubi(:,2)+8*counts_ubi(:,1);
counts_map=zeros(1,16);
counts_map(indices+1)=counts_ubi(:,5);
%deal with ordering of pixels in mapubi:  
counts=permute(reshape(counts_map,[2 2 2 2]),[4 2 3 1]);
return



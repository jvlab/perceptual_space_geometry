function syms=btc_whichsym(vec,ipn,tol,dict)
% syms=btc_whichsym(vec,ipn,tol,dict) determines which symmetry motions
% preserve the block probabilities
%
% vec: a 10-vector 
% ipn: 1 for the symmetries that leave a vec invariant, -1 for those that flip its sign, defaults to 1
% tol: a tolerance for matching, defaults to 0.0001
% dict: the dictionary of 2x2 coordinates, from btc_define (or, btc_define([] if not passed)
%
% syms: a structure, values of fields indicate (0:no, 1: yes) whether each of the symmetry
% transformations preserves the block probabilities
%   rot90: rotation by 90 deg
%   rot180: rotation by 180 deg
%   hflip: flip on horizontal axis
%   vflip: flip on vertical axis
%   hvi_nwse: flip on diagonal, NW to SE axis
%   hvi_nesw: flip on diagonal, NW to SE axis
%   coninv:  exchange white and dark
%
%   See also:  BTC_DEFINE, BTC_ROTCODE,  BTC_HVI, BTC_HFLIP, BTC_VFLIP.
%
if (nargin<=1) ipn=1; end
if (nargin<=2) tol=0.0001; end
if (nargin<=3) dict=btc_define([]); end
codel=dict.codel;
%
syms=[];
for isym=1:6
    switch isym
        case 1
            flipcodes=btc_rotcode(codel,1,dict);
            symname='rot90';
        case 2
            flipcodes=btc_rotcode(codel,2,dict);
            symname='rot180';
        case 3
            flipcodes=btc_hflip(codel,dict);
            symname='hflip';
        case 4
            flipcodes=btc_vflip(codel,dict);
            symname='vflip';
        case 5
            flipcodes=btc_hvi(codel,dict);
            symname='hvi_nwse';
        case 6
            flipcodes=btc_hvirev(codel,dict);
            symname='hvi_nesw';
    end
    letcode=btc_vec2letcode(vec,dict);
    newcode=[];
    for k=1:length(codel)
        newcode.(flipcodes(k))=letcode.(codel(k));
    end
    newvec=btc_letcode2vec(newcode,dict);
    if (max(abs(newvec-ipn*vec))<tol)
        syms.(symname)=1;
    else
        syms.(symname)=0;
    end
end
%contrast inversion
newvec=vec;
oddpointers=find(mod(dict.order,2)==1);
newvec(oddpointers)=-vec(oddpointers);
if (max(abs(newvec-ipn*vec))<tol)
    syms.coninv=1;
else
    syms.coninv=0;
end
return

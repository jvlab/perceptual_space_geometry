function [proj,ptens]=btc_sym_proj_tensor(transforms,dict)
% [proj,ptens]=btc_sym_proj_tensor(transforms,dict) finds the projection 
% and the tensor that symmetrizes by one or more transformations
% in the dihedral group, possibly augmented by contrast inversions 
%
% transforms:  a cell array of strings, each of which is composed of the following:
%   H: horizontal mirror (flip top to bottom)*
%   V: vertical mirror (flip left to right)*
%   B: backslash mirror (flip upper right to lower left)*
%   S: slash mirror (flip upper left to lower right)*
%   R: rotate right*
%   L: rotate left 
%   I: inversion (180 deg rotation)*
%   E: identity (or empty)
%   C: contrast inversion
%   * these are the conventions of btc_symmetries.doc.   
%  for the dihedral group, transforms={'H','R'} (or many other combinations
%    that generate the group)
%  for rotation only, transforms={'R'}
%  for contrast symmetry only, transforms={'C'}
%  for the dihedral group and contrast symmetry, transforms={'H','R','C'}
%
%   a combination xy is interpreted as "first apply x, then apply y", so HB=R
%   case is irrelevant.
% dict: output of btc_define (created if not supplied)
%
% proj:  a projection matrix that symmetrizes a vector v by v_sym=proj*v
% ptens: a tensor (a sum of elementary symmetrizing matrices), for which a
%    matrix m is symmetrized by
%    m_sym=reshape(ptens*m(:),[length(dict.codel),length(dict.codel)])
%
%   See also: BTC_HFLIP, BTC_VFLIP, BTC_HVI, BTC_ROTCODE, BTC_D4MATRIX, BTC_MAT_SYMMETRIZE.
%
tol=10^-10;
if (nargin<=1) dict=btc_define([]); end
n=length(dict.codel);
proj=eye(n);
ptens=eye(n*n);
t=transforms;
%
for k=1:length(transforms)
    proj_term=eye(n);
    ptens_term=eye(n*n);
    npow=1;
    psym=btc_d4matrix(transforms{k},dict);
    ppow=psym;
    while max(max(abs(ppow-eye(n))))>tol
        proj_term=proj_term+ppow;
        ptens_term=ptens_term+kron(ppow',ppow'); %we're computing proj^n'MP^n
        ppow=ppow*psym;
        npow=npow+1;
    end
    %transforms{k}
    %npow
    proj_term=proj_term/npow;
    proj=proj*proj_term;
    ptens_term=ptens_term/npow;
    ptens=ptens_term*ptens;
end
return

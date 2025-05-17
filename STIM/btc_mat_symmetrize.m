function [m_sym,ztens]=btc_mat_symmetrize(m,transforms,dict)
% [m_sym,ztens]=btc_mat_symmetrize(m,transforms,dict) symmetrizes a matrix
% m by one or more transformations in the dihedral group, possibly
% augmented by contrast inversions 
%
% m: a square matrix of size=length(dict.codel)
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
% m_sym:  the matrix, symmetrized by repeated application of the group elements
% ztens: a tensor (a sum of elementary symmetrizing matrices), for which
%    m_sym=reshape(ztens*m(:),[length(dict.codel),length(dict.codel)])
%
%   See also: BTC_HFLIP, BTC_VFLIP, BTC_HVI, BTC_ROTCODE, BTC_D4MATRIX, BTC_SYM_PROJ_TENSOR.
%
if (nargin<=2) dict=btc_define([]); end
[zproj,ztens]=btc_sym_proj_tensor(transforms,dict);
n=length(dict.codel);
m_sym=reshape(ztens*m(:),n,n);
return

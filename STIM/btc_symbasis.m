function [m,desc,opts_used]=btc_symbasis(opts)
% m=btc_symbasis(opts) calculates a matrix for transforming the standard
% btc coordinates into ones that partially diagonalize the metric, based on symmetry
%
%  CAUTION:  eigenvectors calculated in the new basis need to be transformed back
%  to btc coordinates
%
% assumes standard names for coordinates
%
% opts.real_complex: 0 for real m, 1 for complex (defaults to 0)
%
% m: the matrix, 10 x 10, rows are orthonormal; m'*m=m*m'=identity
% desc: cell array, describing the entries
% opts_used: options used
%
%   See also:  BTC_EIVECS_STATS
%
if (nargin<1) opts=[]; end
opts=filldefault(opts,'real_complex',0);
opts_used=opts;
%
m=zeros(10);
desc=cell(10,1);
%invariant coords
desc{1}='gamma     invariant';
desc{2}='beta_hv   invariant';
desc{3}='beta_diag invariant';
desc{4}='theta     invariant';
desc{5}='alpha     invariant';
m( 1,1)=1;
m( 2,[2:3])=1/sqrt(2);
m( 3,[4:5])=1/sqrt(2);
m( 4,[6:9])=1/2;
m( 5,10)=1;
%coords that invert with horiz or vert flip
desc{6}='beta_diag hv_flip_invert';
desc{7}='theta     hv_flip_invert';
m( 6,[4:5])=1/sqrt(2)*[1 -1];
m( 7,[6:9])=1/2*[1 -1 1 -1];
%coords that invert with diag flip
desc{8}='beta_hv diag_flip_invert';
m( 8,[2:3])=1/sqrt(2)*[1 -1];
%coords that transform according to 2-d rep
switch opts.real_complex
    case 0
        m( 9,[6:9])=1/sqrt(2)*([1 0 -1 0]);
        m(10,[6:9])=1/sqrt(2)*([0 1 0 -1]);
        desc{9}='theta 2-dimensional cos';
        desc{10}='theta 2-dimensional sin';
    case 1
        m( 9,[6:9])=1/2*([1 i -1 -i]);
        m(10,:)=conj(m(9,:));
        desc{9}='theta 2-dimensional';
        desc{10}='theta 2-dimensional conj';
end
return


function [exptdirs,ou]=btc_token_edirs(plane,opts,dict)
% [exptdirs,ou]=btc_token_edirs(plane,opts,dict) creates a structure of directions
% and conditions used in a 4-AFC segmentation experiment for tokens
%
% note this is customized for tokens, refer to btc_edirs for standard-check coordinates
%
% plane: a 2-letter specification of the experiment, 
% always graphed with the first (horizontal) coordinate as the SECOND letter
%
% opts.layout='std' for 8 spokes: 5 vals along each axis, 2 along the diags
% opts.layout='rad' for 12 spokes, 3 vals along each direction
%    if empty, defaults to 'std' for all but tg and at, 'rad' for tg and ag
% opts.maxdist_axes=maximum values on the two coordinate axes, in plane LETTER order
% opts.nvals_eachdir: number of values along each tested direction
% opts.maxdist_diags=maximum distances on the diagonals, in plane LETTER order
%     maxdist_diags(:,1) and (:,2):  maximum values for each kind of diagonal
%     maxdist_diags(1,:) is for the diagonal closest to horizontal
%     maxdist_diags(end,:) is for the diagonal closest to vertical
%   NOTE THAT IN maxdist_diags
%     (:,1) is vertical, and corresponds to the FIRST letter in the plane specifier
%     (:,2) is horizontal, and corresponds to the SECOND letter in the plane specifier
% opts.reldist_axes: relative distances to the max, for each condition along the axes
% opts.reldist_diags: relative distances to the max, for each condition along the diagonals
% dict:  coordinate dictionary, from btc_define (or, btc_define([] if not passed)
%
% exptdirs.uvecs: unit vectors in the tested directions
% exptdirs.mvecs: maximum vectors in the tested directions
% exptdirs.allvecs{idir}:  all vectors in the tested directions
%   NOTE THAT IN exptdirs.uvecs,mvecs,allvecs,
%     (:,1) is horizontal, and corresponds to the SECOND letter in the plane specifier
%     (:,2) is vertical, and corresponds to the FIRST letter in the plane specifier
% exptdirs.allvecs_concat: concatenated allvecs{idir}
% exptdirs.opts: options used
%
%   See also:  BTC_TOKEN_XLS2DS, BTC_EDIRS, BTC_SOID_TEST, BTC_PAIRS_NEEDED.
%
if (nargin<=2) dict=btc_define([]); end
if (nargin<=1) opts=[]; end
%%%%editable section for defaults and special cases
%list of planes that default to rad layout
def_radlist={'at','tg'};
%default distances along axes
def_maxdist_axes.gamma    =0.25; %not used as of 27 Aug 15
def_maxdist_axes.beta_hv  =0.90;
def_maxdist_axes.beta_diag=0.90;
def_maxdist_axes.theta    =1.00; %not used as of 27 Aug 15
def_maxdist_axes.alpha    =1.00; %not used as of 27 Aug 15
%special cases for diagonal distances, madj=0
def_diag_special.std=[];
def_diag_special.std.bc=[0.5500 0.5500];
def_diag_special.std.de=[0.7000 0.7000];
def_diag_special.rad=[];
%
%use defaults unless other values are supplied
%layout
def_layout='std';
if ~isempty(strmatch(plane,def_radlist,'exact')) def_layout='rad'; end
opts=filldefault(opts,'layout',def_layout);
%
%look up table of maximum values
maxdist_axes=zeros(1,2);
for ic=1:2
    icn=find(plane(ic)==dict.codel);
    maxdist_axes(ic)=getfield(def_maxdist_axes,dict.name_order_aug{icn});
end
opts=filldefault(opts,'maxdist_axes',maxdist_axes);
%
switch opts.layout
    case 'std'
        nvals_eachdir=repmat([5 2],1,4);
        reldist_axes=[1:5]/5;
        reldist_diags=[0.7 1.0];
    case 'rad'
        nvals_eachdir=repmat(3,1,12);
        reldist_axes=[1:3]/3;
        reldist_diags=[1:3]/3;
end
opts=filldefault(opts,'nvals_eachdir',nvals_eachdir);
opts=filldefault(opts,'reldist_axes',reldist_axes);
opts=filldefault(opts,'reldist_diags',reldist_diags);
ndirs=size(opts.nvals_eachdir,2);
ndkinds=(ndirs/4-1);
%
% set up defaults in diagonal directions as arith sequence, unless there are exceptions
%   for std: maxdist_diags=mean(opts.maxdist_axes) on each coord
%   for rad: maxdist_diags=[1 2;2 1]*diag(opts.maxdist_axes)/3
%
% this is tricky since maxdist_axes and maxdist_diags are in plane letter
% order, so that the entry that is most heavily weighted for the second
% coordinate is most nearly horizontal.  So
% maxdist_diags  (1,:) is most horizontal
% maxdist_diags(end,:) is most vertical
% maxdist_diags(:,2) is the x-axis coordinate
% maxdist_diags(:,1) is the y-axis coordinate
% 
maxdist_diags=[[1:ndkinds];fliplr([1:ndkinds])]'*diag(opts.maxdist_axes)/(ndkinds+1);
%check for exceptions
diag_special=getfield(def_diag_special,opts.layout);
if (isfield(diag_special,plane))
    maxdist_diags=getfield(diag_special,plane);
end
opts=filldefault(opts,'maxdist_diags',maxdist_diags);
%
%calculate table of directions
%
exptdirs=[];
exptdirs.plane=plane;
exptdirs.ndirs=ndirs;
exptdirs.ndiagonal_kinds=ndkinds;
%create max vectors, working counterclockwise from x-axis
exptdirs.allvecs_concat=[];
for idir=1:ndirs
    angle=2*pi*(idir-1)/ndirs;
    %set ca and sa to be -1, 0, o+1
    ca=sign(round(10000*cos(angle)));
    sa=sign(round(10000*sin(angle)));
    if mod(idir-1,(ndkinds+1))==0 %axis
        mvec=[ca*opts.maxdist_axes(2),sa*opts.maxdist_axes(1)];
        fdists=opts.reldist_axes';
    else
        which_diag=mod(idir-1,(ndkinds+1));
        if ~(sa==ca) which_diag=ndkinds+1-which_diag; end
        mvec=[ca*opts.maxdist_diags(which_diag,2),sa*opts.maxdist_diags(which_diag,1)];
        fdists=opts.reldist_diags';
    end
    mvec_mag=sqrt(sum(mvec.^2));
    uvec=mvec/mvec_mag;
    exptdirs.uvecs(idir,:)=uvec; %unit vector in this direction
    exptdirs.maxvecs(idir,:)=mvec; %maximum vector in this direction
    exptdirs.allvecs{idir}=fdists*mvec;
    exptdirs.allvecs_concat=[exptdirs.allvecs_concat;exptdirs.allvecs{idir}];
end
exptdirs.opts=opts;
%
ou=opts;
return

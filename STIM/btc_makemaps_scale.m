function [imgs,optsused,errs,metro]=btc_makemaps_scale(method,opts,dict,auxopts,mtcopts)
% [imgs,optsused,errs,metro]=btc_makemaps_scale(method,opts,dict,auxopts,mtcopts)
% makes one or more,
% textures based on a "method" (from btc_augcoords), with provision checks of arbitrary scale
% 
% here, scale=n means that the local image statistics refer to nxn subblocks of imgs
%
% input arguments are identical to those of btc_makemaps_scale except for additional fields in opts
%
% method: structure returned by btc_augcoords or mtc_augcoords, see btc_makemaps
% opts: as in btc_makemaps with the following additional fields and exceptions:
%    opts.blocked: defaults to 1, makes textures with blocks of given size
%    opts.staggered: defaults to 1, makes textures with blocks of given size, offset by single pixels
%        blocked and staggered can be superimposed
%        blocked*staggered must divide image size
%    opts.firstrow, opts.firstcol, opts.lastcol: ignored, if opts.scale_method='staggered'
% dict: as in btc_makemaps.
% auxopts: as in btc_makemaps
% mtcopts: as in btc_makemaps
%
% output:
% imgs: the images requested, [0: ng-1] if opts.scale_method='staggered', or [0:(ng-1)*(opts.staggered^2)] for opts.staggered>1
% optsused: as in btc_makemaps, but cell array of size {opts.staggered, opts.staggered} 
% errs: as in makemaps, but cell array of size {opts.staggered, opts.staggered} 
% metro: as in makemaps, but cell array of size {opts.staggered, opts.staggered} 
%
% number of gray levels determined from dimension of p2x2 in methods.
% if ng=2, methods may have been generated either by BTC_AUGCOORDS or MTC_AUGCOORDS,
%   so other fields of method are inspected for the non-Pickard methods
%
%   See also BTC_AUGCOORDS, BTC_DEFINE, BTC_TEST, BTC_SURV, GENMRFM,
%   GENMRFM_1D, MAPUBI, BTC_DIMRF, BTC_NOPICK, DONUT_METRO, MTC_FLFS_MAKEMAPS, RFP_BTC_MAKE, BTC_MAKEMAPS_SCALE_TEST.
%
%

if (nargin<=4)
    mtcopts=[];
end
mtcopts=mtc_defopts(mtcopts);
if (nargin<=3)
    auxopts=btc_auxopts([],[]);
end
if ~isfield(auxopts,'metro_opts')
    auxopts_std=btc_auxopts;
    auxopts.metro_opts=auxopts_std.metro_opts;
    auxopts.metro_show=auxopts_std.metro_show;
end
%
if (nargin<=2)
    dict=btc_define();
end
if (nargin<=1)
    opts=[];
end
ng=size(method.p2x2,1); %number of gray levels
%
opts=filldefault(opts,'nmaps',1);
opts=filldefault(opts,'show',[]);
opts=filldefault(opts,'area',[256 256]);
opts=filldefault(opts,'firstrow',[]);
opts=filldefault(opts,'firstcol',[]);
opts=filldefault(opts,'lastcol',[]);
opts=filldefault(opts,'verbose',0);
opts=filldefault(opts,'onaxis',0); %1 if on axis (and bypasses need for Metropolis)
if length(opts.area)==1
    opts.area=repmat(opts.area,1,2);
end
imgs=zeros([opts.area opts.nmaps]);
imgs_staggered=zeros([opts.area/opts.blocked, opts.nmaps]);
opts=filldefault(opts,'blocked',1);
opts=filldefault(opts,'staggered',1);
%
nsteps=opts.staggered;
optsused=cell(nsteps,nsteps);
errs=cell(nsteps,nsteps);
metro=cell(nsteps,nsteps);
opts_staggered=opts;
area_staggered=opts.area/opts.blocked/nsteps;
opts_staggered.area=area_staggered+1;
for istep=1:nsteps
    for jstep=1:nsteps
        [im,optsused{istep,jstep},errs{istep,jstep},metro{istep,jstep}]=btc_makemaps(method,opts_staggered,dict,auxopts,mtcopts);
        im_big=repblk(im,[nsteps nsteps 1]);
        imgs_staggered=imgs_staggered+im_big(istep+[0:(opts.area(1)/opts.blocked-1)],jstep+[0:(opts.area(2)/opts.blocked-1)],:);
    end
end
imgs=repblk(imgs_staggered,[opts.blocked opts.blocked 1]);
return


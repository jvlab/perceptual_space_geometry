function [angles,mults,angles_stats,mults_stats,rayfit,opts_used]=psg_raystats(coords,sa,rays,jit_rms,opts)
% [angles,mults,angles_stats,mults_stats,rayfit,opts_used]=psg_raystats(coords,sa,rays,jit_rms,opts)
% computes angles between rays and gains along each ray from a set of coordinates in a psg experiment
% along with surrogate values obtained by jittering the coordinates
%
% method is repeated calls to psg_rayfit, psg_rayangles, psg_raymults, 
%  two modes, determined by opts.if_bid (i psg_rayfit)
%    np=2 (if_bid=0): each polarity considered separately [neg, pos]
%    np=1 (if_bid=1): bidirectional calculation, neg and positive rays considered as part of a line segment
%
% coords: coordinate array [nstims nd], typically a field of d as returned by psg_readcoords
% sa: the setup structure returned from psg_readcoord_data (sa.typenames needed for labeling)
% rays: a ray structure, typically returned by psg_findrays
% jit_rms: root-mean-squared jitter on each coordinate axis, 
%     values typically obtained from psg_lljit_crit, with jit_type='shell'
% opts: options
%   opts.if_bid: 0 for unidirectional rays, 1 for bidirectional
%   opts.if_origin:
%      0 to not consider origin in fit but allow offset (need at least two points per ray)
%      1 to include as regressor and allow offset (default for if_bid=1)
%     -1 to include and force ray to emanate from the origin (default for if_bid=0)
%   opts.nsurrs: number of surrogates, defaults to 100.
%   opts.if_log: whether to log, applieas only to non-surrogate, defaults to 0
%   opts.p: p-value used for bootstrap calc of confidence limits, defaults
%      to 'extreme', as jit_rms already takes the p-value into account
%
% angles: see psg_rayangles
% mults: see psg_raymults.
% % rayfit: fitted coordinate array, corresponding to coords (see psg_rayfit)
%     * origin is always mapped to orgin, regardless of opts.if_origin
% opts_used: options used
%   opts_used.psg_rayfit: options returned by psg_rayfit with no jitter
%   opts_used.psg_rayangles: options returned by psg_rayangles with no jitter
%   opts_used.psg_raymults: optoins returned by psg_raymults with no jitter
%   opts_used.surrs.angles: structure array of the opts.nsurrs surrogates for angles
%   opts_used.surrs.mults: structure array of the opts.nsurrs surrogates for mults
% 
%   See also:  PSG_READCOORD_DATA, PSG_RAYFIT, FILLDEFAULT, PSG_RAYANGS,
%     PSG_RAYMULTS, PSG_LLJIT_CRIT, BOOTSA.
%
if (nargin<5)
    opts=struct;
end
opts=filldefault(opts,'if_log',0);
opts=filldefault(opts,'if_bid',0);
opts=filldefault(opts,'if_origin',-1+2*opts.if_bid);
opts=filldefault(opts,'nsurrs',100);
opts=filldefault(opts,'p','extreme');
opts_sub=opts;
opts_sub=rmfield(opts_sub,'nsurrs');
opts_sub=rmfield(opts_sub,'p');
opts_sub2=opts_sub;
opts_sub2=rmfield(opts_sub2,'if_bid');
opts_sub2=rmfield(opts_sub2,'if_origin');
%
opts_used=opts;
%
% original data
%
[rayfit,ray_ends,opts_used.psg_rayfit]=psg_rayfit(coords,rays,rmfield(opts_sub,'if_log'));
[angles,opts_used.psg_rayangles]=psg_rayangles(ray_ends,sa,rays,opts_sub2);
[mults ,opts_used.psg_raymults ]=psg_raymults(ray_ends,sa,rays,opts_sub2);
%
%surrogates
%
jits=jit_rms*randn(size(coords,1),size(coords,2),opts.nsurrs);
%
a=struct; %for angles
angles_fns=fieldnames(angles);
m=struct; %for mults
mults_fns=fieldnames(mults);
for isurr=1:opts.nsurrs
    [rayfit_surr,ray_ends_surr]=psg_rayfit(coords+jits(:,:,isurr),rays,opts);
    %
    angles_surr=psg_rayangles(ray_ends_surr,sa,rays,setfield(opts,'if_log',0));
    for ifn=1:length(angles_fns)
        fn=angles_fns{ifn};
        a(isurr).(fn)=angles_surr.(fn);
    end
    %
    mults_surr=psg_raymults(ray_ends_surr,sa,rays,setfield(opts,'if_log',0));
    for ifn=1:length(mults_fns)
        fn=mults_fns{ifn};
        m(isurr).(fn)=mults_surr.(fn);
    end
end
opts_used.surrs.angles=a;
opts_used.surrs.mults=m;
%
[bbias,bdebiased,bvar,bsem,blower,bupper]=bootsa(angles,a,opts.p);
angles_stats.sem=bsem;
angles_stats.clo=blower;
angles_stats.chi=bupper;
%
[bbias,bdebiased,bvar,bsem,blower,bupper]=bootsa(mults,m,opts.p);
mults_stats.sem=bsem;
mults_stats.clo=blower;
mults_stats.chi=bupper;
%
return

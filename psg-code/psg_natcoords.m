function [natcoords,opts_used]=psg_natcoords(coords,sa,rays,opts)
% [natcoords,opts_used]=psg_tangents(coords,sa,rays,opts)
% determines natural coordinates in a psg dataset based on tangents
% to the trajectory of each stimulus-axis locus and similar quantiites.
%
% tangents are computed from deriviatives on a spline
%
% coords: coordinate array [nstims nd], typically a field of d as returned by psg_get_coordsets or psg_readcoords
% sa: metadata structure, typcially returned by psg_get_coordsets
% rays: a ray structure, typically returned by psg_findrays
%    nrays=rays.nrays, number of rays
% opts: options
%   opts.if_origin:
%      0 to assume origin is at zero
%      1 to look for origin via rays.whichray==0 (default, same as psg_rayfit for if_bid=1)
%   opts.minpts: minimum number of points on a locus (defaults to 5)
%   opts.dfrac: fraction of minimum separation of trajectory values to use for numerical derivative
%
% note that for btc stimuli, "origin" refers to the coordinates of the random texture, which may not be zero
%
% natcoords: struture with fields:
%    labels (nrays,1:) cell array of coordinate labels
%    tangents (nrays, nd): tangent vectors to each locus; length is the speed
%       (derivatives taken with respect to rays.mult)
%    t_traj {nrays,1}: list of parameter values along the ray, includes mult and endpoint
%    coords_traj {nrays,1}: (length(t{iray}),nd), coords along the ray
%    fit_line (nrays, nd): comparable to tangents, but is rms best fit along entire trajectory
%    ray_ends( nrays, nd): rms best fit to bidirectional ray, from psg_rayfit, used to calculate fit_line
%    avail: cell array of names of fields with natural coordinates, typically {'fit_line','tangent'} if all can be computed
%
%   Normalization:
%   natcoords.tangents and natcoords.fit_line are normalized the same way:  their length is the
%   distance moved for a unit amount of distance along the stimulus axis.
%   (tangents: determined by a derivative and extrapolating; fit_line: determined by rms best fit)
%   Distance is determined by the distance to the trajectory endpoint, ray.endpt.
%   Coordinates used in experiment are determined by multiplying that endpoint by rays.mult, 
%   and then augmenting via maximum-entropy.
%
% opts_used: options used
%
%   See also:  FILLDEFAULT, PSG_VISUALIZE_DEMO, PSG_FINDRAYS, PSG_RAYFIT, PSG_NATCOORDS_DEMO, SPLINE.
%
if (nargin<3)
    opts=struct;
end
opts=filldefault(opts,'if_origin',1);
opts=filldefault(opts,'minpts',5);
opts=filldefault(opts,'dfrac',0.1);
%
origin_ptr=find(rays.whichray==0);
if isempty(origin_ptr)
    origin=zeros(1,nd);
    if opts.if_origin~=0
        warning('origin requested but not found');
    end
    opts.origin=[];
else
    origin=coords(origin_ptr,:);
    opts.origin=origin; 
end
opts_used=opts;
%
npts=size(coords,1);
nd=size(coords,2);
nrays=rays.nrays;
%
natcoords=struct();
natcoords.avail=cell(0);
%
%invoke psg_rayfit to find best-fitting line
opts_rayfit=struct;
opts_rayfit.if_origin=opts.if_origin; 
opts_rayfit.if_bid=1; %always bidirectional
[rayfit,ray_ends,opts_rayfit_used]=psg_rayfit(coords,rays,opts_rayfit);
natcoords.ray_ends=ray_ends;
%
endpt_dists=ones(nrays,1); %total correlatoin at endpoint of a ray
if isfield(rays,'endpt')
    for iray=1:nrays
        endpt_dists(iray)=sqrt(sum(rays.endpt(iray,:).^2));
    end
else
    warning('endpt not present in ray structure');
end
natcoords.fit_line=natcoords.ray_ends./repmat(endpt_dists,1,nd);
natcoords.avail{end+1}='fit_line';
%
natcoords.tangents=NaN(nrays,nd);
natcoords.labels=cell(nrays,1);
natcoords.t_traj=cell(nrays,1);
natcoords.coords_traj=cell(nrays,1);
%
%
%find tangents and adjust fit_line
%
for iray=1:nrays
    nz_ptrs=find(rays.whichray==iray);
    if opts.if_origin
        cs=coords([origin_ptr;nz_ptrs(:)],:);
    else
        cs=[zeros(1,nd);coords(nz_ptrs,:)];
    end
    t_vals_nz=rays.mult(nz_ptrs(:));
    if isfield(rays,'endpt')
        t_vals_nz=t_vals_nz*endpt_dists(iray);
    end
    t_max_ptr=find(t_vals_nz==max(t_vals_nz));
    t_vals=[0;t_vals_nz];
    natcoords.labels{iray}=sa.typenames{nz_ptrs(t_max_ptr)}; %label is typenames for largest value on axis
    %extract coordinates on this axis, and sort by multiplier
    vc=sortrows([t_vals,cs]);
    t=vc(:,1);
    coords_traj=vc(:,2:end);
    natcoords.t_traj{iray}=t;
    natcoords.coords_traj{iray}=coords_traj;
    %compute tangent by taking deriv on each coordinate, via spline
    if size(cs,1)>=opts.minpts
        dt=opts.dfrac*min(abs(diff(t)));
        for id=1:nd
            coords_near_z=spline(t,coords_traj(:,id),dt*[-.5 .5]);
            natcoords.tangents(iray,id)=diff(coords_near_z)/dt;
        end
    end
end
if any(~isnan(natcoords.tangents(:)))
    natcoords.avail{end+1}='tangents';
end
return
end

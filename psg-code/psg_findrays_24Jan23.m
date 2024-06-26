function [rays,opts_used]=psg_findrays(stim_coords,opts)
% [rays,opts_used]=psg_findrays(stim_coords,opts) finds and describes the rays in a set of coordinates
%
% stim_coords: the coordinates, size [npts nd], e.g., sa.btc_specoords from psg_read_coorddata
%   Note that sa.btc_augcoords could also be used, but these values include the induced correlations
%   so, for example, increasing values of gamma or beta will not lie on a ray
% opts: options structure, can be omitted
%   opts.tol: tolerance for matching points or angle wrapping, defaults to 10^-5
%   opts.res_ring: resolution for finding points with equal magnitudes for rings, defaults to 10^-2
%     res_ring is lax because coordinates are derived from stimulus descriptions, which have few significant digits
%   opts.min_ring: minimum number of points in a ring, defaults to 4
%     for example, to make the first ray found to be ray number x and the second ray found
%     to have ray number y, then opts.permute_raynums=[x y ...]
%   opts.permute_raynums: a permutation to apply to ray numbers and endpoints for nonzero rays
%
% rays: structure
%   rays.nrays: number of rays, other than the 0 ray
%   rays.whichray, size [npts,1]: whichray(ipt) is the ray that stim_coords(ipt,:) is asssigned to
%   rays.mult, size [npts,1]: the multiplier
%   rays.endpt, size [nrays,1]: 
% opts_used: options used
%
% 24Jan23: add rings. Note that code is not fully tested since existing
%   datasets all have 'rand' condition at end, so fnz is [1:nstims-1], and
%   does not re-order the coordinates
%
%   See also:  PSG_READ_COORDDATA, FILLDEFAULT, PSG_VISUALIZE_DEMO, PSG_PLOTCOORDS, 
%   PSG_PLANECYCLE,PSG_PCAOFFSET.
%
if (nargin<2)
    opts=struct;
end
opts=filldefault(opts,'tol',10^(-5));
opts=filldefault(opts,'res_ring',10^(-2));
opts=filldefault(opts,'min_ring',4);
opts=filldefault(opts,'permute_raynums',[]);
opts_used=opts;
npts=size(stim_coords,1);
nd=size(stim_coords,2);
%
rays=struct;
%
stim_coords(isnan(stim_coords))=0;
%
unassigned=[1:npts];
mags=sqrt(sum(stim_coords.^2,2));
zero_length=find(mags<=opts.tol);
rays.whichray=zeros(npts,1);
rays.mult=zeros(npts,1);
rays.endpt=[];
unassigned=setdiff(unassigned,zero_length);
iray=0;
while ~isempty(unassigned)
    iray=iray+1;
    v=stim_coords(min(unassigned),:);
    projs=stim_coords(unassigned,:)*v'/(v*v');
    diffs=stim_coords(unassigned,:)-projs*v;
    mags=sqrt(sum(diffs.^2,2));
    matches=find(mags<=opts.tol);
    %
    rays.whichray(unassigned(matches))=iray;
    maxproj=max(projs(matches));
    maxproj_ind=min(find(projs==maxproj));
    rays.endpt(iray,:)=stim_coords(unassigned(maxproj_ind),:);
    rays.mult(unassigned(matches))=projs(matches)/maxproj;
    %
    unassigned=setdiff(unassigned,unassigned(matches));
end
rays.nrays=length(unique(setdiff(rays.whichray,0)));
fnz=find(rays.whichray~=0); %all points that are on a ray (typically, all points except the origin)
if ~isempty(opts.permute_raynums)
    rays.whichray(fnz)=opts.permute_raynums(rays.whichray(fnz));
    rays.endpt(opts.permute_raynums,:)=rays.endpt;
end
%find rings
%find progression of angles for all points that are on rays
stim_coords_mod=stim_coords;
stim_coords_mod(find(isnan(stim_coords)))=0;
[recon_pcplane,cyclic_order]=psg_planecycle(stim_coords_mod(fnz,:),opts); 
%find groups of points that are at the same distance from the origin
mult_vals=round(abs(rays.mult(fnz(cyclic_order))./opts.res_ring));
mult_vals_unique=flipud(unique(mult_vals)); %mult_vals_unique begins with largest value
nrings=0;
rings=cell(0);
for iring=1:length(mult_vals_unique)
    thisring=find(mult_vals==mult_vals_unique(iring));
    if length(thisring)>=opts.min_ring
        rings{iring}.mult_val=opts.res_ring*mult_vals_unique(iring);
        rings{iring}.coord_ptrs=fnz(cyclic_order(thisring))';
    end
end
%
%alternate code that looks at one ring at a time -- 
%simpler, but the rings may start at different positions because of
%degeneracy of PCA
%
% %find groups of points that are at the same distance from the origin
% mult_vals=round(abs(rays.mult(fnz)./opts.res_ring));
% mult_vals_unique=flipud(unique(mult_vals)); %mult_vals_unique begins with largest value
% nrings=0;
% rings=cell(0);
% for iring=1:length(mult_vals_unique)
%     thisring=find(mult_vals==mult_vals_unique(iring));
%     if length(thisring)>=opts.min_ring
%         %if this were applied to all rings together, then all rings would start in the same place
%         [recon_pcplane,cyclic_order]=psg_planecycle(stim_coords_mod(fnz(thisring),:),opts);
%         rings{iring}.mult_val=opts.res_ring*mult_vals_unique(iring);
%         rings{iring}.coord_ptrs=fnz(thisring(cyclic_order))';
%     end
% end
rays.rings=rings;
return

function [rays,opts_used]=psg_findrays(coords,opts)
% [rays,opts_used]=psg_findrays(coords,opts) finds and describes the rays in a set of coordinatas
%
% coords: the coordinates, size [npts nd], e.g., sa.btc_augcoords or sa.btc_specoords
%   from psg_read_coorddata
% opts: options structure, can be omitted
%   opts.tol: tolerance, defaults to 10^-5
%   opts.permute_raynums: a permutation to apply to ray numbers and endpoints for nonzero rays
%     for example, to make the first ray found to be ray number x and the second ray found
%     to have ray number y, then opts.permute_raynums=[x y ...]
%
% rays: structure
%   rays.nrays: number of rays, other than the 0 ray
%   rays.whichray, size [npts,1]: whichray(ipt) is the ray that coords(ipt,:) is asssigned to
%   rays.mult, size [npts,1]: the multiplier
%   rays.endpt, size [nrays,1]: 
% opts_used: options used
%
%   See also:  PSG_READ_COORDDATA, FILLDEFAULT, PSG_VISUALIZE_DEMO, PSG_PLOTCOORDS.
%
if (nargin<2)
    opts=struct;
end
opts=filldefault(opts,'tol',10^(-5));
opts=filldefault(opts,'permute_raynums',[]);
opts_used=opts;
npts=size(coords,1);
nd=size(coords,2);
%
rays=struct;
%
coords(isnan(coords))=0;
%
unassigned=[1:npts];
mags=sqrt(sum(coords.^2,2));
zero_length=find(mags<=opts.tol);
rays.whichray=zeros(npts,1);
rays.mult=zeros(npts,1);
rays.endpt=[];
unassigned=setdiff(unassigned,zero_length);
iray=0;
while ~isempty(unassigned)
    iray=iray+1;
    v=coords(min(unassigned),:);
    projs=coords(unassigned,:)*v'/(v*v');
    diffs=coords(unassigned,:)-projs*v;
    mags=sqrt(sum(diffs.^2,2));
    matches=find(mags<=opts.tol);
    %
    rays.whichray(unassigned(matches))=iray;
    maxproj=max(projs(matches));
    maxproj_ind=min(find(projs==maxproj));
    rays.endpt(iray,:)=coords(unassigned(maxproj_ind),:);
    rays.mult(unassigned(matches))=projs(matches)/maxproj;
    %
    unassigned=setdiff(unassigned,unassigned(matches));
end
rays.nrays=length(unique(setdiff(rays.whichray,0)));
if ~isempty(opts.permute_raynums)
    fnz=find(rays.whichray~=0);
    rays.whichray(fnz)=opts.permute_raynums(rays.whichray(fnz));
    rays.endpt(opts.permute_raynums,:)=rays.endpt;
end
return



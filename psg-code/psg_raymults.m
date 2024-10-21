function [mults,opts_used]=psg_raymults(ray_ends,sa,rays,opts)
% [mults,opts_used]=psg_raymults(ray_ends,sa,rays,opts) computes multipliers for
% lengths of rays from a set of coordinates in a btc experiment
%
%  the "multiplier" is the length, in coordinate units, of the linear extrapolation of a 
%  best-fitting ray to a correlation strength of 1
%
%  two modes, determined by np=size(ray_ends,2) (if_bid=2-np))
%    np=2 (if_bid=0): each polarity considered separately [neg, pos]
%    np=1 (if_bid=1): bidirectional calculation, neg and positive rays considered as part of a line segment
%     and as part of a single regression
%
% ray_ends: regression coefs at ends of each ray [nrays nd 2-if_bid], typically from psg_rayfit
%    dim 3 is neg then pos for unidirectional rays (if_bid=0)
%    dim 3 is length 1 for bidirectional rays (if_bid=1)
% sa: the setup structure returned from psg_readcoord_data
%    (sa.typenames needed for labeling but dummy names used if sa is empty)
% rays: a ray structure, typically returned by psg_findrays
% opts: options
%  opts.if_log=1: to log results
%
% mults: each field is dim 1 is ray, dim 2 is [neg, pos] if np=2, trivial if np=1,
%   dist_gain: distance  for a unit change in correlation strength
%   dist_raw:  %raw distances to largest correlation strength available
%   dist_endpts: %distances to endpoints (dist_mult=dist_raw./repmat(dist_endpts,1,np)
% opts_used: options used
%
% 16Oct24: allow for sa to be empty 
%
%   See also:  PSG_READCOORD_DATA, PSG_RAYFIT, FILLDEFAULT, PSG_VISUALIZE_DEMO, PSG_FINDRAYS, PSG_RAYANGLES.
%
if (nargin<4)
    opts=struct;
end
opts=filldefault(opts,'if_log',0);
if isempty(sa) %dummy typenames if sa not supplied
    sa=struct;
    sa.typenames=cell(length(rays.whichray),1);
    for k=1:length(sa.typenames)
        sa.typenames{k}=sprintf('type %2.0f',k);
    end
end
%
nrays=size(ray_ends,1);
nd=size(ray_ends,2);
%if_bid=2-size(ray_ends,3);
np=size(ray_ends,3); %number of polarities
%
opts_used=opts;
dotprods=zeros(nrays,nrays,np,np);
dists=zeros(nrays,np); %raw distances to endpoints (np=2) or distances between endpoints (np=1)
labels=cell(nrays,np);
for iray=1:nrays
    for ip=1:np
        dists(iray,ip)=sqrt(sum(ray_ends(iray,:,ip).^2)); 
    end %ip
end %iray
mults.dist_raw=dists;
%normalize by endpoints
dist_endpts=sqrt(sum(rays.endpt.^2,2));
mults.dist_endpts=dist_endpts;
mults.dist_gain=dists./repmat(dist_endpts,1,np);
sign_select=cell(np,1);
if np==1
    sign_select{1}=find(rays.mult~=0);
else
    sign_select{1}=find(rays.mult<0);
    sign_select{2}=find(rays.mult>0);
end
%
if (opts.if_log)
    disp(sprintf('distances along rays in a model with %3.0f dimensions',nd));
end
for iray=1:nrays %section retains logic of psg_rayangs
    str=[];
    for ip=1:np
        ipoints=intersect(find(rays.whichray==iray),sign_select{ip});
        if isempty(ipoints) | (ipoints(end)<1)
            ilab='?';
        else
            ilab=sa.typenames{ipoints(end)};
        end
        %
        labels{iray,ip}=sprintf('%12s',ilab);
        str=cat(2,str,labels{iray,ip},sprintf(': %7.3f  (max coord dist: %7.3f)  ',mults.dist_gain(iray,ip),mults.dist_endpts(iray)));
    end %ip
    if (opts.if_log)
        disp(str);
    end
end %iray
opts_used.labels=labels;
return

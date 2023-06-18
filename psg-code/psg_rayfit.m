function [rayfit,ray_ends,opts_used]=psg_rayfit(coords,rays,opts)
% [rayfit,ray_ends,opts_used]=psg_rayfit(coords,rays,opts) creates a
% coordinate structure of best-fitting uni- or bidirectional rays
%   bidirectional rays include neg and postive multipliers
%   unidirectional rays: one multiplier for each sign
%
% coords: coordinate array [nstims nd], typically a field of d as returned by psg_readcoords
% rays: a ray structure, typically returned by psg_findrays
% opts: options
%   opts.if_bid: 0 for unidirectional rays, 1 for bidirectional
%   opts.if_origin:
%      0 to not consider origin in fit but allow offset (need at least two points per ray)
%      1 to include as regressor and allow offset (default for if_bid=1)
%     -1 to include and force ray to emanate from the origin (default for if_bid=0)
%   notes:
%      * origin is not necessarily at zero, but is given by coords(find(rays.whichray==0),:)
%      * rayfit always maps origin to origin
%
% note that "origin" refers to the coordinates of the random texture, which
% may not be zero
%
% rayfit: fitted coordinate array, corresponding to coords
%     * origin is always mapped to orgin, regardless of opts.if_origin
% ray_ends: array of size [nrays nd 2-if_bid], if if_bid=0, dim 3 is neg then pos
%   these are the regression coefficients, i.e., the distance from a contrast of 0 to a contrast
%   at the end of the ray
% opts_used: options used
%
% 17Jun22: protect from regressing if rays are too short
%
%   See also:  FILLDEFAULT, PSG_VISUALIZE_DEMO, PSG_FINDRAYS, PSG_RAYANGLES.
%
if (nargin<3)
    opts=struct;
end
opts=filldefault(opts,'if_bid',0);
opts=filldefault(opts,'if_origin',-1+2*opts.if_bid);
%
npts=size(coords,1);
nd=size(coords,2);
nrays=rays.nrays;
%
rayfit=zeros(npts,nd);
ray_ends=zeros(nrays,nd,2-opts.if_bid);
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
    rayfit(origin_ptr,:)=origin;
    opts.origin=origin; 
end
opts_used=opts;
%
%regression fit to each ray or bidirectional ray
%
if (opts.if_bid)
    sign_lims=0; %fit the two directions together
else
    sign_lims=[1 2]; %fit each direction separately
end
sign_msg={' ',' in neg direction',' in pos direction'};
for iray=1:nrays
    for isign=sign_lims
        allpoints_no=find(rays.whichray==iray); %all points on this ray, excluding the origin
        switch isign
            case 0
                is=1;
            case 1
                allpoints_no=intersect(allpoints_no,find(rays.mult<0));
                is=1;
            case 2
                allpoints_no=intersect(allpoints_no,find(rays.mult>0));
                is=2;
        end
        allpoints_no=allpoints_no(:); %a column
        %regress on each coordinate
        switch opts.if_origin
            case 0 %ignore origin
                if length(allpoints_no)<=1
                    warning(sprintf('ray %2.0f%s has only %2.0f points; fit option needs at least 2',...
                        iray,sign_msg{isign+1},length(allpoints_no)));
                    if length(allpoints_no)>0
                        rayfit(allpoints_no,:)=coords(allpoints_no,:);
                        ray_ends(iray,:,is)=coords(allpoints_no(end),:)./rays.mult(allpoints_no(end));
                    end
                else
                    for id=1:nd
                        x=[rays.mult(allpoints_no),ones(length(allpoints_no),1)]; %include an offset term for the ray
                        y=coords(allpoints_no,id);
                        b=regress(y,x);
                        rayfit(allpoints_no,id)=x*b;
                        ray_ends(iray,id,is)=b(1);
                    end
                end
            case 1
                allpoints=[allpoints_no(:);origin_ptr]; %include origin as a regressor
                if length(allpoints)<=1
                    warning(sprintf('ray %2.0f%s has only %2.0f points; fit option needs at least 2',...
                        iray,sign_msg{isign+1},length(allpoints)));
                    if length(allpoints_no)>0
                        rayfit(allpoints_no,:)=coords(allpoints_no,:);
                        ray_ends(iray,:,is)=coords(allpoints_no(end),:)./rays.mult(allpoints_no(end));
                    end
                else
                    for id=1:nd
                        x=[rays.mult(allpoints),ones(length(allpoints),1)]; %include an offset term for the ray
                        y=coords(allpoints,id);
                        b=regress(y,x);
                        rayfit(allpoints_no,id)=x(1:end-1,:)*b; %do not fit origin
                        ray_ends(iray,id,is)=b(1);
                    end
                end
            case -1
                if length(allpoints_no)<=1
                    warning(sprintf('ray %2.0f%s has only %2.0f points; fit option needs at least 2',...
                        iray,sign_msg{isign+1},length(allpoints_no)));
                    if length(allpoints_no)>0
                        rayfit(allpoints_no,:)=coords(allpoints_no,:);
                        ray_ends(iray,:,is)=coords(allpoints_no(end),:)./rays.mult(allpoints_no(end));
                    end
                else
                    for id=1:nd
                        x=rays.mult(allpoints_no); %no offset term
                        y=coords(allpoints_no,id)-origin(id);
                        b=regress(y,x);
                        rayfit(allpoints_no,id)=origin(id)+x*b; %force to emanate from through the origin
                        ray_ends(iray,id,is)=b;
                    end
                end
        end %opts.if_origin
    end %isign
end %iray
%
return

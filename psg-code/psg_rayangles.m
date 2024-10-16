function [angles,opts_used]=psg_rayangles(ray_ends,sa,rays,opts)
% [angles,opts_used]=psg_rayangles(ray_ends,sa,rays,opts) computes angles between rays
% from a set of coordinates in a psg experiment
%
%  two modes, determined by np=size(ray_ends,2) (if_bid=2-np))
%    np=2 (if_bid=0): each polarity considered separately [neg, pos]
%    np=1 (if_bid=1): bidirectional calculation, neg and positive rays considered as part of a line segment
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
% angles: each field is [isign, jsign point to neg and pos if if_bid=0, or are =1 if_bid=1]
%   dotprods: dot products, (iray,jray,isign,jsign)
%   cosangs:  cosine of angles (iray,jray,isign,jsign)
%   ang_degs: anlges in degress (iray,jray, isign, jsign)
% opts_used: options used
%
% 15Oct24: immprove documentation, labels created even if no logging; improve label justification
% 16Oct24: allow for sa to be empty 
%
%   See also:  PSG_READCOORD_DATA, PSG_RAYFIT, FILLDEFAULT, PSG_VISUALIZE_DEMO, PSG_FINDRAYS, PSG_PLOTANGLES
%   PSG_RAYMULTS.
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
cosangs=zeros(nrays,nrays,np,np);
dotprods=zeros(nrays,nrays,np,np);
mags=zeros(nrays,np);
labels=cell(nrays,nrays,np.^2);
for iray=1:nrays
    for ip=1:np
        mags(iray,ip)=sqrt(sum(ray_ends(iray,:,ip).^2));
    end %ip
end %iray
sign_select=cell(np,1);
if np==1
    sign_select{1}=find(rays.mult~=0);
else
    sign_select{1}=find(rays.mult<0);
    sign_select{2}=find(rays.mult>0);
end
%note that we flip the coord sign for a negative ray (if_bid=0,ip=1)
for iray=1:nrays
    for ip=1:np
        isign=1;
        if (ip==1 & np==2)
            isign=-1;
        end
        ivec=isign*ray_ends(iray,:,ip);
        for jray=1:nrays
            for jp=1:np
                jsign=1;
                if (jp==1 & np==2)
                    jsign=-1;
                end
                jvec=jsign*ray_ends(jray,:,jp);
                dotprods(iray,jray,ip,jp)=sum(ivec.*jvec);
                cosangs(iray,jray,ip,jp)=dotprods(iray,jray,ip,jp)./(mags(iray,ip)*mags(jray,jp));
            end %jp
        end %jray
    end %ip
end %iray
angles.dotprods=dotprods;
angles.cosangs=cosangs;
angles.ang_degs=acos(max(min(cosangs,1),-1))*(180/pi);
%
if (opts.if_log)
    disp(sprintf('angles (deg) between rays in a model with %3.0f dimensions',nd));
end
for iray=1:nrays
    for jray=1:nrays
        str=[];
        for ip=1:np
            for jp=1:np
                ipoints=intersect(find(rays.whichray==iray),sign_select{ip});
                ilab=sa.typenames{ipoints(end)};
                jpoints=intersect(find(rays.whichray==jray),sign_select{jp});
                jlab=sa.typenames{jpoints(end)};
                %
                ijp=(ip-1)*np+jp;
                labels{iray,jray,ijp}=cat(2,sprintf('%12s',ilab),' ',sprintf('%12s',jlab)); %15Oct24
                str=cat(2,str,labels{iray,jray,ijp},sprintf(': %7.3f  ',angles.ang_degs(iray,jray,ijp)));
            end %jp
        end %ip
        if opts.if_log
            if (jray<iray) | ((jray==iray) & np==2) %display only unique data
                disp(str);
            end
            if (jray==iray)
                disp(' ');
            end
        end
    end %jray
end %iray
opts_used.labels=labels;
return

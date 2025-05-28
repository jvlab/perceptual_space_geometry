function [fnew,s_diag_all,u_full,v_full,s_full,coords_all]=hlid_coords_svd(f,resps_use,maxdim,maxdim_use,if_submean,stims_nonan,stims_nan)
% [fnew,s_diag_all,u_full,v_full,s_full,coords_all]=hlid_coords_svd(f,resps_use,maxdim,maxdim_use,if_submean,stims_nonan,stims_nan)
% is a utility to create coordinates by svd, and metadata fields for a coordinate file
%
% assumes that stimuli without responses have already been eliminated from resps_use
% and installs NaN's in u(stims_nan,:)
%
% f: starting structure of fields to save
% resps_use: response array, each row is a stimulus, each column is a roi or similar
% maxdim: maximum dimension of coords to be created
% maxdim_use: maximum number of dimensions available from data
% if_submean: 1 to subtract the mean, otherwise 0
% stims_nonan: list of stimuli from original set that have responses, defaults to [1:size(resps_use,1)]
% stims_nan: complement of stims_nan in [1:size(resps_use,1)]
%
% fnew: f, with additional fields added
% s_diag_all:  diagonal of S in resps_use=U*S*V', up to maxdim_use
% u_full: full U in resps_use=U*S*V', including any values beyond maxdim_use
% v_full full V in resps_use=U*S*V', including any values beyond maxdim_use
% s_full: full S in resps_use=U*S*V', including any values beyond maxdim_use
% coords_all: all coords written
%
% 19Dec24: fix display bug: s needs to be squared to compute variance
% 06May25: include stims_nan and stims_nonan
%
%   See also:  HLID_RASTIM2COORDS_DEMO, HLID_RASTIM2COORDS_POOL, HLID_CSV2COORDS_DEMO, HLID_ORN_MERGE, HLID_ORN_MERGE2.
%
dim_text='dim'; %leadin for fields of d
if (nargin<=5)
    stims_nonan=[1:size(resps_use,1)];
    stims_nan=setdiff([1:size(resps_use,1)],stims_nonan);
end
%
[u,s,v]=svd(resps_use); %resps_use=u*s*v', with u and v both orthogonal, so u*s=resps_use*v
s_diag_all=diag(s(1:maxdim_use,1:maxdim_use));
s_full=s;
u_full=u;
v_full=v;
var_total=sum(s_diag_all.^2); %^2 added 19Dec24
s=s(1:maxdim_use,1:maxdim_use);
u=u(:,1:maxdim_use);
v=v(:,1:maxdim_use);
if maxdim_use<maxdim %add zero eigenveectors and eigenvalues, should only happen if there are NaN stimuli
    u(stims_nonan,:)=u;
    u(stims_nan,:)=NaN;
    s(:,maxdim_use+1:maxdim)=0; %higher eigenvals and eigenvecs are zero
    s(maxdim_use+1:maxdim,:)=0;
    u(:,maxdim_use+1:maxdim)=0;
    v(:,maxdim_use+1:maxdim)=0;
end
%
s_diag=diag(s);
disp('fraction of variance explained with each component')
disp(sprintf('%6.4f ',s_diag.^2/var_total)) %^2 added 19Dec24
coords_all=u*s;
for idim=1:maxdim
    f.(cat(2,dim_text,sprintf('%1.0f',idim)))=coords_all(:,1:idim);
end
f.coord_opts.method='svd';
f.coord_opts.maxdim=maxdim;
f.coord_opts.if_submean=if_submean;
f.coord_opts.aux.desc='resps, or resps-repmat(mean(resps,1),nstims,1), = u*s*transpose(v) ';
f.coord_opts.aux.u=u;
f.coord_opts.aux.s=s;
f.coord_opts.aux.v=v;
%
fnew=f;
%
return
end



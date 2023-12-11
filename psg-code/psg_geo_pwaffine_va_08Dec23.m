function [d,transform,u]=psg_geo_pwaffine_va(ref,adj,vcut,acut,opts)
%
% [d,transform,u]=psg_geo_pwaffine_va(ref,adj,vcut,acut,opts) finds the best
% piecewise affine model, given a known cutpoint and cut direction
%
% ref: reference dataset, [npts dim_xy], dim_xy >= dim_x
% adj: dataset to adjust, [npts dim_x]
% vcut: unit row vector of length dim_x, orthogonal to cut plane
% acut: cut value
% opts.tol_cut: tolerance for cutpoints (defaults to 10^-7);
%
% d: residuals, normalized for squared dev of ref
% transform: transform structure, see =psg_geo_pwaffine
% u: orthonormal basis, for which vcut is first row
%     coordinates in the new orthonormal basis are given by post-multiplying by u'
%     Note that u depends on vcut but not acut
%
%    See also:  PSG_GEO_PWAFFINE, REGRESS, EXTORTHB.
%
if (nargin<=4)
    opts=struct;
end
opts=filldefault(opts,'tol_cut',10^-7);
dim_x=size(adj,2); %data to adjust
dim_xy=size(ref,2); %refernce data (to fit), already augmented
npts=size(adj,1);
n_pw=2; %number of regions
%
[basis,q]=extorthb(vcut'); %i_dir: vectors are rows but extorthb wants columns
u=q'; %this is the new orthonormal basis, as rows
%coordinates in the new orthonormal basis are given by post-multiplying by u'=q
x_prime=adj*u';
xpa=x_prime(:,1)-acut;
%
%set up regions
%
insides=zeros(npts,n_pw);
insides(:,1)=double(xpa>opts.tol_cut);
insides(:,2)=double(xpa<-opts.tol_cut);
%create augmented regressor matrix
x_prime_aug=zeros(npts,dim_x+2);
for i_pw=1:n_pw
    x_prime_aug(insides(:,i_pw)>0,i_pw)=xpa(insides(:,i_pw)>0);
end
empties=find(all(insides==0,1));
nonempties=setdiff([1:dim_x+2],empties);
%
x_prime_aug(:,3:dim_x+1)=x_prime(:,2:dim_x);
x_prime_aug(:,dim_x+2)=1;
s_aug=zeros(dim_x+2,dim_xy);
s_nz=zeros(length(nonempties),dim_xy);
%do the regression on nonzero regressors
for icol=1:dim_xy
    s_nz(:,icol)=regress(ref(:,icol),x_prime_aug(:,nonempties));
end
s_aug(nonempties,:)=s_nz;
%unpack the results
h=s_aug(dim_x+2,:);
T=zeros(dim_x,dim_xy,n_pw);
c=zeros(n_pw,dim_xy);
for i_pw=1:n_pw
    T(:,:,i_pw)=u'*s_aug([i_pw 3:(dim_x+1)],:); %omit a row
    c(i_pw,:)=h-acut*s_aug(i_pw,:);
end
%compute, display, and analyze residuals
transform=struct;
transform.b=1;
transform.T=T;
transform.c=c;
transform.vcut=vcut;
transform.acut=acut;
adj_model=psg_pwaffine_apply(transform,adj);
d_num=sum(sum((ref-adj_model).^2,1));
d_den=sum(sum((ref-repmat(mean(ref,1),npts,1)).^2,1));
d=d_num/d_den;
return
end

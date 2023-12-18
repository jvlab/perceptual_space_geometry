function [y,sign_vecs,sign_inds]=psg_pwaffine_apply(transform,x)
% [y,sign_vec,sign_inds]=psg_pwaffine_apply(transform,x) applies a piecewise affine transformation
%
%  See psg_piecewise_notes.doc for details on algorithm
%
% x: original coordinates, size=[npts,dim_x]
% transform:
%   transform.b: scalar, equal to 1 (scale absorbed into T)
%   transform.T: stack of matrices, size [dim_x dim_y 2^ncuts]
%      typically ncuts=1, but in general
%      let sign_vecs=sign(x*vcut'-a), of size [npts,ncuts] (with equality going to 1)
%      consider each row of sign_vecs:
%      sign_ind=1       for sign_vec=[+ + .... +]
%      sign_ind=2       for sign_vec=[- + .... +]
%      sign_ind=3       for sign_vec=[+ - .... +]
%      sign_ind=4       for sign_vec=[- - .... +]
%         ....
%      sign_ind=2^ncuts for sign_vec=[- - .... -]
%        use transform.T(:,:,sign_ind for transformation
%   transform.c: stack of offsets, size [ncuts dim_y], use (sign_ind,:)
%   transform.vcut: unit vectors, stack of rows, size [ncuts dim_x], orthog to cut planes
%   transform.acut: vector of length ncuts, the cutpoints
%
% y: transformed coordinates, size=[npts,dim_y]
% sign_vecs: array [dim_x ncuts] of signs (+1,-1), each row as above
% sign_inds: column of length dim_x, sign indices for each row of x
%
% 15Dec23: begin multiple cutpoints
%
%   See also: PSG_GEOMODELS_TEST, PSG_GEO_PWAFFINE.
%
ncuts=size(transform.vcut,1);
cuts=x*transform.vcut'; %column of values of cutpoints
npts=size(x,1);
sign_vecs=sign(cuts-repmat(transform.acut(:)',npts,1));
sign_vecs(sign_vecs==0)=1; %take care of ties
%
sign_inds=1+(1-sign_vecs)/2*(2.^[0:ncuts-1]');%first part maps +1 to 1, -1 to 2
%
n_pw=size(transform.T,3); %number of pieces
ypw=zeros(npts,size(transform.T,2),n_pw);
%ypw(:,:,n_pw) are the alternative values of y in each piece
for ipw=1:n_pw
    ypw(:,:,ipw)=x*transform.T(:,:,ipw)+repmat(transform.c(ipw,:),npts,1);
end
y=zeros(npts,size(transform.T,2));
for ipw=1:n_pw
    y(sign_inds==ipw,:)=ypw(sign_inds==ipw,:,ipw);
end
%%specific to a single cut
%y=ypw(:,:,1);
%y(cuts<transform.acut,:)=ypw(cuts<transform.acut,:,2);
return
end

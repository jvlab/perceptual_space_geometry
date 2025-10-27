%tetra_knit_test:  test whether one can knit together planar datasets to make a tetrahedron
%
%   See also: PSG_CONSENSUS.
%
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
if (if_frozen~=0) 
    rng('default');
    if (if_frozen<0)
        rand(1,abs(if_frozen));
    end
else
    rng('shuffle');
end
%
equilat=[0 0;1 0;0.5 sqrt(3)/2];
z=NaN(4,2,4);
z([1 2 3],:,1)=equilat;
z([1 2 4],:,2)=equilat;
z([1 3 4],:,3)=equilat;
z([2 3 4],:,4)=equilat;
opts_pcon=struct;
opts_pcon.allow_scale=getinp('1 to allow scaling','d',[0 1],0);
opts_pcon.initialize_set='pca';
disp(' distance ranges for various setups')
[con,zn,tx,de,ou]=procrustes_consensus(z,opts_pcon);
disp(' z: two-d setup');
dists=sqrt(cootodsq(con(1:4,:)))+diag(NaN(1,4)); disp([min(dists(:)) max(dists(:))]);
disp(sprintf(' niters: %3.0f',length(de.ts_cum)));
%
%pad to 3d with zeros
%
z3d=z;
z3d(:,3,:)=0;
[con_3d,zn_3d,tx_3d,de_3d,ou_3d]=procrustes_consensus(z3d,opts_pcon);
disp(' z3d: zero-padded to 3d, pc initialization');
dists=sqrt(cootodsq(con_3d(1:4,:)))+diag(NaN(1,4)); disp([min(dists(:)) max(dists(:))]);
disp(sprintf(' niters: %3.0f',length(de_3d.ts_cum)));
%
z3d=z;
z3d(:,3,:)=0;
[con_3dnr,zn_3dnr,tx_3dnr,de_3dnr,ou_3dnr]=procrustes_consensus(z3d,setfield(opts_pcon,'if_initpca_rot',0));
disp(' z3d: zero-padded to 3d, pc initialization without rotation');
dists=sqrt(cootodsq(con_3dnr(1:4,:)))+diag(NaN(1,4)); disp([min(dists(:)) max(dists(:))]);
disp(sprintf(' niters: %3.0f',length(de_3dnr.ts_cum)));
%
opts_pcon.initial_guess=randn(size(z3d,1),size(z3d,2));
[con_3d_ri,zn_3d_ri,tx_3d_ri,de_3d_ri,ou_3d_ri]=procrustes_consensus(z3d,setfield(opts_pcon,'initialize_set',0));
disp(' z3d_ri: zero-padded to 3d, random initzialiation');
dists=sqrt(cootodsq(con_3d_ri(1:4,:)))+diag(NaN(1,4)); disp([min(dists(:)) max(dists(:))]);
disp(sprintf(' niters: %3.0f',length(de_3d_ri.ts_cum)));
%
opts_pcon.initial_guess=randn(size(z3d,1),size(z3d,2));
opts_pcon=rmfield(opts_pcon,'initial_guess');
[con_3d_rg,zn_3d_rg,tx_3d_rg,de_3d_rg,ou_3d_rg]=procrustes_consensus(z3d,setfield(opts_pcon,'initialize_set',0));
disp(' z3d_rg: zero-padded to 3d, random initialiation internally generated');
dists=sqrt(cootodsq(con_3d_rg(1:4,:)))+diag(NaN(1,4)); disp([min(dists(:)) max(dists(:))]);
disp(sprintf(' niters: %3.0f',length(de_3d_rg.ts_cum)));
%
%pad to 3d and permute dimensions
%
z3p=z3d;
z3p(:,[3 1 2],2)=z3d(:,:,2);
z3p(:,[2 3 1],3)=z3d(:,:,3);
z3p(:,[1 3 2],4)=z3d(:,:,4);
[con_3p,zn_3p,tx_3p,de_3p,ou_3p]=procrustes_consensus(z3p,opts_pcon);
disp(' z3p: zero-padded to 3d, dimensions permuted, pc initialization');
dists=sqrt(cootodsq(con_3p(1:4,:)))+diag(NaN(1,4)); disp([min(dists(:)) max(dists(:))]);
disp(sprintf(' niters: %3.0f',length(de_3p.ts_cum)));
%
opts_pcon.initial_guess=randn(size(z3p,1),size(z3p,2));
[con_3p_ri,zn_3p_ri,tx_3p_ri,de_3p_ri,ou_3p_ri]=procrustes_consensus(z3p,setfield(opts_pcon,'initialize_set',0));
disp(' z3p_ri: zero-padded to 3d, dimensions permuted, random initialization');
dists=sqrt(cootodsq(con_3p_ri(1:4,:)))+diag(NaN(1,4)); disp([min(dists(:)) max(dists(:))]);
disp(sprintf(' niters: %3.0f',length(de_3p_ri.ts_cum)));
%
%pad to 8 dimensions, random initialization
%
z8=zeros(4,8,4);
z8(:,[1 2],1)=z(:,[1 2],1);
z8(:,[3 4],2)=z(:,[1 2],2);
z8(:,[5 6],3)=z(:,[1 2],3);
z8(:,[7 8],4)=z(:,[1 2],4);
opts_pcon.initial_guess=randn(size(z8,1),size(z8,2));
[con_8,zn_8,tx_8,de_8,ou_8]=procrustes_consensus(z8,setfield(opts_pcon,'initialize_set',0));
disp(' z8: spread to 8d, random initialization');
dists=sqrt(cootodsq(con_8(1:4,:)))+diag(NaN(1,4)); disp([min(dists(:)) max(dists(:))]);
disp(sprintf(' niters: %3.0f',length(de_8.ts_cum)));
%
%pad to 8 dimensions, add additional points with jitter
%
z8a=[z8;z8+0.01*randn(size(z8));z8+0.01*randn(size(z8))];
[con_8a,zn_8a,tx_8a,de_8a,ou_8a]=procrustes_consensus(z8a,opts_pcon);
disp(' z8a: spread to 8d, random augmentation, pc initialiation');
dists=sqrt(cootodsq(con_8a(1:4,:)))+diag(NaN(1,4)); disp([min(dists(:)) max(dists(:))]);
disp(sprintf(' niters: %3.0f',length(de_8a.ts_cum)));








%psg_geomodels_apply_test: test applying geometric models
%
%   See also: PSG_GEOMODELS_APPLY, PSG_PWAFFINE_APPLY, PERSP_APPLY.
%
if ~exist('ref_dim') ref_dim=4; end
if ~exist('adj_dim') adj_dim=3; end
if ~exist('npts') npts=1000; end
if ~exist('ncuts') ncuts=2; end
if ~exist('tol') tol=10^-7; end
%
string_ans={'OK','NOT OK'};
%
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
%
if (if_frozen~=0)
    rng('default');
    if (if_frozen<0)
        rand(1,abs(if_frozen));
    end
else
    rng('shuffle');
end
%
%create a full model with 2 cutpoints, projective
%
T_stack=randn(adj_dim,ref_dim,2^ncuts); %transformation is from adj_dim to ref_dim
vcut=randn(ncuts,adj_dim);
vcut=vcut./repmat(sqrt(sum(vcut.^2,2)),1,adj_dim);
acut=randn(ncuts,1);
b=randn(1);
c=randn(2^ncuts,ref_dim);
p_list=randn(adj_dim,2^ncuts);
%
xvals=randn(npts,adj_dim);
%
%consider just affine
%
Taffine=struct;
Taffine.T=T_stack(:,:,1);
Taffine.b=1;
Taffine.c=c(1,:);
Taffine.p=zeros(adj_dim,1);
%
Tpwaffine=Taffine;
Tpwaffine.T=T_stack;
Tpwaffine.c=c;
Tpwaffine.vcut=vcut;
Tpwaffine.acut=acut;
Tpwaffine.p=zeros(adj_dim,2^ncuts);
%
Tprojective=Taffine;
Tprojective.p=p_list(:,1);
%
Tpwprojective=Tpwaffine;
Tpwprojective.p=p_list;
%
yvals=struct;
%
%affine
yvals.affine=psg_geomodels_apply('affine',xvals,Taffine);
%
%projective
yvals.projective_trivial=psg_geomodels_apply('projective',xvals,Taffine);
yvals.projective=psg_geomodels_apply('projective',xvals,Tprojective);
yvals.projective_direct=persp_apply(Tprojective.T,Tprojective.c,Tprojective.p,xvals);
%
%piecewise affine
yvals.pwaffine_nocuts=psg_geomodels_apply('pwaffine',xvals,Taffine);
yvals.pwaffine_nocuts_direct=psg_pwaffine_apply(Taffine,xvals);
yvals.pwaffine=psg_geomodels_apply('pwaffine',xvals,Tpwaffine);
[yvals.pwaffine_direct,sign_vecs_pwaffine,sign_inds_pwaffine,yalts_pwaffine]=psg_pwaffine_apply(Tpwaffine,xvals);
%
%piecewise projective
yvals.pwprojective_trivial=psg_geomodels_apply('pwprojective',xvals,Tpwaffine);
yvals.pwprojective_nocuts=psg_geomodels_apply('pwprojective',xvals,Tprojective);
yvals.pwprojective=psg_geomodels_apply('pwprojective',xvals,Tpwprojective);
[yvals.pwprojective_direct,sign_vecs_pwprojective,sign_inds_pwprojective,yalts_pwprojective]=psg_pwprojective_apply(Tpwprojective,xvals);
%
%compare affine with pwaffine, but no cuts, should match
[ifok,matches]=psg_geomodels_apply_util(yvals,'affine','pwaffine_nocuts','same',tol);
%
%compare affine with direct pwaffine no cuts, should match
[ifok,matches]=psg_geomodels_apply_util(yvals,'affine','pwaffine_nocuts_direct','same',tol);
%
%compare affine with each piece of a nontrivial piecewise affine, only first piece should match affine
for ia=1:2^ncuts
    if (ia==1)
        samediff='same';
    else
        samediff='diff';
    end
     [ifok,matches]=psg_geomodels_apply_util(setfield(yvals,'pwaffine_piece',yalts_pwaffine(:,:,ia)),'affine','pwaffine_piece',samediff,tol);
end
%
%compare affine with projective with trivial projection, should match
[ifok,matches]=psg_geomodels_apply_util(yvals,'affine','projective_trivial','same',tol);
%
%compare affine with projective, should not match
[ifok,matches]=psg_geomodels_apply_util(yvals,'affine','projective','diff',tol);
%
%compare affine with pwaffine, should not match unless ncuts=0
if (ncuts==0)
    samediff='same';
else
    samediff='diff';
end
[ifok,matches]=psg_geomodels_apply_util(yvals,'affine','pwaffine',samediff,tol);
%
%compare piecewise affine with piecewise affine direct, should match
[ifok,matches]=psg_geomodels_apply_util(yvals,'pwaffine','pwaffine_direct','same',tol);
%
%compare piecewise affine with piecewise projective with trivial projection, should match
[ifok,matches]=psg_geomodels_apply_util(yvals,'pwaffine','pwprojective_trivial','same',tol);
%compare piecewise affine with piecewise projective, should not match
[ifok,matches]=psg_geomodels_apply_util(yvals,'pwaffine','pwprojective','diff',tol);
%
%compare projective with direct projective
[ifok,matches]=psg_geomodels_apply_util(yvals,'projective','projective_direct','same',tol);
%
%compare projective with pwprojective, no cuts
[ifok,matches]=psg_geomodels_apply_util(yvals,'projective','pwprojective_nocuts','same',tol);
%compare projective with each piece of a nontrivial piecewise projective, only first piece should match affine
for ia=1:2^ncuts
    if (ia==1)
        samediff='same';
    else
        samediff='diff';
    end
     [ifok,matches]=psg_geomodels_apply_util(setfield(yvals,'pwprojective_piece',yalts_pwprojective(:,:,ia)),'projective','pwprojective_piece',samediff,tol);
end
%
%compare projective with pwprojective, should not match unless ncuts=0
if (ncuts==0)
    samediff='same';
else
    samediff='diff';
end
[ifok,matches]=psg_geomodels_apply_util(yvals,'projective','pwprojective',samediff,tol);
%compare piecewise projective with piecewise projective direct, should match
[ifok,matches]=psg_geomodels_apply_util(yvals,'pwprojective','pwprojective_direct','same',tol);

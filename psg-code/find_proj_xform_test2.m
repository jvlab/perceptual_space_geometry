%find_proj_xform_test2: find a projective transformation in n-dimensional space
%
%
% differs from find_proj_xform_test in that 
%  * also tests fmin method, with and without multiple start options
%  * also checks persp_apply
%  * correctly gives discrepancies as rms (no square root applied in proj_xform_test)
%  * different ways to choose values for the data points
%  * 'shrinkage' >1 moves the data outside the viewpoint, fmin fits become unreliable;
%  * 'shrinkage' < 1, fmin fits are more reliable
% 
%  Note that fmin method does not work if the projective component proj_true_hom(1:end-1,end) is significant 
%
% nd dimensions, nd+1 homogeneous dimensions; number of points must be at
% least nd+2 for problem to be determined.
%
% Tests the 'oneshot' method (ref below) and also the 'fmin' method
% Method: Estimating Projective Transformation Matrix (Collineation, Homography)
% Zhengyou Zhang
% November 1993; Updated May 29, 2010
% Microsoft Research Techical Report MSR-TR-2010-63
%
% This uses Method 2 of the above reference, but the data are in rows
%
%   See also: REGRESS, PERSP_XFORM_FIND, FIND_PROJ_XFORM_TEST, PERSP_APPLY.
%
% generate some test data
%
rng('default');
if ~exist('jit') jit=0.01; end; %random jitter
if ~exist('npts') npts=25; end; %number of points
if ~exist('nd') nd=3; end %dimension
if ~exist('shrinkage') shrinkage=1; end %to ensure that the data are inside the viewpoint
%
proj_true_hom=randn(nd+1,nd+1);
%
viewpoint_recip=(proj_true_hom(1:end-1,end)/proj_true_hom(end,end));
viewpoint_mag=1./sqrt(sum(viewpoint_recip.^2));
%put the data points inside a sphere of radius viewpoint_mag*shrinkage
x_gau=randn(npts,nd);
x_gau_shell=x_gau./repmat(sqrt(sum(x_gau.^2,2)),1,nd);
x_rads=rand(npts,1).^(1/nd); %to make points uniformly distributeed inside a unit sphere
x_orig=x_gau_shell.*repmat(x_rads,1,nd);
x_orig_hom=[x_orig*viewpoint_mag*shrinkage ones(npts,1)];
% 
y_orig_hom=x_orig_hom*proj_true_hom;
%
%generate data to test the algorithm:  remove the homogeneous multiplier
%
x_orig=x_orig_hom(:,1:nd)./repmat(x_orig_hom(:,end),1,nd);
y_orig_unjit=y_orig_hom(:,1:nd)./repmat(y_orig_hom(:,end),1,nd);
y_orig=y_orig_unjit+jit*randn(npts,nd);
%
%augment the data
x_aug=[x_orig,ones(npts,1)];
y_aug=[y_orig,ones(npts,1)];
%
%create the "A" matrix 
%
A=zeros((nd+1)*npts,(nd+1)^2-1+npts);
for ipt=1:npts
    Mi=zeros(nd+1,(nd+1)^2);
    for id=1:nd+1
        Mi(id,(nd+1)*(id-1)+[1:(nd+1)])=x_aug(ipt,:);
    end
    Arows=(nd+1)*(ipt-1)+[1:(nd+1)];
    A(Arows,1:(nd+1)^2)=Mi;
    if (ipt>1)
        A(Arows,(nd+1)^2+ipt-1)=y_aug(ipt,:)';
    end  
end
b=zeros((nd+1)*npts,1);
b(1:nd+1)=y_aug(1,:);
%
S=regress(b,A);
mults_found=S((nd+1)^2+1:end); %multipliers for y2, y3, ...
proj_found_hom=reshape(S(1:(nd+1)^2),nd+1,nd+1);
disp('orig homogeneous projection matrix')
disp(proj_true_hom);
disp('found homogeneous projection matrix');
disp(proj_found_hom)
disp('quotient')
disp(proj_found_hom./proj_true_hom);
%
disp('using persp_xform_find with if_cycle=0')
[p,y_fit,ou]=persp_xform_find(x_orig,y_orig);
disp('p')
disp(p);
disp('quotient');
disp(p./proj_true_hom);
disp('rms diff of y_fit and y_orig');
disp(sqrt(mean((y_fit(:)-y_orig(:)).^2)));
%verify transformation with persp_apply
q=p/p(end,end); %normalize
a=q(1:end-1,1:end-1);b=q(end,1:end-1);c=q(1:end-1,end);y_fit_check=persp_apply(a,b,c,x_orig);
disp('rms diff of y_fit and y_orig transformed by persp_apply');
disp(sqrt(mean((y_fit(:)-y_fit_check(:)).^2)));
%
disp('using persp_xform_find with if_cycle=1')
[p_cyc,y_fit_cyc,ou_cyc]=persp_xform_find(x_orig,y_orig,setfield([],'if_cycle',1));
disp('p_cyc')
disp(p_cyc);
disp('quotient');
disp(p_cyc./proj_true_hom);
disp('rms diff of y_fit and y_orig');
disp(sqrt(mean((y_fit_cyc(:)-y_orig(:)).^2)));
disp(sprintf('best point: %4.0f',ou_cyc.oneshot.best_point));
%verify transformation with persp_apply
q=p_cyc/p_cyc(end,end); %normalize
a=q(1:end-1,1:end-1);b=q(end,1:end-1);c=q(1:end-1,end);y_fit_cyc_check=persp_apply(a,b,c,x_orig);
disp('rms diff of y_fit and y_orig transformed by persp_apply');
disp(sqrt(mean((y_fit_cyc(:)-y_fit_cyc_check(:)).^2)));
%
disp('using persp_xform_find, fmin method, initialize at zero')
[p_fmin,y_fit_fmin,ou_fmin]=persp_xform_find(x_orig,y_orig,setfield([],'method','fmin'));
disp('p_fmin')
disp(p_fmin);
disp('quotient');
disp(p_fmin./proj_true_hom);
disp('rms diff of y_fit and y_orig');
disp(sqrt(mean((y_fit_fmin(:)-y_orig(:)).^2)));
%verify transformation with persp_apply
q=p_fmin/p_fmin(end,end); %normalize
a=q(1:end-1,1:end-1);b=q(end,1:end-1);c=q(1:end-1,end);y_fit_fmin_check=persp_apply(a,b,c,x_orig);
disp('rms diff of y_fit and y_orig transformed by persp_apply');
disp(sqrt(mean((y_fit_fmin(:)-y_fit_fmin_check(:)).^2)));
%
disp('using persp_xform_find, fmin method, multiple start points')
[p_fmin_multi,y_fit_fmin_multi,ou_fmin_multi]=persp_xform_find(x_orig,y_orig,setfields([],{'method','if_multi'},{'fmin',1}));
disp('p_fmin_multi')
disp(p_fmin_multi);
disp('quotient');
disp(p_fmin_multi./proj_true_hom);
disp('rms diff of y_fit and y_orig');
disp(sqrt(mean((y_fit_fmin_multi(:)-y_orig(:)).^2)));
%verify transformation with persp_apply
q=p_fmin_multi/p_fmin_multi(end,end); %normalize
a=q(1:end-1,1:end-1);b=q(end,1:end-1);c=q(1:end-1,end);y_fit_fmin_multi_check=persp_apply(a,b,c,x_orig);
disp('rms diff of y_fit and y_orig transformed by persp_apply');
disp(sqrt(mean((y_fit_fmin_multi(:)-y_fit_fmin_check(:)).^2)));
if isfield(ou_fmin_multi.fmin,'x0_best')
    disp('best starting point');
    disp(ou_fmin_multi.fmin.x0_best');
else
    warning('multi option not available in active version of persp_xform_find');
end
%
disp('using persp_xform_find, best method, if_cycle and multiple start points')
[p_best,y_fit_best,ou_best]=persp_xform_find(x_orig,y_orig,setfields([],{'method','if_cycle','if_multi'},{'best',1,1}));
disp('p_best')
disp(p_best);
disp('quotient');
disp(p_best./proj_true_hom);
disp('rms diff of y_fit and y_orig');
disp(sqrt(mean((y_fit_best(:)-y_orig(:)).^2)));
%verify transformation with persp_apply
q=p_best/p_best(end,end); %normalize
a=q(1:end-1,1:end-1);b=q(end,1:end-1);c=q(1:end-1,end);y_fit_best_check=persp_apply(a,b,c,x_orig);
disp('rms diff of y_fit and y_orig transformed by persp_apply');
disp(sqrt(mean((y_fit_best(:)-y_fit_fmin_check(:)).^2)));
%
disp(sprintf('summary of rms diffs, nd=%2.0f, npts=%3.0f, jit=%6.3f, shrinkage=%4.3f',nd,npts,jit,shrinkage));
disp(sprintf('%9.5f: %30s',sqrt(mean((y_fit(:)-y_orig(:)).^2)),'oneshot (default)'));
disp(sprintf('%9.5f: %30s',sqrt(mean((y_fit_cyc(:)-y_orig(:)).^2)),'oneshot, if_cyc=1'));
disp(sprintf('%9.5f: %30s',sqrt(mean((y_fit_fmin(:)-y_orig(:)).^2)),'fmin'));
disp(sprintf('%9.5f: %30s',sqrt(mean((y_fit_fmin_multi(:)-y_orig(:)).^2)),'fmin, if_multi=1'));
disp(sprintf('%9.5f: %30s',sqrt(mean((y_fit_best(:)-y_orig(:)).^2)),'best, if_cyc=1, if_multi=1'));


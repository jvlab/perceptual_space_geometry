%find_proj_xform_test: find a projective transformation in n-dimensional space
%
% nd dimensions, nd+1 homogeneous dimensions; number of points must be at
% least nd+2 for problem to be determined.
%
% Method: Estimating Projective Transformation Matrix (Collineation, Homography)
% Zhengyou Zhang
% November 1993; Updated May 29, 2010
% Microsoft Research Techical Report MSR-TR-2010-63
%
% This uses Method 2 of the above reference, but the data are in rows
%
%   See also: REGRESS, PERSP_XFORM_FIND.
%
% generate some test data
%
rng('default');
if ~exist('jit') jit=0.01; end; %random jitter
if ~exist('npts') npts=25; end; %number of points
if ~exist('nd') nd=3; end %dimension
x_orig_hom=randn(npts,nd+1)-1; % 
proj_true_hom=randn(nd+1,nd+1);
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
disp(mean((y_fit(:)-y_orig(:)).^2));
%
disp('using persp_xform_find with if_cycle=1')
[p_cyc,y_fit_cyc,ou_cyc]=persp_xform_find(x_orig,y_orig,setfield([],'if_cycle',1));
disp('p_cyc')
disp(p_cyc);
disp('quotient');
disp(p_cyc./proj_true_hom);
disp('rms diff of y_fit and y_orig');
disp(mean((y_fit_cyc(:)-y_orig(:)).^2));
disp(sprintf('best point: %4.0f',ou_cyc.oneshot.best_point));

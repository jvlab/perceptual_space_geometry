%psg_procrustes_regr_test: test Procrustes and regression routines
%
% this works but seems to not provide robust estimates to the projective
% transformation: if one artifically makes the adj file equal to a
% projective transformation of ref, then there's a zero residual -- but
% estimating the projective transformation seems not so good.
%
%   See also:  PROCRUSTES, REGRESS, PERSP_XFORM_FIND, PSG_GET_COORDSETS.
%
if ~exist('ref_dim') ref_dim=3; end
if ~exist('if_regr_const') if_regr_const=1; end
if ~exist('if_center') if_center=1; end
if ~exist('nshuff_proc') nshuff_proc=100; end
if ~exist('nshuff_regr') nshuff_regr=100; end
if ~exist('if_cycle') if_cycle=1; end %option for projective fit
if ~exist('opts_read') opts_read=struct; end
if ~exist('opts_rays') opts_rays=struct; end
if ~exist('opts_qpred') opts_qpred=struct; end
if ~exist('opts_persp') opts_persp=struct; end
%
if_builtin=getinp('use built-in datasets','d',[0 1]);
if if_builtin
    if ~exist('ref_file') ref_file='./psg_data/bgca3pt_coords_MC_sess01_10.mat'; end
    if ~exist('adj_file') adj_file='./psg_data/bgca3pt_coords_BL_sess01_10.mat'; end
else
    disp('dataset 1 will be reference, dataset 2 will be adjusted to fit.');
    nsets=2;
    opts_read=filldefault(opts_read,'if_log',1);
    [sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred,nsets);
    ref_file=sets{1}.label_long;
    adj_file=sets{2}.label_long;
end
disp(sprintf('reference dataset: %s',ref_file));
disp(sprintf('dataset to be adjusted: %s',adj_file));
%
ref_dim=getinp('dimension of reference dataset','d',[1 7],ref_dim);
adj_dim=getinp('dimension of dataset to be adjusted (no greater than reference)','d',[1 ref_dim],ref_dim);
if if_builtin
    ref=getfield(load(ref_file),sprintf('dim%1.0f',ref_dim));
    adj=getfield(load(adj_file),sprintf('dim%1.0f',adj_dim));
else
    ref=ds{1}{ref_dim};
    adj=ds{2}{adj_dim};
end
%
if_regr_const=getinp('1 to include constant term in regression','d',[0 1],if_regr_const);
if_center=getinp('1 to center the data','d',[0 1],if_center);
%
tstring=sprintf(' ref %s dim %1.0f adj %s dim %1.0f regr const %1.0f center %1.0f',...
    ref_file,ref_dim,adj_file,adj_dim,if_regr_const,if_center);
%
npts=size(ref,1);
if if_center
    ref=ref-repmat(mean(ref,1),npts,1);
    adj=adj-repmat(mean(adj,1),npts,1);
end
%
[d,adj_proc,transform]=procrustes(ref,adj);
disp('transform');
disp(transform);
adj_check=transform.b*adj*transform.T+transform.c;
%
ref_del=ref-adj_proc; %residuals after Procrustes
%
%recalculate normalized measure of deviation (d) from Procrustes
%
d_den=sum(sum((ref-repmat(mean(ref,1),npts,1)).^2,1));
d_num=sum(sum((ref_del).^2,1));
d_check=d_num/d_den;
disp(sprintf('d (procrustes): %12.7f   d_check: %12.7f   diff: %12.7f',d,d_check,d-d_check))
%
% compare d with shuffle surrogates
%
if nshuff_proc>0
    d_shuff_proc=zeros(1,nshuff_proc);
    for ishuff=1:nshuff_proc
        d_shuff_proc(ishuff)=procrustes(ref,adj(randperm(npts),:));
    end
    disp(sprintf('d is >= d_shuff for %5.0f of %5.0f shuffles, d=%12.7f, min(shuffle)=%12.7f',...
        sum(double(d>=d_shuff_proc)),nshuff_proc,d,min(d_shuff_proc)));
end
%
% regression fit from scratch
%
b_regr=zeros(if_regr_const+adj_dim,ref_dim);
if (if_regr_const==1)
    adj_aug=[ones(npts,1),adj];
else
    adj_aug=adj;
end
for iref_dim=1:ref_dim
    b_regr(:,iref_dim)=regress(ref(:,iref_dim),adj_aug);
end
adj_regr=adj_aug*b_regr; %reconstitute based on regression
%
d_regr_num=sum(sum((ref-adj_regr).^2,1));
d_regr=d_regr_num/d_den;
disp(sprintf('d (regression): %12.7f',d_regr))
%
% compare d with shuffle surrogates in which we scramble the residuals from the Procrustes fits
% to see if the regression is significant
%
if nshuff_regr>0
    resids=ref-adj_regr;
    d_shuff_regr=zeros(1,nshuff_regr);
    for ishuff=1:nshuff_regr
        b_shuf=zeros(if_regr_const+adj_dim,ref_dim);
        ref_shuf=ref-adj_regr+adj_regr(randperm(npts),:); %reference with scrambled resituals
        for iref_dim=1:ref_dim
            b_shuff(:,iref_dim)=regress(ref_shuf(:,iref_dim),adj_aug);
        end
        adj_shuf=adj_aug*b_shuff;
        d_shuff_regr(ishuff)=sum(sum((ref_shuf-adj_shuf).^2,1))/d_den;
    end
    disp(sprintf('d is >= d_shuff for %5.0f of %5.0f shuffles, d=%12.7f, min(shuffle)=%12.7f',...
        sum(double(d_regr>=d_shuff_regr)),nshuff_regr,d_regr,min(d_shuff_regr)));
end
%
%projective fit from scratch
%
[persp,adj_proj,ou_persp]=persp_xform_find(adj,ref,setfield(opts_persp,'if_cycle',if_cycle));
d_proj_num=sum(sum((ref-adj_proj).^2,1));
d_proj=d_proj_num/d_den;
disp(sprintf('d (projective): %12.7f',d_proj))
%
%projective fit to Procrustes residuals -- not "fair" because a sum of a
%linear transformatoin and a projective transformation need NOT be a
%projective transformation
%
[persp_del,adj_del,ou_persp_del]=persp_xform_find(adj,ref_del,setfield(opts_persp,'if_cycle',if_cycle));
d_del_num=sum(sum((ref_del-adj_del).^2,1));
d_del=d_del_num/d_den;
disp(sprintf('d (projective, fitting Procrustes residuals): %12.7f',d_del))
%
% plots
%
if (ref_dim==2) | (ref_dim==3)
    figure;
    set(gcf,'Position',[50 50 1000 800]);
    switch ref_dim
        case 2
            hr=plot(ref(:,1),ref(:,2),'k.','MarkerSize',8);
            hold on;
            hap=plot(adj_proc(:,1),adj_proc(:,2),'r.','MarkerSize',8);
            hac=plot(adj_check(:,1),adj_check(:,2),'gx','MarkerSize',6);
            har=plot(adj_regr(:,1),adj_regr(:,2),'m.','MarkerSize',8);
            for k=1:npts
                plot([ref(k,1) adj_proc(k,1)],[ref(k,2) adj_proc(k,2)],'r');
                plot([ref(k,1) adj_regr(k,1)],[ref(k,2) adj_regr(k,2)],'m');
            end
        case 3
            hr=plot3(ref(:,1),ref(:,2),ref(:,3),'k.','MarkerSize',8);
            hold on;
            hap=plot3(adj_proc(:,1),adj_proc(:,2),adj_proc(:,3),'r.','MarkerSize',8);
            hac=plot3(adj_check(:,1),adj_check(:,2),adj_check(:,3),'gx','MarkerSize',6);
            har=plot3(adj_regr(:,1),adj_regr(:,2),adj_regr(:,3),'m.','MarkerSize',8);
            for k=1:npts
                plot3([ref(k,1) adj_proc(k,1)],[ref(k,2) adj_proc(k,2)],[ref(k,3) adj_proc(k,3)],'r');
                plot3([ref(k,1) adj_regr(k,1)],[ref(k,2) adj_regr(k,2)],[ref(k,3) adj_regr(k,3)],'m');
            end
            view(3);
            zlabel('dim 3');
    end
    xlabel('dim 1');
    ylabel('dim 2');
    axis equal;
    legend([hr;hap;hac;har],{'ref','adjusted (Procrustes)','check (Procrustes)','regression'});
    box on;
    grid on;
    title(sprintf('d (procrustes) %9.6f (regression) %9.6f',d,d_regr));
    %
    axes('Position',[0.01,0.01,0.01,0.01]); %for text
    text(0,0,tstring,'Interpreter','none','FontSize',10);
    axis off;
end

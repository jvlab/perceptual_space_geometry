%psg_pwprojective_test:  test psg_projective_apply, especially w.r.t. for continuity
%
% Uses psg_geo_pwaffine_va to determine a piecewise affine transformation
% that satisfies continuity conditions for T, and then uses 
% psg_geo_pwprojective_pset to set up projective params that satisfy continuity conditions
% 
% See also:  PSG_GEO_PW_VA_TEST, PSG_GEO_PWPROJECTIVE_PSET, PSG_PWPROJECTIVE_APPLY, PSG_GEO_PWAFFINE_VA.
%
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
if ~exist('npts') npts=100; end
if ~exist('ncuts_max') ncuts_max=2; end
%
if ~exist('dx_list') dx_list=[3]; end
if ~exist('dy_list') dy_list=[2 3]; end
%
x_all=randn(npts,max(dx_list));
y_all=randn(npts,max(dy_list));
%
if ~exist('plotvals')
    plotvals=cell(1,2);
    nplotvals=zeros(1,2);
    for ip=1:2
        plotvals{ip}=[-10:1:10]; %range to plot for first and second coord -- can make unequal to debug
        nplotvals(ip)=length(plotvals{ip});
    end
end
if ~exist('proj_safe') proj_safe=3; end %to ensure that projection denominator is small in comparison with coordinates
%
vcut_unnorm=randn(ncuts_max,max(dx_list));
acut_all=randn(1,ncuts_max);
p0_all=randn(max(dx_list),1)/max([max(abs(plotvals{1})),max(abs(plotvals{2}))])/proj_safe;
w_all=randn(1,ncuts_max);
w_all=zeros(1,ncuts_max);
for dx_ptr=1:length(dx_list)
    dim_x=dx_list(dx_ptr);
    for dy_ptr=1:length(dy_list)
        dim_y=dy_list(dy_ptr);
        for ncuts=1:ncuts_max
            if (ncuts<=dim_x)
                labstring=sprintf('dim_x %2.0f dim_y %2.0f ncuts %2.0f',dim_x,dim_y,ncuts);
                disp(labstring);
                vcut=vcut_unnorm(1:ncuts,1:dim_x);
                vcut=vcut./repmat(sqrt(sum(vcut.^2,2)),[1 dim_x]);
                acut=acut_all(1,1:ncuts);
                %find a piecewise affine transforms that are continuous at boundary
                x=x_all(:,1:dim_x);
                y=y_all(:,1:dim_y);
                [d_pwaffine,transform_pwaffine]=psg_geo_pwaffine_va(y,x,vcut,acut);
                %
                p0=p0_all(1:dim_x);
                w=w_all(1,ncuts);
                %
                p_list=psg_geo_pwprojective_pset(vcut,p0,w);
                transform_pwprojective=setfield(transform_pwaffine,'p',p_list);
                %
                %now examine in each coordinate plane
                %
                dimpairs=nchoosek([1:dim_x],2);
                npairs=size(dimpairs,1);
                figure;
                set(gcf,'Position',[100 100 1400 800]);
                set(gcf,'NumberTitle','off');
                set(gcf,'Name',labstring);
                ncols=max(max(dy_list),3);
                nrows=max(npairs,3);
                for idp=1:npairs
                    x_plot=zeros(prod(nplotvals),dim_x);
                    [xm1,xm2]=meshgrid(plotvals{1},plotvals{2});
                    x_plot(:,dimpairs(idp,1))=xm1(:);
                    x_plot(:,dimpairs(idp,2))=xm2(:);
                    [y_plot,sign_vecs,sign_inds,ypw]=psg_pwprojective_apply(transform_pwprojective,x_plot);
                    for iyd=1:dim_y
                        subplot(nrows,ncols,iyd+(idp-1)*ncols);
                        surf(xm1,xm2,reshape(y_plot(:,iyd),nplotvals(2),nplotvals(1))); %matlab surf plots as transpose
                        set(gca,'ZLim',max(abs(y_plot(:)))*[-1 1]);
                        xlabel(sprintf('x_%1.0f',dimpairs(idp,1)));
                        ylabel(sprintf('x_%1.0f',dimpairs(idp,2)));
                        zlabel(sprintf('y_%1.0f',iyd));
                    end
                end %x dimension pair
                axes('Position',[0.01,0.05,0.01,0.01]); %for text
                text(0,0,labstring,'Interpreter','none');
                axis off;
                %
            end %ncuts ok
        end %ncuts
    end %dy_ptr
end %dx_ptr

%psg_geo_pw_va_test: test psg_geo_pw[affine|projective]_va
%
% This tests function [d,transform,u,opts_used]=psg_geo_pwaffine_va(y,x,vcut,acut,opts) 
% and [d,transform,u,opts_used]=psg_geo_pprojective_va(y,x,vcut,acut,pstruct,opts) 
% Differs from psg_geo_pwaffine_va_test in that this also checks that the transform found is correct
%   To test this, random data (y_all) is fit, and then, using the fit, new data
%   are created from x_all, and then refit.
%   This is intended to be applicable to testing piecewise projective fitting.
%
% from psg_geo_pwaffine_va:
% u: basis used for analysis. vcut is first ncuts rows; remaining rows are orthogonal to vcut
%    * The coordinates in the analysis basis are given by post-multiplying x by uinv.
%    * The first ncuts columns of inv(u) are rows of vcut.
%    * Note that u depends on vcut but not acut
%    *  If if_orth=1 or ncuts>1, the last (dim_x-ncuts) rows of u are orthonormal, and these are orthogonal to the first ncuts rows.
%    The first ncuts rows of u may not be orthogonal to each other, as they are unit vectors orthogonal to the cutplanes.
%    *  If if_orth=0 and ncuts=1, the last (dim_x-1) rows of u are orthogonal to the first row, but may not be orthogonal to each other.
%
%   See also:  PSG_GEO_PWAFFINE_VA, PSG_GEO_PWPROJECTIVE_VA, PSG_GEO_PWAFFINE_VA_TEST, PSG_GEOMODELS_APPLY.
%
%
%  11Jul24: add tests of psg_geo_pwprojectiv_va
%  15Jul24: add tests that the first ncuts rows of u and vcut have the same span
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
if ~exist('d_max') d_max=10; end
if ~exist('npts') npts=100; end
if ~exist('ncuts_max') ncuts_max=3; end
if ~exist('opts') opts=struct; end
if ~exist('tol') tol=10^-10; end
%
x_all=randn(npts,d_max);
y_all=randn(npts,d_max);
vcut_unnorm=randn(ncuts_max,d_max);
acut_all=randn(1,ncuts_max);
%
if ~exist('dx_list') dx_list=[1 2 3 6]; end
if ~exist('dy_list') dy_list=[1 2 3 6]; end
%
plist_lib=randn(max(dx_list),2^ncuts_max);
%
okstring={'BAD','OK-Y','OK-N'};
for dx_ptr=1:length(dx_list)
    dim_x=dx_list(dx_ptr);
    for dy_ptr=1:length(dy_list)
        dim_y=dy_list(dy_ptr);
        disp(' ');
        if_fitok_pwaffine=1; %reset for series of cuts
        if_fitok_pwproj0=1; %all proj params=0 (trivial)
        if_fitok_all=1; %specify all proj params
        if_fitok_same=1; %all proj params the same
        for ncuts=1:ncuts_max
            if (ncuts<=dim_x)
                vcut=vcut_unnorm(1:ncuts,1:dim_x);
                vcut=vcut./repmat(sqrt(sum(vcut.^2,2)),[1 dim_x]);
                acut=acut_all(1,1:ncuts);
                x=x_all(:,1:dim_x);
                y=y_all(:,1:dim_y);
                for if_orth=0:1
                    if (ncuts==1) | (if_orth==1)
                        [d,transform,u,ou]=psg_geo_pwaffine_va(y,x,vcut,acut,setfield(opts,'if_orth',if_orth));
                        transform_fields=fieldnames(transform);
                        y_fit=psg_geomodels_apply('pwaffine',x,transform);
                        [d_refit,transform_refit]=psg_geo_pwaffine_va(y_fit,x,vcut,acut,setfield(opts,'if_orth',if_orth));
                        %compare the fitted pwaffine transform with actual
                        devs=struct;
                        for ifn=1:length(transform_fields)
                            fn=transform_fields{ifn};
                            devs.(fn)=max(abs(transform.(fn)(:)-transform_refit.(fn)(:)));
                            if devs.(fn)>=tol
                                if_fitok_pwaffine=0;
                            end
                        end
                        devs.d=d_refit;
                        if devs.d>=tol
                            if_fitok_pwaffine=0;
                        end
                        if (if_fitok_pwaffine==0)
                            disp(devs);
                        end
                        %compare fitted pwproj assuming trivial projection with actual
                        pstruct=struct;
                        pstruct.mode='zero';
                        [d_pwproj0,transform_pwproj0]=psg_geo_pwprojective_va(y_fit,x,vcut,acut,pstruct,setfield(opts,'if_orth',if_orth));
                        devs=struct;
                        for ifn=1:length(transform_fields)
                            fn=transform_fields{ifn};
                            devs.(fn)=max(abs(transform.(fn)(:)-transform_pwproj0.(fn)(:)));
                            if devs.(fn)>=tol
                                if_fitok_pwproj0=0;
                            end
                        end
                        devs.d_pwproj0=d_pwproj0;
                        if devs.d_pwproj0>=tol
                            if_fitok_pwproj0=0;
                        end
                        if (if_fitok_pwproj0==0)
                            disp(devs);
                        end
                        %
                        %compare fitted pwproj with a full list of projection params
                        %
                        plist_all=plist_lib(1:dim_x,1:2^ncuts);
                        transform_all=transform;
                        transform_all.p=plist_all;
                        y_fit_all=psg_geomodels_apply('pwprojective',x,transform_all);
                        pstruct_all=struct;
                        pstruct_all.mode='all';
                        pstruct_all.vals=plist_all;
                        [d_all,transform_all_refit]=psg_geo_pwprojective_va(y_fit_all,x,vcut,acut,pstruct_all,setfield(opts,'if_orth',if_orth));
                        transform_fields_proj=transform_fields;
                        transform_fields_proj{end+1}='p';
                        devs=struct;
                        for ifn=1:length(transform_fields_proj)
                            fn=transform_fields_proj{ifn};
                            devs.(fn)=max(abs(transform_all.(fn)(:)-transform_all_refit.(fn)(:)));
                            if devs.(fn)>=tol
                                if_fitok_all=0;
                            end
                        end
                        devs.d_all=d_all;
                        if devs.d_all>=tol
                            if_fitok_all=0;
                        end
                        if (if_fitok_all==0)
                            disp(devs);
                        end
                        %compare fitted pwproj whenprojection params are all the same
                        %
                        plist_same=repmat(plist_lib(1:dim_x,1),[1 2^ncuts]);
                        transform_same=transform;
                        transform_same.p=plist_same;
                        y_fit_same=psg_geomodels_apply('pwprojective',x,transform_same);
                        pstruct_same=struct;
                        pstruct_same.mode='same';
                        pstruct_same.vals=plist_same(:,1); 
                        [d_same,transform_same_refit]=psg_geo_pwprojective_va(y_fit_same,x,vcut,acut,pstruct_same,setfield(opts,'if_orth',if_orth));
                        transform_fields_proj=transform_fields;
                        transform_fields_proj{end+1}='p';
                        devs=struct;
                        for ifn=1:length(transform_fields_proj)
                            fn=transform_fields_proj{ifn};
                            devs.(fn)=max(abs(transform_same.(fn)(:)-transform_same_refit.(fn)(:)));
                            if devs.(fn)>=tol
                                if_fitok_same=0;
                            end
                        end
                        devs.d_same=d_same;
                        if devs.d_same>=tol
                            if_fitok_same=0;
                        end
                        if (if_fitok_same==0)
                            disp(devs);
                        end
                        %
                        %now check u
                        %
                        %are the first ncut rows of u within the span of vcut?
                        vproj=vcut'*(inv(vcut*vcut'))*vcut;
                        spandiff=u(1:ncuts,:)-u(1:ncuts,:)*vproj;
                        spanv_check=all(abs(spandiff(:))<tol);
                        spanv_string=okstring{1+spanv_check};
                        %
                        %is vcut within the span of the first ncut rows of u?
                        utop=u(1:ncuts,:);
                        utop_proj=utop'*inv(utop*utop')*utop;
                        spandiff=vcut-vcut*vproj;
                        spanu_check=all(abs(spandiff(:))<tol);
                        spanu_string=okstring{1+spanu_check};
                        %
                        uinv=inv(u);
                        %do the first ncuts columns of uinv match the rows of vcut?
                        %should always be yes
                        col_check=all(all(abs(uinv(:,[1:ncuts])-vcut')<tol));
                        col_string=okstring{1+col_check};
                        %
                        uup=u*u';
                        botrows=[ncuts+1:dim_x];
                        toprows=[1:ncuts];                       
                        %are the last (dim_x-ncuts) rows of u orthogonal?
                        %should be yes unless if_orth=0
                        if dim_x>ncuts
                            lastrow_check=all(all(abs(uup(botrows,botrows)-eye(length(botrows)))<tol));
                            if (if_orth)
                                lastrow_string=okstring{1+lastrow_check};
                            else %need not be orthogonal
                                lastrow_string=okstring{3-2*lastrow_check};
                            end
                        else
                            lastrow_string='N/A';
                        end
                        %are the last (dim_x-ncuts) rows of u orthogonal to the first ncuts rows?
                        %should always be yes
                        if dim_x>ncuts
                            mix_check=all(all(abs(uup(toprows,botrows))<tol));
                            mix_string=okstring{1+mix_check};
                        else
                            mix_string='N/A';
                        end
                        %are the first ncuts rows of u orthogonal to each other? 
                        %should be no (unless there is only one row)
                        if ncuts>1
                            first_check=all(all(abs(uup(toprows,toprows)-eye(ncuts))<tol));
                            first_string=okstring{3-2*first_check};
                        else
                            first_string='N/A';
                        end
                        %
                        disp(sprintf(' dim_x %2.0f dim_y %2.0f ncuts %2.0f if_orth %1.0f d %12.7f span_check [v, utop] [%4s %4s] col_check %4s lastrow_check %4s mix_check %4s first_check %4s',...
                            dim_x,dim_y,ncuts,if_orth,d,spanv_string,spanu_string,col_string,lastrow_string,mix_string,first_string));
                    end
                end
            end %ncuts large enough?
        end %ncuts
        disp(sprintf('if_fitok_pwaffine: %2.0f',if_fitok_pwaffine));
        disp(sprintf('if_fitok_pwproj0: %2.0f',if_fitok_pwproj0));
        disp(sprintf('if_fitok_all: %2.0f',if_fitok_all));
        disp(sprintf('if_fitok_same: %2.0f',if_fitok_same));
    end %dim_y
end %dim_x

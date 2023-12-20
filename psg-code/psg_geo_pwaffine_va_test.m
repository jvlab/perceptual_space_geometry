%psg_geo_pwaffine_va_test: test psg_geo_pwaffine_va
%
%   See also:  PSG_GEO_PWAFFINE_VA.
%
% This tests function [d,transform,u,opts_used]=psg_geo_pwaffine_va(y,x,vcut,acut,opts) and especially u
% from psg_geo_pwaffine_va:
% u: basis used for analysis. vcut is first ncuts rows; remaining rows are orthogonal to vcut
%    * The coordinates in the analysis basis are given by post-multiplying x by uinv.
%    * The first ncuts columns of inv(u) are the rows of vcut.
%    * Note that u depends on vcut but not acut
%    *  If if_orth=1 or ncuts>1, the last (dim_x-ncuts) rows of u are orthonormal, and these are orthogonal to the first ncuts rows.
%    The first ncuts rows of u may not be orthogonal to each other, as they are unit vectors orthogonal to the cutplanes.
%    *  If if_orth=0 and ncuts=1, the last (dim_x-1) rows of u are orthogonal to the first row, but may not be orthogonal to each other.
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
okstring={'BAD','OK-Y','OK-N'};
for dx_ptr=1:length(dx_list)
    dim_x=dx_list(dx_ptr);
    for dy_ptr=1:length(dy_list)
        dim_y=dy_list(dy_ptr);
        disp(' ');
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
                        disp(sprintf(' dim_x %2.0f dim_y %2.0f ncuts %2.0f if_orth %1.0f d %12.7f col_check %4s lastrow_check %4s mix_check %4s first_check %4s',...
                            dim_x,dim_y,ncuts,if_orth,d,col_string,lastrow_string,mix_string,first_string));
                    end
                end
            end %ncuts large enough?
        end %ncuts
    end %dim_y
end %dim_x

function [ha,ou]=btc_soid_plotqf(qf,opts)
% [ha,ou]=btc_soid_plotqf(qf,opts) plots the quadratic form fitted to a btc eperiment
%
% qf:  the quadratic form.  third dimension is > length 1, then
%    qf(:,:,1) is plotted as a point
%    qf(:,:,2) is plotted as a lower CL
%    qf(:,:,3) is plotted as a upper CL
%    qf(:,:,4) is plotted as an x (indicating median)
%  
% opts: plot options, see below
%    opts.coords: coordinate labels (a string of chars)
%    view: if empty, will yield standard 3-d view, otherwize [az el]
%
%  See also:  BTC_SOID_DEMO, BTC_PAIRSNEEDED, BTC_EDIRS, BTC_DEFINE.
%
% ha: pointer to the axis
% ou: options used
ncoords=size(qf,1);
if (nargin<=1) opts=[]; end
opts=filldefault(opts,'fontsize',8);
opts=filldefault(opts,'zero_color','k');
opts=filldefault(opts,'pos_color','r');
opts=filldefault(opts,'neg_color','b');
opts=filldefault(opts,'view',[]);
opts=filldefault(opts,'hashsize',0.05);
if (ncoords==3)
    opts=filldefault(opts,'cords','atg');
end
if (ncoords==10)
    opts=filldefault(opts,'cords','gbcdetuvwa');
end
ou=opts;
%
for ix=1:ncoords
    for iy=1:ncoords
        ixm=ix-0.5;
        iym=iy-0.5;
        btc_soid_plotqf3(ixm,iym,qf(ix,iy,1),'.',opts); hold on;
        btc_soid_plotqf3(ixm*[1 1],iym*[1 1],[0 qf(ix,iy,1)],'-',opts); hold on;
        if (size(qf,3)>=3) %error bar indicators
            iebc=[2:3];
            for ieb=iebc
                btc_soid_plotqf3(ixm+opts.hashsize*[-1 1],iym*[1 1],qf(ix,iy,ieb)*[1 1],'-',opts); hold on;
                btc_soid_plotqf3(ixm*[1 1],iym+opts.hashsize*[-1 1],qf(ix,iy,ieb)*[1 1],'-',opts); hold on;
            end
            btc_soid_plotqf3(ixm*[1 1],iym*[1 1],squeeze(qf(ix,iy,iebc)),':',opts); hold on;
        end
        if (size(qf,3)>=4)
            %make a box
            btc_soid_plotqf3(ixm+opts.hashsize*[-1 1 1 -1 -1],iym+opts.hashsize*[-1 -1 1 1 -1],repmat(qf(ix,iy,4),1,5),'-',opts);
            hold on;
        end
    end
end
for ixy=1:ncoords
    hpl=plot3((ixy-0.5)*[1 1],[0 ncoords],[0 0],'k--');
    set(hpl,'Color',[0.5 0.5 0.5]);
    hpl=plot3([0 ncoords],(ixy-0.5)*[1 1],[0 0],'k--');
    set(hpl,'Color',[0.5 0.5 0.5]);
end
set(gca,'XLim',[0 ncoords]);
set(gca,'XTick',[1:ncoords]-0.5);
set(gca,'XTickLabel',opts.coords');
set(gca,'YLim',[0 ncoords]);
set(gca,'YTick',[1:ncoords]-0.5);
set(gca,'YTickLabel',opts.coords');
if (isempty(opts.view))
    view(3);
else
    set(gca,'View',opts.view);
end
ha=gca;
return

function btc_soid_plotqf3(x,y,z,lt,opts)
%plot a point or a line, with color determined by z-value: opts.pos_color,opts.neg_color,opts.zero_color
%assumes that nonneg and nonpos values are consecutive subsets
if max(z(:))==min(z(:));
    c=opts.zero_color;
    if max(z(:))>0;c=opts.pos_color; end
    if max(z(:))<0;c=opts.neg_color; end
    plot3(x,y,z,cat(2,c,lt)); hold on;
else
    %if some are nonneg some are nonpos, assume that values are consecutive, and plot each subset
    kp=find(z(:)>=0);
    plot3(x(kp),y(kp),z(kp),cat(2,opts.pos_color,lt)); hold on;
    kn=find(z(:)<=0);
    plot3(x(kn),y(kn),z(kn),cat(2,opts.neg_color,lt)); hold on;
end
return

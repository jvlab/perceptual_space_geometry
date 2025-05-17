%btc_meabb_asymp: show asymptotics for btc_meabb, the induced colors_aux
% maximum-entropy fourth-order correlation from second-order correls
%
%   See also:  BTc_MEABB.
%
if ~exist('bclist') bclist=[-.5:.01:.5]; end
n=length(bclist);
a=zeros(n,n);
a_quadratic=zeros(n,n);
a_quartic=zeros(n,n);
for ib=1:n
    for ic=1:n
        b=bclist(ib);
        c=bclist(ic);
        a(ib,ic)=btc_meabb(b,c);
        a_quadratic(ib,ic)=b^2+c^2;
        a_quartic(ib,ic)=b^2*c^2;
    end  
end
figure;
set(gcf,'Position',[100 100 1400 850]);
%
subplot(2,3,1);
surf(bclist,bclist,a);set(gca,'ZLim',[0 .5]);title('a');xlabel('b'),ylabel('c');
%
subplot(2,3,2);
surf(bclist,bclist,a_quadratic);set(gca,'ZLim',[0 .5]);title('quadratic approx');xlabel('b'),ylabel('c');
%
subplot(2,3,3);
surf(bclist,bclist,a_quadratic+a_quartic);set(gca,'ZLim',[0 .5]);title('quartic approx');xlabel('b'),ylabel('c');
%
subplot(2,3,5);
surf(bclist,bclist,a-a_quadratic);set(gca,'ZLim',[-0.025 .025]);title('a-quadratic approx');xlabel('b'),ylabel('c');
%
subplot(2,3,6);
surf(bclist,bclist,a-a_quadratic-a_quartic);set(gca,'ZLim',[-0.025 .025]);title('a-quartic approx');xlabel('b'),ylabel('c');




function y=lognz(x)
%
% y=lognz(x) produces log(x) where x>0, and 0 where x<0
% useful for calculating entropies
%
xx=reshape(x,[1,prod(size(x))]);
y=zeros(size(xx));
y(find(xx>0))=log(xx(find(xx>0)));
y=reshape(y,size(x));
return


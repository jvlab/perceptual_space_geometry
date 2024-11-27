function a=multi_shuff_enum_novec(counts)
% a=multi_shuff_enum_novec(counts) creates an array whose rows are all the ways
% of ordering counts(1) 1's, counts(2) 2's, etc.
% 
% counts may only have non-negative integer entries
%
% faster vectorized version: multi_shuff_enum.m
% 
% dim(a,1) is sum(counts)!/prod(counts(k)!), a multinomial coefficient
% 
% tic;q=multi_shuff_enum_novec([2 2 2 2 2 2]);toc,size(q),factorial(12)/2^6,
% Elapsed time is 8.486275 seconds.
%      7484400          12
%      7484400
% tic;q=multi_shuff_enum_novec([2 2 2 2 2]);toc,size(q),factorial(10)/2^5,
% Elapsed time is 0.177920 seconds.
%       113400          10
%       113400
% tic;q=multi_shuff_enum_novec([1 1 1 1 1 1 1]);toc,size(q),factorial(7)
% Elapsed time is 0.011072 seconds.
%         5040           7
%         5040
% tic;q=multi_shuff_enum_novec([3 3 3 3]);toc,size(q),factorial(12)/6^4
% Elapsed time is 0.419715 seconds.
%       369600          12
%       369600
% 
%  See also:  NCHOOSEK.
%
if length(counts)==1
    a=ones(1,counts(1));
    return
else
    c=counts;
    lc=length(c);
    nck=nchoosek([1:sum(c)],sum(c(1:end-1)));
    a=repmat(lc,size(nck,1),sum(c));
    for r=1:size(nck,1)
        a(r,nck(r,:))=lc-1;
    end
    %do a recursion if length(c)>=3:
    if (lc>=3)
        b=multi_shuff_enum_novec(counts(1:end-1));
        %duplicate each row of a by the number of rows in b
        %replace every occurrence of lc-1 in this, by elements in b
 %      ar=reprow(a,size(b,1));
        ar=reshape(repmat(reshape(a,prod(size(a)),1),1,size(b,1))',size(a).*[size(b,1) 1]);
        br=repmat(b,size(a,1),1);
        for k=1:size(ar,1)
            ar(k,find(ar(k,:)==lc-1))=br(k,:);
        end
        a=ar;
    end
end
return

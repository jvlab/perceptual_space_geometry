function a=multi_shuff_enum(counts,opts)
% a=multi_shuff_enum(counts) creates an array whose rows are all the ways
% of ordering counts(1) 1's, counts(2) 2's, etc.
% 
% counts may only have non-negative integer entries
%
% for unvectorized version see multi_shuff_enum_novec.m
% 
% dim(a,1) is sum(counts)!/prod(counts(k)!), a multinomial coefficient
% 
% tic;q=multi_shuff_enum([2 2 2 2 2 2]);toc,size(q),factorial(12)/2^6,
% Elapsed time is 1.530700 seconds.
%      7484400          12
%      7484400
% tic;qnv=multi_shuff_enum_novec([3 3 3 3]);toc,size(q),factorial(12)/6^4
% Elapsed time is 0.605726 seconds.
%       369600          12
%       369600
% max(abs(q(:)-qnv(:)))
%    0
% tic;q=multi_shuff_enum([2 2 2 2 2]);toc,size(q),factorial(10)/2^5,
% Elapsed time is 0.099007 seconds.
%       113400          10
%       113400
% tic;q=multi_shuff_enum([1 1 1 1 1 1 1]);toc,size(q),factorial(7)
% Elapsed time is 0.004967 seconds.
%         5040           7
%         5040
% tic;q=multi_shuff_enum([3 3 3 3]);toc,size(q),factorial(12)/6^4
% Elapsed time is 0.082875 seconds.
%       369600          12
%       369600
% 
%  See also:  NCHOOSEK.
%
if nargin<2
    opts=struct;
end
opts=filldefault(opts,'if_log',0);
opts=filldefault(opts,'level',0);
if (opts.if_log)
    disp(sprintf('multi_shuff_enum: level %1.0f, counts: %s',opts.level,sprintf('%2.0f ',counts)));
end
if length(counts)==1
    a=ones(1,counts(1));
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
        b=multi_shuff_enum(counts(1:end-1),setfield(opts,'level',opts.level+1));
        %duplicate each row of a by the number of rows in b
        %replace every occurrence of lc-1 in this, by elements in b
 %      ar=reprow(a,size(b,1))'; %working in transpose so that (:) goes along rows of original a
        ar=(reshape(repmat(reshape(a,prod(size(a)),1),1,size(b,1))',size(a).*[size(b,1) 1]))';
        rows=size(ar,2); %transpose
        br=(repmat(b,size(a,1),1))';
        ar(ar==(lc-1))=br(:);
        a=reshape(ar,[sum(c) rows])'; %transpose back
    end    
end
if opts.if_log
    disp(sprintf('size(a): %10.0f %2.0f, values: %1.0f to %1.0f',size(a),min(a(:)),max(a(:))));
end
return


% test routine for loglik_beta, demonstrating balance of continuous and
% discrete normalization
%
%   See also:  LOGLIK_BETA.
%
if ~exist('n') n=6; end
if ~exist('alist') alist=[1 4 16]; end %beta-param
lls=zeros(n+1,2,length(alist));
for aptr=1:length(alist)
    a=alist(aptr);
    disp(sprintf('beta parameter: %4.0f',a));
    for k=0:n
        lls(k+1,1,aptr)=loglik_beta([a a],[k n],setfields([],{'qvec','hvec'},{0.5,1})); %discrete only
        lls(k+1,2,aptr)=loglik_beta([a a],[k n],setfields([],{'qvec','hvec'},{0.5,0})); %continuous only
        disp(sprintf(' [k,n]= [%3.0f %3.0f], log likelihoods: discrete %7.5f continuous %7.5f',...
            k,n,lls(k+1,:,aptr)));
    end
end
    
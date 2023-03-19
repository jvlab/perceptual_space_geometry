function excepts=loglik_beta_except(a,b,p)
% excepts=loglik_beta_except(a,b,p) returns one or more exception flags for
% evaluating the beta-distribution values 
% prod((p^(a-1)(1-p)^(b-1)) for a vector of probabilities
%
% a,b: parameters of the beta distribution, must be >0
% p: vector of probabilities,all in [0,1]
%
% excepts: cell array, empty if no exceptions,
%  contains 'big' if a<1 with any p=0 or b<1 with any p=1
%  contains 'small' if a>1 with any p=0 or b>1 with any p=1
%
%   See also:  LOGLIK_BETA.
%
excepts=cell(0);
if (any(p==0) & a<1) | (any(p==1) & b<1)
    excepts{end+1}='big';
end
if (any(p==0) & a>1) | (any(p==1) & b>1)
    excepts{end+1}='small';
end
return

function [cvecs_new,alpha_range]=btc_split_alpha(cvecs,dict)
% [cvecs,alpha_range]=btc_split_alpha(cvecs,dict) tries to adjust the alpha parameter 
% of a pair to be mixed to keep the probabilities as non-negative as possible
% that can be mixed to crate a target symmetric texture
%
% cvecs:  a set of row vectors (typically size=[2 10] of the coordinates of the component textures
% dict:  output of btc_define (created if not supplied)
%
% opts: opts.nsteps: number of steps for a gamma search (defaults to 9, relevant to strat 2 and 4)
%       opts.ifshow: show intermediate calcs (defaults to 0)
%
% cvecs_new: two new row vectors with same mean as cvecs, for which alpha has been
%   adjusted to be in the middle of the range in which all probabilities are non-negative
%
% alpha_range(ivec,ilohi): size is [2 2], the low and high permissible values of alpha for each
%       component (NaN if no value of alpha leads to nonzero probabilities)
%
if (nargin<1)
    dict=[];
end
if isempty(dict)
    dict=btc_define;
end
%
ia=strmatch('alpha',dict.name_order_aug,'exact'); %position of alpha
tol=10^-6;
%
%assume no values with non-negative probabilities
%
alpha_range=repmat(NaN,[2 2]);
cvecs_new=cvecs;
%
alpha_tot=sum(cvecs(:,ia)); %assume alpha is in final position
%
% exploit that the p2x2's depend linearly on alpha.  So we can determine them as a function
% of alpha from two test cases
p=zeros(2,2,2,2,2,2); %first 4 indices are p2x2.  Fifth index is cvec 1 or 2. 
% Sixth dim is whether alpha=0 or 1 in cvec 1
for alpha_val=0:1
    cvecs_test=cvecs;
    cvecs_test(1,ia)=alpha_val;
    cvecs_test(2,ia)=alpha_tot-alpha_val;
    for ispec=1:2
        corrs=btc_vec2corrs(cvecs_test(ispec,:),dict);
        p(:,:,:,:,ispec,alpha_val+1)=getp2x2_corrs(corrs);
    end
end
p=reshape(p,[32 2]);
intercepts=p(:,1);
slopes=p(:,2)-p(:,1);
a01vals=[-intercepts./slopes,(1-intercepts)./slopes]; %these give the ranges of alpha in the FIRST cvec
% for which each probability (16 for the first component, 16 for the second) are all non-negagtive)
aminmaxvals=[min(a01vals,[],2),max(a01vals,[],2)]; %order the intervals so the low end is first
arange=[max(aminmaxvals(:,1)),min(aminmaxvals(:,2))]; %find the intersection of the intervals
if arange(1)<=arange(2)
    alpha_range(1,:)=arange;
    alpha_range(2,:)=alpha_tot-arange;
    cvecs_new(:,ia)=mean(alpha_range,2);
end
return

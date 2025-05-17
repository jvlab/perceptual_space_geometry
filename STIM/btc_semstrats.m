function [cvecs,vec_made,vec_symok,split_strat,aux,opts_used]=btc_semstrats(vec_targ,dict,opts) 
% [cvecs,vec_made,vec_symok,split_strat,aux,opts_used]=btc_semstrat(vec_targ,dict,opts) determines btc coordinates
% of component textures that can be mixed to crate a target "semantic" texture, with parameters symmetric w.r.t. side-to-side flip
%
%  see btc_semanticdirs_notes.docx for further info
%
% vec_targ: target vector.  gamma should be 0, beta-diags should be equal, t=u, v=w
% dict:  output of btc_define (created if not supplied)
% opts: opts.ifshow: show intermediate calcs (defaults to 0)
%
% cvecs:  a set of row vectors (typically size=[2 10] of the coordinates of the component textures
%         the average of cvecs should correspond to val_targ
% vec_made: the coordinate vector made
% vec_symok: val_targ with symmetry conditions enforced
% split_strat:
%   1->equal weights, gamma=0, non-Pickard thetas set to zero
%   possibly for future:
%    allow for unequal gammas in splits
%    allow for nonzero gamma in val_targ
%    determine maximum-entropy alpha, and range of feasible alphas
% aux: auxiliary values
%   aux.wts: vector of weights for the components (typically [0.5 0.5])
%   
%
%    Note wts*cvecs=vec_made
% opts_used: options used
%
%   See also:  BTC_MAKESEM_DEMO, BTC_SYMSTRAT, BTC_SYMSTRATS,BTC_DEFINE,BTC_VEC2CORRS,GETCORRS_P2X2.
%
if (nargin<2)
    dict=[];
end
if isempty(dict)
    dict=btc_define;
end
if (nargin<3)
    opts=[];
end
aux=[];
nbtc=length(dict.codel);
opts=filldefault(opts,'ifshow',0);
%
opts_used=opts;
split_strat=0;
cvecs=zeros(2,nbtc);
vec_made=zeros(1,nbtc);
let_targ=btc_vec2letcode(vec_targ,dict);
let_symok=let_targ;
let_symok.g=0;
targ.de=(let_targ.d+let_targ.e)/2;
targ.tu=(let_targ.t+let_targ.u)/2;
targ.vw=(let_targ.v+let_targ.w)/2;
%
let_symok.d=targ.de;
let_symok.e=targ.de;
let_symok.t=targ.tu;
let_symok.u=targ.tu;
let_symok.v=targ.vw;
let_symok.w=targ.vw;
vec_symok=btc_letcode2vec(let_symok,dict);
%
%try splits with gamma=0, non-Pickard theta=0
%
let=cell(0);
let{1}=let_symok;
let{2}=let_symok;
bxc=let_targ.b*let_targ.c;
let{1}.d=2*targ.de-bxc;
let{1}.e=bxc;
let{1}.t=0;
let{1}.u=2*targ.tu;
let{1}.v=0;
let{1}.w=2*targ.vw;
%
let{2}.d=bxc;
let{2}.e=2*targ.de-bxc;
let{2}.t=2*targ.tu;
let{2}.u=0;
let{2}.v=2*targ.vw;
let{2}.w=0;
%
wts=[0.5 0.5];
aux.wts=wts;
cvecs=[btc_letcode2vec(let{1},dict);btc_letcode2vec(let{2},dict)];
vec_made=[0.5 0.5]*cvecs;
corrs=cell(0);
p2x2=cell(0);
pa=zeros(2,2);
%
if_even=double(getp2x2_corrs(btc_vec2corrs([1:nbtc]==nbtc))>0);
parity=2*if_even(:)-1;
%
p2x2_az=cell(0); %probabilities with alpha=0
alpha_range_comp=zeros(2,2);
for icomp=1:2
    corrs{icomp}=btc_vec2corrs(cvecs(icomp,:),dict);
    p2x2{icomp}=getp2x2_corrs(corrs{icomp});
    corrs{icomp}=getcorrs_p2x2(p2x2{icomp},1); %this adds on the Pickard values
    p2x2_vec=p2x2{icomp}(:);
    pa(icomp,1)=min(p2x2_vec(parity==1));
    pa(icomp,2)=min(p2x2_vec(parity==-1));
    %
    %find the range of alpha that keeps all probabilities >=0
    p2x2_az{icomp}=getp2x2_corrs(setfield(corrs{icomp},'alpha',0));
    p2x2_az_vec=p2x2_az{icomp}(:);
    pa_az(icomp,1)=min(p2x2_az_vec(parity==1));
    pa_az(icomp,2)=min(p2x2_az_vec(parity==-1));
    %
    alpha_range_comp(icomp,1)=-8*pa_az(icomp,1);     %minimum value of alpha to raise probs with positive parities to zero
    alpha_range_comp(icomp,2)= 8*pa_az(icomp,2);     %maximum value of alpha that does not depress negative-parities below zero
    %
    aux.corrs{icomp}=corrs{icomp};
    aux.p2x2{icomp}=p2x2{icomp};
end
aux.wts=wts;
aux.alpha_range_comp=alpha_range_comp;
aux.alpha_range=[max(alpha_range_comp(:,1)),min(alpha_range_comp(:,2))];
if (opts.ifshow)
    disp('cvecs');
    disp(cvecs)
    %
    for icomp=1:2
        disp(sprintf(' component %1.0f',icomp));
        disp(sprintf(' minimum probability              of all %7.5f, of even %7.5f, of odd %7.5f',min(p2x2{icomp}(:)),pa(icomp,:)));
        disp(sprintf(' minimum probability with alpha=0 of all %7.5f, of even %7.5f, of odd %7.5f',min(p2x2_az{icomp}(:)),pa_az(icomp,:)));
        disp(sprintf(' alpha range: %7.5f %7.5f',alpha_range_comp(icomp,:)));
        disp(corrs{icomp}.cig_conds);
    end
    %
    disp('vec_made')
    disp(vec_made)
    disp('vec_made-vec_symok')
    disp(vec_made-vec_symok)  
end
return

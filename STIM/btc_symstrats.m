function [cvecs,split_strat,opts_used]=btc_symstrats(val_targ,dict,opts) 
% [cvecs,split_strat,opts_used]=btc_symstrats(val_targ,dict) determines btc coordinates of component textures
% that can be mixed to crate a target symmetric texture via "best" strategy
%
% val_targ: val_targ(1)=gamma, (2)=beta_hv, (3)=beta_diag, (4)=theta, (5)=alpha
% dict:  output of btc_define (created if not supplied)
% opts: opts.nsteps: number of steps for a gamma search (defaults to 9, relevant to strat 2 and 4)
%       opts.ifshow: show intermediate calcs (defaults to 0)
%       opts.split_alpha: enable asymmetric alpha values
%
% cvecs:  a set of row vectors (typically size=[2 10] of the coordinates of the component textures
%         the average of cvecs should correspond to val_targ
% split_strat: determined from calls to btc_symstrat
%  If not feasible (i.e., if cvecs requires negative probabilities), then split_strat=0.
%   1: 'gamma_beta_card': component textures have equal values of gamma and beta_card
%   2: 'gamma_beta_card_zero': component textures have equal values of gamma, each has one beta_hv=0
%   to add:
%   3: 'gamma_beta_card': component textures have equal values of beta_card,
%       but gamma adjusted to find largest min(prob)
%   4: 'gamma_beta_card_zero': component textures one beta_hv=0
%       but gamma adjusted to find largest min(prob)
%   possibly for future: allow for adjusting alpha
%   possibly for future: allow for mixtures that are not 1:1
%
% opts_used: options used
%
%   See also:  BTC_MAKESPEC_DEMO, BTC_SYMSTRAT, BTC_SEMSTRATS.
%
table=btc_symstrat();
nstrats=length(table);
if (nargin<2)
    dict=[];
end
if isempty(dict)
    dict=btc_define;
end
if (nargin<3)
    opts=[];
end
opts=filldefault(opts,'nsteps',9);
opts=filldefault(opts,'ifshow',0);
opts=filldefault(opts,'split_alpha',1);
%
opts_used=opts;
split_strat=0;
%
pmc_best=-Inf;
for istrat=1:nstrats
    [cvecs_try{istrat},ou{istrat}]=btc_symstrat(val_targ,istrat,dict,opts);
    for ispec=1:2
        corrs{ispec}=btc_vec2corrs(cvecs_try{istrat}(ispec,:),dict);
        p2x2{ispec}=getp2x2_corrs(corrs{ispec});
        minp2x2(ispec)=min(p2x2{ispec}(:));
    end
    if min(minp2x2)>pmc_best
        cvecs=cvecs_try{istrat};
        pmc_best=min(minp2x2);
        opts_used=ou{istrat};
        if pmc_best>=0
            split_strat=istrat;
        end
    end
end
return

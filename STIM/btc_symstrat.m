function [cvecs,opts_used]=btc_symstrat(val_targ,split_strat,dict,opts) 
% [cvecs,opts_used]=btc_symstrat(val_targ,split_strat,dict) determines btc coordinates of component textures
% that can be mixed to crate a target symmetric texture
%
% if called with no arguments, returns a description of the available strategies
%
% val_targ: val_targ(1)=gamma, (2)=beta_hv, (3)=beta_diag, (4)=theta, (5)=alpha
% split_strat
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
% dict:  output of btc_define (created if not supplied)
% opts: opts.nsteps: number of steps for a gamma search (defaults to 9, relevant to strat 2 and 4)
%       opts.ifshow: show intermediate calcs (defaults to 0)
%       opts.split_alpha: enable asymmetric alpha values
%
% cvecs:  a set of row vectors (typically size=[2 10] of the coordinates of the component textures
%         the average of cvecs should correspond to val_targ
% opts_used: options used
%
%  Does not check whether cvecs require negative probabilities
%
%   See also:  BTC_SYMSTRAT_DEMO, BTC_MAKESYM_DEMO, BTC_MAKEPICKARD, BTC_SPLIT_ALPHA, BTC_SYMSTRATS,  BTC_SEMSTRATS.
%
table=[];
table{1}.label='each texture has beta_cards equal';
table{1}.makepickard='gamma_beta_card';
table{2}.label='each texture has one beta_card=0, other nonzero';
table{2}.makepickard='gamma_beta_card_zero';
table{3}.label='each texture has beta_cards equal, gammas can be unequal';
table{3}.makepickard='gamma_beta_card';
table{4}.label='each texture has one beta_card=0, other nonzero, gammas can be unequal';
table{4}.makepickard='gamma_beta_card_zero';
if (nargin==0)
    cvecs=table;
    return
end
if (nargin<3)
    dict=[];
end
if isempty(dict)
    dict=btc_define;
end
if (nargin<4)
    opts=[];
end
opts=filldefault(opts,'nsteps',9);
opts=filldefault(opts,'ifshow',0);
opts=filldefault(opts,'split_alpha',1);
%
tol=10^-6;
%
opts_used=opts;
strat_name=table{split_strat}.makepickard;
if (split_strat==1) | (split_strat==2)
    nsteps=1;
else
    nsteps=opts.nsteps;
end
cvecs_try=cell(1,nsteps);
gamma_range=min(1-tol-val_targ(1),1-tol+val_targ(1));
bhv=strmatch('beta_hv',dict.name_order_aug);
for istep=1:nsteps
    % stepping through gamma, calling btc_makepickard for each half
    if (nsteps>1)
        if (strcmp(strat_name,'gamma_beta_card'))
            gamma_dev=gamma_range*(istep-1)/(nsteps-1);
        end
        if (strcmp(strat_name,'gamma_beta_card_zero'))
            gamma_dev=gamma_range*(2*(istep-1)/(nsteps-1)-1); %have to look at positive and negative values
            % since the beta_hv's are not symmetric
        end
    else
        gamma_dev=0;
    end
    gamma_val(1)=val_targ(1)-gamma_dev;
    gamma_val(2)=val_targ(1)+gamma_dev;
    cvecs_try{istep}=zeros(2,10);
    for ipick=1:2
        %ipick=1 is the texture that satisfies Pickard condition 1 and therefore is made with augcoords variant 2
        %ipick=2 is the texture that satisfies Pickard condition 2 and therefore is made with augcoords variant 1
        cvecs_try{istep}(ipick,strmatch('gamma',dict.name_order_aug,'exact'))=gamma_val(ipick);% gamma
        if (strcmp(strat_name,'gamma_beta_card'))
            cvecs_try{istep}(ipick,bhv)=val_targ(2); %beta_card
        end
        if (strcmp(strat_name,'gamma_beta_card_zero'))
            cvecs_try{istep}(ipick,bhv(ipick))=0;
            cvecs_try{istep}(ipick,bhv(3-ipick))=2*val_targ(2);
        end
        vin(1)=gamma_val(ipick);
        vin(2)=val_targ(2);
        vout=btc_makepickard(vin,strat_name,dict);
        inpickard=find(dict.inpickard(ipick,:));
        notpickard=find(~dict.inpickard(ipick,:));
        cvecs_try{istep}(ipick,intersect(inpickard,strmatch('beta_diag',dict.name_order_aug,'exact')))=...
            vout(1); %computed beta-diag
        cvecs_try{istep}(ipick,intersect(inpickard,strmatch('theta',dict.name_order_aug,'exact')))=...
            vout(2); %computed theta
        cvecs_try{istep}(3-ipick,intersect(inpickard,strmatch('beta_diag',dict.name_order_aug,'exact')))=...
            2*val_targ(3)-vout(1); %complementary beta_diag
        cvecs_try{istep}(3-ipick,intersect(inpickard,strmatch('theta',dict.name_order_aug,'exact')))=...
            2*val_targ(4)-vout(2); %computed theta
        cvecs_try{istep}(ipick,strmatch('alpha',dict.name_order_aug,'exact'))=val_targ(5);
        %
    end
    if (opts.ifshow==1)
        disp(sprintf(' before split_alpha: step %3.0f: gamma_dev=%7.4f, cvecs',istep,gamma_dev));
        disp(cvecs_try{istep});
    end
    if (opts.split_alpha)
        %
        [cvecs_try{istep},alpha_range]=btc_split_alpha(cvecs_try{istep},dict);
        %
        if (opts.ifshow==1)
            disp(sprintf(' after split_alpha: step %3.0f: gamma_dev=%7.4f, cvecs',istep,gamma_dev));
            disp(cvecs_try{istep});
            disp('alpha_range')
            disp(alpha_range)
        end
    end
    for ispec=1:2
        corrs{ispec}=btc_vec2corrs(cvecs_try{istep}(ispec,:),dict);
        p2x2{ispec}=getp2x2_corrs(corrs{ispec});
        minp2x2(istep,ispec)=min(p2x2{ispec}(:));
    end
end %istep
minp_both=min(minp2x2,[],2);
mp_ptr=min(find(minp_both==max(minp_both)));
if (opts.ifshow==1)
    disp('minp2x2')
    disp(minp2x2);
    disp(sprintf(' pointer to maximum of minp2x2: %4.0f',mp_ptr));
end
cvecs=cvecs_try{mp_ptr};
return


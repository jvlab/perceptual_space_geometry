% psg_umi_triplike_demo.m: analyze rank-choice probability data to
% determine max-likelihood Dirichlet prior and compatibility with symmetry and ultrametric inequality
% 
% For background, see .../jv/ey07977/psg_umi_notes.doc.
% 
% results are in the structure r:
%  r.dirichlet: fitting of choice probability distribution
%  r.su: log likelihood ratios for symmetry and ultrametric-ness
%
% private calculation (ipg=1): a and h calculated only from the triplets that meet threshold
% global calculation (ipg=2): a and h calculated from all triplets
%
% nfit_min is number of triplets required for the threshold to be used; set to 3 if not provided
%
% if_fixa=0: (default): a is fitted allowing h to vary, and also fitted with h drawn from h_fixlist
% if_fixa=1:            a is fitted allowing h to vary, but assigned a fixed value for all of h in h_fixlist
%
% 19Feb23: if_fast option to avoid recomputing triplet likelihoods for global case 
% 12Mar23: adapt to function version of psg_umi_triplike_plot[a]
% 20Mar23: add auto mode (if_auto=1; must set data_fullname and optionally set the structure auto)
% 20Mar23: add if_fixa
% 21Mar23: move ipg loop to outside and speed up by checking whether new threshold requires recalculation
%
% See also:  PSG_UMI_TRIPLIKE, PSG_TRIAD_STATS, PSG_UMI_STATS, PSG_TRIPLET_CHOICES, 
% PBETABAYES_COMPARE, LOGLIK_BETA, LOGLIK_BETA_DEMO2, PSG_READ_CHOICEDATA, PSG_UMI_TRIPLIKE_PLOT, PSG_UMI_TRIPLIKE_PLOTA.
%
if ~exist('if_auto') if_auto=0; end
auto=struct;
auto=filldefault(auto,'if_fast',1);
auto=filldefault(auto,'if_del',1);
auto=filldefault(auto,'if_plot',1);
auto=filldefault(auto,'if_plota',1);
auto=filldefault(auto,'if_fixa',0);
auto=filldefault(auto,'a_fixval',[]);
%
rng('default');
%for fitting dirichlet params
if ~exist('a_limits') a_limits=[2^-10 2^3]; end
if ~exist('h_limits') h_limits=[0 1]; end
if ~exist('h_init') h_init=0; end %initial value for h
if ~exist('h_fixlist') h_fixlist=[0 .001 .01 .1 .2 .4]; end %fixed nonzero values of h for fitting
h_fixlist=unique([0 h_fixlist(:)']);
nhfix=length(h_fixlist);
%
if if_auto
    if_fixa=auto.if_fixa;
    a_fixval=auto.a_fixval;
else
    a_fixval=[];
    if_fixa=getinp('1 to use fixed value for a','d',[0 1],0);
    if if_fixa
        a_fixval=getinp('value','f',[0 Inf]);
    end
end
%
opts_loglik=struct;
opts_loglik.qvec=0.5; %discrete part always is at 0.5
%
if ~exist('opts_triplike')
    opts_triplike=struct;
    opts_triplike.if_fast=[];
end
%
if ~exist('nfit_min') nfit_min=3;end %minimum number of data points (triads for dirichlet, or triplets for sym, umi) required to attempt a fit
%
if ~exist('sim_nstims') sim_nstims=25; end
if ~exist('sim_density') sim_density=0.5; end
if ~exist('sim_meantrials') sim_meantrials=5; end
%
if ~exist('data_fullname_def') data_fullname_def='./psg_data/bc6pt_choices_MC_sess01_10.mat'; end
%
thr_types={'min','max','avg'};
nthr_types=length(thr_types);
%
if if_auto
    if_fast=auto.if_fast;
else
    if_fast=getinp('1 for accelerated global computation and standard private, -1 to omit private','d',[-1 1]);
end
ipg_min=1;
if (if_fast==-1)
    ipg_min=2;
end
%
%read or simulate rank choice probabilities
%
if if_auto
    if_del=auto.if_del;
    if_plot=auto.if_plot;
    if_plota=auto.if_plota;
    underscore_sep=min(strfind(data_fullname,'_choices'));
    setup_fullname=cat(2,data_fullname(1:underscore_sep-1),'9.mat');
    [data,sa,opts_read_used]=psg_read_choicedata(data_fullname,setup_fullname);
    nstims=length(unique(data(:,[1:3])));
else
    if_del=getinp('1 to delete large variables','d',[0 1],0);
    if_plot=getinp('1 for detailed plots','d',[0 1]);
    if_plota=getinp('1 for summary (asymptotic) plots','d',[0 1]);
    switch getinp('0 to for random unstructured rank choice probabilities, 1 to read (-1 to skip reordering of stimuli)','d',[-1 1],1)
        case 1
            [data,sa,opts_read_used]=psg_read_choicedata([],[],setfields(struct(),{'data_fullname_def','if_log'},{data_fullname_def,1}));
            nstims=length(unique(data(:,[1:3])));
            data_fullname=opts_read_used.data_fullname;
        case -1
            data_fullname=getinp('full path and file name of data file','s',[],data_fullname_def);
            data=getfield(load(data_fullname),'responses');
            disp('data loaded but stimulus types not aligned with setup file.');
            sa=struct();
            nstims=length(unique(data(:,[1:3])));
        case 0
            data=nchoosek([1:sim_nstims],3);
            data=[data;data(:,[2 1 3]);data(:,[3 1 2])]; %permute but keep second column less than third
            data=data(rand(3*nchoosek(sim_nstims,3),1)<=sim_density,:);
            data=data(randperm(size(data,1)),:);
            for k=1:size(data,1)
                data(k,[2:3])=data(k,1+randperm(2));
            end
            sim_ntrials=poissrnd(sim_meantrials,size(data,1),1);
            nz=find(sim_ntrials>0);
            sim_ntrials=sim_ntrials(nz);
            data=data(nz,:);
            sim_probs=rand(size(data,1),1); %flat distribution of probabilities
            sim_ncloser=binornd(sim_ntrials,sim_probs);
            data=[data,sim_ncloser,sim_ntrials];
            sa=struct;
            nstims=sim_nstims;
            data_fullname='synthetic data';
    end
end
%
%report number of stimulus types
%
col_closer=4;
col_trials=5;
ntrials_found=sum(data(:,col_trials));
ntriads_found=size(data,1);
ustims_found=unique(reshape(data(:,[1:3]),[3*ntriads_found,1]));
nstims_found=length(ustims_found);
disp(sprintf('number of unique stimuli found: %3.0f; range from %3.0f to %3.0f',nstims_found,min(ustims_found),max(ustims_found)));
disp(sprintf('number of trials found: %6.0f',ntrials_found));
disp(sprintf('number of triads found: %6.0f',ntriads_found));
%
triplets=nchoosek([1:nstims],3); %triplets: unordered subsets of 3
ntriplets=nchoosek(nstims,3); %ntriplets: number of unordered subsets of 3
%
% ncloser: [ntriplets,3]: N(d(a,b)<d(a,c)), N(d(b,c)<d(b,a)), N(d(c,a)<d(c,b))
% ntrials: [ntriplets,3]: total trials in above
[ncloser,ntrials]=psg_triplet_choices(nstims,data); %extract triplets and sort
%
disp(sprintf('number of trials   after sorting: %6.0f',sum(ntrials(:)))) 
disp(sprintf('number of triads   after sorting: %6.0f',sum(ntrials(:)>0)));
disp(sprintf('number of triplets after sorting: %6.0f',sum(sum(ntrials,2)>0)));
if (sum(ntrials(:))~=ntrials_found)
    warning('mismatch with number of trials in data')
end
disp(sprintf('number of trials per triad range from %6.0f to %6.0f',min(ntrials(:)),max(ntrials(:))));
%
r=struct;
%find Dirichlet and Dirichlet/point mass parameters thresholding on a range of values for ntrials
%report min and max of ntrials, number available with each thresholding,
%best-fitting Dirichlet params, and log likelihood per triplet
%
r.h_fixlist=h_fixlist;
%
r.dirichlet=struct();
r.dirichlet.columns_tallies={'min_trials_per_triad','ntriads','ntrials'};
r.dirichlet.columns_a={'a','loglike_per_trial'};
r.dirichlet.columns_ah={'a','h','loglike_per_trial'};
r.dirichlet.a_fixval=a_fixval;
%
ithr=0;
%
for thr=min(ntrials(:)):max(ntrials(:))
    triads_use=find((ntrials(:)>=thr));
    ntriads_use=length(triads_use);
    ntrials_use=sum(ntrials(triads_use));
    if (ntriads_use>=nfit_min)
        ithr=ithr+1;
        r.dirichlet.tallies(ithr,:)=[thr,ntriads_use,ntrials_use];
        data_use=[ncloser(triads_use) ntrials(triads_use)];
        %fixed  values of h
        for ihfix=1:nhfix
            if if_fixa==0
                [fit_a,nll_a,exitflag_a]=fminbnd(@(x) -loglik_beta(x,data_use,setfield(opts_loglik,'hvec',h_fixlist(ihfix))),...
                    a_limits(1),a_limits(2)); %optimize assuming discrete part
            else
                fit_a=a_fixval;
                nll_a=-loglik_beta(fit_a,data_use,setfield(opts_loglik,'hvec',h_fixlist(ihfix)));
            end
            r.dirichlet.a(ithr,:,ihfix)=[fit_a,-nll_a/ntrials_use];
        end
        ah_init=[r.dirichlet.a(ithr,1,1);h_init]; %optimize with discrete part, using a_only fit as starting point
        [fit_ah,nll_ah,exitflag_ah,output_ah]=fminsearch(@(x) -loglik_beta(x(1),data_use,setfield(opts_loglik,'hvec',x(2))),ah_init);
        r.dirichlet.ah(ithr,:)=[fit_ah(:)',-nll_ah/ntrials_use];
    end
end
%now analyze for umi and symmetry
r.su.thr_types=thr_types;
r.su.columns_tallies={'thr','ntriplets','ntrials'};
r.su.columns_a={'a','loglike_per_trial'}; %values of a and h determined from the selected trials
r.su.columns_ah={'a','h','loglike_per_trial'}; %values of a and h determined from the selected trials
r.su.columns_sym={'llr_sym_vs_sym+notsym'}; %from likrat.sym of psg_umi_triplike
r.su.columns_umi={'llr_umi_trans_vs_trans'}; %from likrat.umi_trans of psg_umi_triplike
r.su.columns_sym_hfixed=r.su.columns_sym; %from likrat.sym of psg_umi_triplike
r.su.columns_umi_hfixed=r.su.columns_umi; %from likrat.umi_trans of psg_umi_triplike
r.su.global.a=r.dirichlet.a(1,1,:); % values with h fixed
%compute using global a and h from unthresholded Dirichlet and save in r.su.global.ah
if r.dirichlet.ah(1,2)>=0 %use full fit if h>=0
    r.su.global.ah=r.dirichlet.ah(1,1:2);
else %otherwise use best fit with h=0
    r.su.global.ah=[r.dirichlet.a(1,1,1),0];
end
%compute later using private a and h, to go in r.su.private.[a|ah]{ithr_type}
r.su.private.a=cell(1,nthr_types); 
r.su.private.ah=cell(1,nthr_types);
%
r.su.global_private_d1={'mean of sum','variance of sum'};
r.su.global_private_d2={'threshold type'};
r.su.global.sym=cell(2,nthr_types); %d1: mean of sum, d2: var of sum
r.su.global.umi=cell(2,nthr_types);
r.su.global.sym_hfixed=cell(2,nthr_types);
r.su.global.umi_hfixed=cell(2,nthr_types);
%
r.su.private.sym=cell(2,nthr_types);
r.su.private.umi=cell(2,nthr_types);
r.su.private.sym_hfixed=cell(2,nthr_types);
r.su.private.umi_hfixed=cell(2,nthr_types);
%
ncomps=3;
flipconfigs=int2nary([0:2^ncomps-1]',2);  %rows are [0 0 0;1 0 0;0 1 0;1 1 0; 0 0 1;1 0 1;0 1 1;1 1 1];
nflips=size(flipconfigs,1); %2^8
%
r.su.llr_d1={'threshold value'};
r.su.llr_d2={'orig data','flip all','flip any'};
r.su.llr_d3={'hfixed'};
nsurr=length(r.su.llr_d2); %three kinds of surrogates: native, flip all, flip any
llr_sym=cell(nsurr,2); %summed log likelihood ratio across trials, and summed variance of total 
llr_umi=cell(nsurr,2);
llr_sym_hfixed=cell(nsurr,2);
llr_umi_hfixed=cell(nsurr,2);
surr_list={1,[1 nflips],[1:nflips]};
if (if_fast~=0)
    %if_fast=1: calculate probabilities for all triplets
    loglik_rat_sym_all=zeros(ntriplets,nflips);
    loglik_rat_umi_all=zeros(ntriplets,nflips);
    loglik_rat_sym_hfixed_all=zeros(ntriplets,nflips,nhfix);
    loglik_rat_umi_hfixed_all=zeros(ntriplets,nflips,nhfix);
    ah=r.su.global.ah;
    ah_fixed=[squeeze(r.dirichlet.a(1,1,:)),h_fixlist(:)];
    for itriplet=1:ntriplets %accumulate likelihood ratios from each set of triplets
        obs_orig(:,1)=ncloser(itriplet,:)';
        obs_orig(:,2)=ntrials(itriplet,:)';
        obs_orig_flip=obs_orig(:,2)-obs_orig(:,1); 
        %
        for iflip=1:nflips %each surrogate
            obs=obs_orig;
            whichflip=find(flipconfigs(iflip,:)==1);
            obs(whichflip,1)=obs_orig_flip(whichflip);
            params.a=ah(1);
            params.h=ah(2);
            likrat=psg_umi_triplike(params,obs,opts_triplike);
            loglik_rat_sym_all(itriplet,iflip)=log(likrat.sym);
            loglik_rat_umi_all(itriplet,iflip)=log(likrat.umi_trans);
            for ihfix=1:nhfix
                params.a=ah_fixed(ihfix,1);
                params.h=ah_fixed(ihfix,2);
                likrat=psg_umi_triplike(params,obs,opts_triplike);
                loglik_rat_sym_hfixed_all(itriplet,iflip,ihfix)=log(likrat.sym);
                loglik_rat_umi_hfixed_all(itriplet,iflip,ihfix)=log(likrat.umi_trans);
            end
        end %iflip
    end %itriplet
    disp(sprintf(' global calculations done.'));
end
for ipg=ipg_min:2 %private and global
    for ithr_type=1:nthr_types %three kinds of thresholds: min, max, average
        if_ok=1;
        thr=0; %threshold
        ithr=1; %threshold pointer
        disp(sprintf('analyzing for symmetry and ultrametric likelihood ratio for threshold type %s',thr_types{ithr_type}));
        nuse_prev=-1; %will allow for reuse if increasing the threshold doesn't change the number of triplets/tents used
        while (if_ok)
            switch thr_types{ithr_type}
                case 'min'
                    triplets_use=find(min(ntrials,[],2)>=thr);
                    thr_val=thr;
                case 'max'
                    triplets_use=find(max(ntrials,[],2)>=thr);
                    thr_val=thr;
                case 'avg'
                    triplets_use=find(sum(ntrials,2)>=thr);
                    thr_val=thr/3; %average not total
            end
            if (length(triplets_use)>=nfit_min)
                ntriplets_use=length(triplets_use);
                if ntriplets_use~=nuse_prev
                    did_or_skipped='did'; %have to calculate
                    nuse_prev=ntriplets_use;
                    ntrials_use=sum(sum(ntrials(triplets_use,:)));
                    r.su.tallies{ithr_type}(ithr,:)=[thr_val ntriplets_use ntrials_use]; %threshold, number of triplets, number of trials
                    %compute private best-fitting a and h
                    data_use=[reshape(ncloser(triplets_use,:),3*ntriplets_use,1) reshape(ntrials(triplets_use,:),3*ntriplets_use,1)];
                    %fit with assuming fixed values of h
                    for ihfix=1:nhfix
                        if (if_fixa==0)
                            [fit_a,nll_a,exitflag_a]=fminbnd(@(x) -loglik_beta(x,data_use,setfield(opts_loglik,'hvec',h_fixlist(ihfix))),...
                                a_limits(1),a_limits(2)); %optimize assuming discrete part
                        else
                            fit_a=a_fixval;
                            nll_a=-loglik_beta(fit_a,data_use,setfield(opts_loglik,'hvec',h_fixlist(ihfix)));
                        end
                        r.su.private.a{ithr_type}(ithr,:,ihfix)=[fit_a,-nll_a/ntrials_use];
                    end
                    ah_init=[r.su.private.a{ithr_type}(ithr,1,1);h_init]; %optimize with discrete part, using a_only fit as starting point
                    [fit_ah,nll_ah,exitflag_ah,output_ah]=fminsearch(@(x) -loglik_beta(x(1),data_use,setfield(opts_loglik,'hvec',x(2))),ah_init);
                    if fit_ah(2)>=0
                        r.su.private.ah{ithr_type}(ithr,:)=[fit_ah(:)',-nll_ah/ntrials_use];
                    else
                        r.su.private.ah{ithr_type}(ithr,:)=[fit_a,0,-nll_a/ntrials_use];
                    end
                    %
                    %fast global option:calculate probabilities for all triplets and later select
                    %
                    if if_fast~=0 & ipg==2
                        loglik_rat_sym=loglik_rat_sym_all(triplets_use,:);
                        loglik_rat_umi=loglik_rat_umi_all(triplets_use,:);
                        loglik_rat_sym_hfixed=loglik_rat_sym_hfixed_all(triplets_use,:,:);
                        loglik_rat_umi_hfixed=loglik_rat_umi_hfixed_all(triplets_use,:,:);
                    else %if_fast==0
                        if (ipg==1) %private
                            ah=r.su.private.ah{ithr_type}(ithr,:); %a and h both fitted
                            ah_fixed=[squeeze(r.su.private.a{ithr_type}(ithr,1,:)),h_fixlist(:)]; %a fitted, h fixed
                        else %global
                            ah=r.su.global.ah;
                            ah_fixed=[squeeze(r.dirichlet.a(1,1,:)),h_fixlist(:)];
                        end
                        loglik_rat_sym=zeros(ntriplets_use,nflips);
                        loglik_rat_umi=zeros(ntriplets_use,nflips);
                        loglik_rat_sym_hfixed=zeros(ntriplets_use,nflips,nhfix);
                        loglik_rat_umi_hfixed=zeros(ntriplets_use,nflips,nhfix);
                        for itriplet=1:ntriplets_use %accumulate likelihood ratios from each set of triplets
                            obs_orig(:,1)=ncloser(triplets_use(itriplet),:)';
                            obs_orig(:,2)=ntrials(triplets_use(itriplet),:)';
                            obs_orig_flip=obs_orig(:,2)-obs_orig(:,1); 
                            %
                            for iflip=1:nflips %each surrogate
                                obs=obs_orig;
                                whichflip=find(flipconfigs(iflip,:)==1);
                                obs(whichflip,1)=obs_orig_flip(whichflip);
                                params.a=ah(1);
                                params.h=ah(2);
                                likrat=psg_umi_triplike(params,obs,opts_triplike);
                                loglik_rat_sym(itriplet,iflip)=log(likrat.sym);
                                loglik_rat_umi(itriplet,iflip)=log(likrat.umi_trans);
                                for ihfix=1:nhfix
                                    params.a=ah_fixed(ihfix,1);
                                    params.h=ah_fixed(ihfix,2);
                                    likrat=psg_umi_triplike(params,obs,opts_triplike);
                                    loglik_rat_sym_hfixed(itriplet,iflip,ihfix)=log(likrat.sym);
                                    loglik_rat_umi_hfixed(itriplet,iflip,ihfix)=log(likrat.umi_trans);
                                end
                            end %iflip
                        end
                    end %if_fast
                    %do statistics
                    for isurr=1:nsurr
                        surr_sel=surr_list{isurr}; %for isurr=1, this is just the original data (1)
                        llr_sym{isurr,1}=sum(mean(loglik_rat_sym(:,surr_sel),2),1);
                        llr_umi{isurr,1}=sum(mean(loglik_rat_umi(:,surr_sel),2),1);
                        llr_sym_hfixed{isurr,1}=reshape(sum(mean(loglik_rat_sym_hfixed(:,surr_sel,:),2),1),[1 1 nhfix]);
                        llr_umi_hfixed{isurr,1}=reshape(sum(mean(loglik_rat_umi_hfixed(:,surr_sel,:),2),1),[1 1 nhfix]);
                        %each triplet contributes independently to the variance
                        %variance for each triplet is normalized by N not N-1, since we have all the values
                        llr_sym{isurr,2}=sum(var(loglik_rat_sym(:,surr_sel),1,2),1);
                        llr_umi{isurr,2}=sum(var(loglik_rat_umi(:,surr_sel),1,2),1);
                        llr_sym_hfixed{isurr,2}=reshape(sum(var(loglik_rat_sym_hfixed(:,surr_sel,:),1,2),1),[1 1 nhfix]);
                        llr_umi_hfixed{isurr,2}=reshape(sum(var(loglik_rat_umi_hfixed(:,surr_sel,:),1,2),1),[1 1 nhfix]);
                        %
                        for imv=1:2% mean and variance
                            if (ipg==1)
                                r.su.private.sym{imv,ithr_type}(ithr,isurr)=llr_sym{isurr,imv};
                                r.su.private.umi{imv,ithr_type}(ithr,isurr)=llr_umi{isurr,imv};
                                r.su.private.sym_hfixed{imv,ithr_type}(ithr,isurr,:)=llr_sym_hfixed{isurr,imv};
                                r.su.private.umi_hfixed{imv,ithr_type}(ithr,isurr,:)=llr_umi_hfixed{isurr,imv};
                            else
                                r.su.global.sym{imv,ithr_type}(ithr,isurr)=llr_sym{isurr,imv};
                                r.su.global.umi{imv,ithr_type}(ithr,isurr)=llr_umi{isurr,imv};
                                r.su.global.sym_hfixed{imv,ithr_type}(ithr,isurr,:)=llr_sym_hfixed{isurr,imv};
                                r.su.global.umi_hfixed{imv,ithr_type}(ithr,isurr,:)=llr_umi_hfixed{isurr,imv};
                            end
                        end %imv
                    end %isurr
                else
                    did_or_skipped='skp';
                    r.su.tallies{ithr_type}(ithr,:)=r.su.tallies{ithr_type}(ithr-1,:);
                    r.su.tallies{ithr_type}(ithr,1)=thr_val; %threshold is new
                    if (ipg==1)
                        r.su.private.a{ithr_type}(ithr,:,:)=r.su.private.a{ithr_type}(ithr-1,:,:);
                        r.su.private.ah{ithr_type}(ithr,:)=r.su.private.ah{ithr_type}(ithr-1,:);
                    end
                    for isurr=1:nsurr
                        for imv=1:2% mean and variance
                            if (ipg==1)
                                r.su.private.sym{imv,ithr_type}(ithr,isurr)=llr_sym{isurr,imv};
                                r.su.private.umi{imv,ithr_type}(ithr,isurr)=llr_umi{isurr,imv};
                                r.su.private.sym_hfixed{imv,ithr_type}(ithr,isurr,:)=llr_sym_hfixed{isurr,imv};
                                r.su.private.umi_hfixed{imv,ithr_type}(ithr,isurr,:)=llr_umi_hfixed{isurr,imv};
                            else
                                r.su.global.sym{imv,ithr_type}(ithr,isurr)=llr_sym{isurr,imv};
                                r.su.global.umi{imv,ithr_type}(ithr,isurr)=llr_umi{isurr,imv};
                                r.su.global.sym_hfixed{imv,ithr_type}(ithr,isurr,:)=llr_sym_hfixed{isurr,imv};
                                r.su.global.umi_hfixed{imv,ithr_type}(ithr,isurr,:)=llr_umi_hfixed{isurr,imv};
                            end
                        end %imv
                    end %isurr
                end %nuse_prev
                disp(sprintf('%s ipg %3.0f ithr_type %3.0f ithr %3.0f thr %3.0f ntriplets_use %6.0f size(loglik_rat_sym) %6.0f %4.0f size(loglik_rat_sym_hfixed) %6.0f %4.0f %4.0f',...
                    did_or_skipped,ipg,ithr_type,ithr,thr,ntriplets_use,size(loglik_rat_sym),size(loglik_rat_sym_hfixed)));
                thr=thr+1; %threshold
                ithr=ithr+1; %pointer
            else
                if_ok=0;
            end
        end %if_ok
    end %thr_type
end %ipg
%
%finish and do plots
%
if (if_del)
    clear *all
end
plot_opts=[];
plot_opts.ipg_min=ipg_min;
plot_opts.data_fullname=data_fullname;
plot_opts.llr_field='su';
if (if_plot)
    psg_umi_triplike_plot(r,plot_opts);
end
if (if_plota)
    psg_umi_triplike_plota(r,plot_opts);
end
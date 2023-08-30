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
% if_fixa=2:            a is fixed for global fits when h is allowed to vary or when h is drawn from h_fixlist
%                       but allowed to vary for private when h can vary and best fit for h>=0
%
% 19Feb23: if_fast option to avoid recomputing triplet likelihoods for global case 
% 12Mar23: adapt to function version of psg_umi_triplike_plot[a]
% 20Mar23: add auto mode (if_auto=1; must set data_fullname and optionally set the structure auto)
% 20Mar23: add if_fixa
% 21Mar23: move ipg loop to outside and speed up by checking whether new threshold requires recalculation
% 04Apr23: automate saving results in a database
% 05Apr23: allow for external control of plot_opts.frac_keep_list
% 02May23: adding 'flip_one' conform surrogate, change defaults
% 04May23: code cleanup (ipg_string)
% 11May23: add psg_choicedata_makeeven
% 15Jun23: add standard error of the mean for orig data and conform surrogate
% 23Jun23: add psg_select_choicedata.  Note that the field name in db includes the selection descriptor, if present
% 27Jul23: add to prompt that when no metadata are present (no reordering of stimuli), then ncloser is assumed for columm 4 of data
% 31Jul23: modify so that when no metadata are present, ncloser is properly determined based on responses_colnames
% 01Aug23: all reads go through psg_read_choicedata
% 11Aug23: modify with if_fixa=2
% 30Aug23: change defaults
%
% See also:  PSG_UMI_TRIPLIKE, PSG_TRIAD_STATS, PSG_UMI_STATS, PSG_TRIPLET_CHOICES, 
% PBETABAYES_COMPARE, LOGLIK_BETA, LOGLIK_BETA_DEMO2, PSG_READ_CHOICEDATA, PSG_CHOICEDATA_MAKEEVEN,
% PSG_UMI_TRIPLIKE_PLOT, PSG_UMI_TRIPLIKE_PLOTA, PSG_UMI_TRIP_TENT_RUN,PSG_INEQ_LOGIC, PSG_CONFORM, PSG_SELECT_CHOICEDATA.
%
if ~exist('if_auto') if_auto=0; end
if ~exist('auto')
    auto=struct;
end
auto=filldefault(auto,'if_fast',1);
auto=filldefault(auto,'if_del',1);
auto=filldefault(auto,'if_plot',1);
auto=filldefault(auto,'if_plota',1);
auto=filldefault(auto,'if_fixa',0);
auto=filldefault(auto,'a_fixval',[]);
auto=filldefault(auto,'if_reorder',1);
auto=filldefault(auto,'if_conform',1);
auto=filldefault(auto,'db_file','.\psg_data\psg_tentlike_db.mat');
auto=filldefault(auto,'opts_conform',struct());
auto=filldefault(auto,'if_makeeven',0);
auto=filldefault(auto,'opts_makeeven',struct());
auto=filldefault(auto,'sel_apply',0);
auto=filldefault(auto,'sel_string',[]);
auto=filldefault(auto,'sel_desc',[]);
auto=filldefault(auto,'sel_opts',struct());
if ~exist('sel_opts') sel_opts=struct();end
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
    if_makeeven=auto.if_makeeven;
    disp('running in auto mode, with settings')
    auto
else
    a_fixval=[];
    if_fixa=getinp('1 to use fixed value for a when h varies, 2 for fixed a when h varies or fixed but global)','d',[0 2],0);
    if if_fixa
        a_fixval=getinp('value','f',[0 Inf]);
    end
    if_makeeven=getinp('1 to reduce all triads to an even number of trials','d',[0 1],0);
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
    if_fast=getinp('1 for accelerated global computation and standard private, -1 to omit private','d',[-1 1],-1);
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
    if (auto.if_reorder>1)
        [data,sa,opts_read_used]=psg_read_choicedata(data_fullname,setup_fullname);
    else
        [data,sa,opts_read_used]=psg_read_choicedata(data_fullname,setup_fullname,setfield([],'nometa',1));
    end
    nstims=length(unique(data(:,[1:3])));
    if_conform=auto.if_conform;
else
    if_del=getinp('1 to delete large variables','d',[0 1],1);
    if_plot=getinp('1 for detailed plots','d',[0 1],0);
    if_plota=getinp('1 for summary and summary (asymptotic) plots','d',[0 1],1);
    switch getinp('0 to generate random unstructured rank choice probabilities, 1 to read data (-1 to skip reordering of stimuli)','d',[-1 1],-1)
        case 1
            [data,sa,opts_read_used]=psg_read_choicedata([],[],setfields(struct(),{'data_fullname_def','if_log'},{data_fullname_def,1}));
            nstims=length(unique(data(:,[1:3])));
            data_fullname=opts_read_used.data_fullname;
        case -1
            [data,sa,opts_read_used]=psg_read_choicedata([],[],setfields(struct(),{'data_fullname_def','if_log','nometa'},{data_fullname_def,1,1}));
            nstims=length(unique(data(:,[1:3])));
            data_fullname=opts_read_used.data_fullname;
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
    if_conform=getinp('1 to analyze conforming surrogates','d',[0 1],1);
end
%apply a selection?
if (if_auto)
    sel_apply=auto.sel_apply;
    sel_string=auto.sel_string;
    sel_desc=auto.sel_desc;
    sel_opts=auto.sel_opts;
else
    sel_apply=getinp('1 to apply a selection criterion','d',[0 1],0);
    sel_string=[];
    sel_desc=[];
end
if (sel_apply)
    [data,sa,sel_string,sel_desc,sel_opts_used]=psg_select_choicedata(data,sa,sel_string,sel_desc,sel_opts);
    nstims=sel_opts_used.nstims;
end
%
if (if_makeeven)
    if if_auto
        opts_makeeven=auto.opts_makeeven;
    end
    if ~exist('opts_makeeven')
        opts_makeeven=struct;
    end
    [data,opts_makeeven_used]=psg_choicedata_makeeven(data,opts_makeeven);
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
r.nstims=nstims;
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
        if (if_fixa==2)
            [fit_h,nll_h,exitflag_h]=fminbnd(@(x) -loglik_beta(a_fixval,data_use,setfield(opts_loglik,'hvec',x)),...
                 0,1); %assume a, optimize discrete part
            r.dirichlet.ah(ithr,:)=[a_fixval,fit_h,-nll_h/ntrials_use];
        else
            ah_init=[r.dirichlet.a(ithr,1,1);h_init]; %optimize with discrete part, using a_only fit as starting point
            [fit_ah,nll_ah,exitflag_ah,output_ah]=fminsearch(@(x) -loglik_beta(x(1),data_use,setfield(opts_loglik,'hvec',x(2))),ah_init);
            r.dirichlet.ah(ithr,:)=[fit_ah(:)',-nll_ah/ntrials_use];
        end
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
%
if (if_conform)
    nconform=1; %only one kind of conforming 
    %set up for flip_one
    if if_auto
        opts_conform=auto.opts_conform;
    end
    if ~exist('opts_conform')
        opts_conform=struct;
    end
    opts_conform=filldefault(opts_conform,'method','flip_one');
    opts_conform=filldefault(opts_conform,'penalty','chi2');
    opts_conform=filldefault(opts_conform,'if_log',1);
    %
    logic_type=struct;
    logic_type.sym='exclude_sym';
    logic_type.umi='exclude_umi_trans';
    logic_types=fieldnames(logic_type);
    %
    nlogic_conform=length(logic_types); %sym and umi
    partitions_conform=struct();
    ncloser_conform=struct();
    if_flip_conform=struct();
    which_flip_conform=struct();
    opts_conform_used=struct();
    %
    for ilc=1:nlogic_conform
        lt=logic_types{ilc};
        disp(sprintf('flip-one analysis: %s, using exclusion logic of %s',lt,logic_type.(lt)));
        partitions_conform.(lt)=psg_ineq_logic(ncomps,logic_type.(lt));
        [ncloser_conform.(lt),if_flip_conform.(lt),opts_conform_used.(lt)]=psg_conform(ncloser,ntrials,partitions_conform.(lt),opts_conform);
        which_flip_conform.(lt)=1+if_flip_conform.(lt)*(2.^[0:ncomps-1]'); %points to list of all possible flips
    end
    r.conform_results=opts_conform_used;
    r.conform_results.which_flip_conform=which_flip_conform;
    r.conform_results.logic_type=logic_type;
    r.su.llr_d2{end+1}=strrep(opts_conform.method,'_',' ');
else
    nconform=0;
end
r.nsurr=nsurr;
r.nconform=nconform;
%
llr_sym=cell(nsurr+nconform,2); %summed log likelihood ratio across trials, and summed variance of total 
llr_umi=cell(nsurr+nconform,2);
llr_sym_hfixed=cell(nsurr+nconform,2);
llr_umi_hfixed=cell(nsurr+nconform,2);
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
    disp(sprintf('overall global calculations done.'));
end
ipg_strings={'private','global'};
for ipg=ipg_min:2 %private and global
    disp(sprintf('%10s calculations',ipg_strings{ipg}));
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
                    for isurr=1:nsurr+nconform
                        if (isurr<=nsurr)
                            surr_sel=surr_list{isurr}; %for isurr=1, this is just the original data (1)
                            llr_sym{isurr,1}=sum(mean(loglik_rat_sym(:,surr_sel),2),1);
                            llr_umi{isurr,1}=sum(mean(loglik_rat_umi(:,surr_sel),2),1);
                            llr_sym_hfixed{isurr,1}=reshape(sum(mean(loglik_rat_sym_hfixed(:,surr_sel,:),2),1),[1 1 nhfix]);
                            llr_umi_hfixed{isurr,1}=reshape(sum(mean(loglik_rat_umi_hfixed(:,surr_sel,:),2),1),[1 1 nhfix]);
                            if (isurr>1)
                                %each triplet contributes independently to the variance
                                %variance for each triplet is normalized by N not N-1, since we have all the values
                                llr_sym{isurr,2}=sum(var(loglik_rat_sym(:,surr_sel),1,2),1);
                                llr_umi{isurr,2}=sum(var(loglik_rat_umi(:,surr_sel),1,2),1);
                                llr_sym_hfixed{isurr,2}=reshape(sum(var(loglik_rat_sym_hfixed(:,surr_sel,:),1,2),1),[1 1 nhfix]);
                                llr_umi_hfixed{isurr,2}=reshape(sum(var(loglik_rat_umi_hfixed(:,surr_sel,:),1,2),1),[1 1 nhfix]);
                            else %isurr=1: original data. Here, goal is for psg_umi_triplike_plota to compute standard error of the mean
                                %which is sqrt(var)/ntriplets_use, but
                                %psg_umi_triplike_plota will find square root and then divide by ntriplets_use
                                %so here we just compute var, normalized by N-1 since it is a sample
                                %here, surr_sel=1
                                llr_sym{isurr,2}=var(loglik_rat_sym(:,surr_sel),0,1);
                                llr_umi{isurr,2}=var(loglik_rat_umi(:,surr_sel),0,1);
                                llr_sym_hfixed{isurr,2}=reshape(var(loglik_rat_sym_hfixed(:,surr_sel,:),0,1),[1 1 nhfix]);
                                llr_umi_hfixed{isurr,2}=reshape(var(loglik_rat_umi_hfixed(:,surr_sel,:),0,1),[1 1 nhfix]);
                            end
                        else %do conform
                            %select the appropriate flip for each triplet
                            loglik_rat_sym_conform=zeros(ntriplets_use,1);
                            loglik_rat_umi_conform=zeros(ntriplets_use,1);
                            loglik_rat_sym_hfixed_conform=zeros(ntriplets_use,nhfix);
                            loglik_rat_umi_hfixed_conform=zeros(ntriplets_use,nhfix);
                            for itptr=1:ntriplets_use
                                it=triplets_use(itptr);
                                loglik_rat_sym_conform(itptr)=loglik_rat_sym(itptr,which_flip_conform.sym(it));
                                loglik_rat_umi_conform(itptr)=loglik_rat_umi(itptr,which_flip_conform.umi(it));
                                loglik_rat_sym_hfixed_conform(itptr,:)=reshape(loglik_rat_sym_hfixed(itptr,which_flip_conform.sym(it),:),[1 nhfix]);
                                loglik_rat_umi_hfixed_conform(itptr,:)=reshape(loglik_rat_umi_hfixed(itptr,which_flip_conform.umi(it),:),[1 nhfix]);
                            end
                            llr_sym{isurr,1}=sum(loglik_rat_sym_conform);
                            llr_umi{isurr,1}=sum(loglik_rat_umi_conform);
                            llr_sym_hfixed{isurr,1}=reshape(sum(loglik_rat_sym_hfixed_conform),[1 1 nhfix]);
                            llr_umi_hfixed{isurr,1}=reshape(sum(loglik_rat_umi_hfixed_conform),[1 1 nhfix]);
                            %variances are calculated as for original data, but surr_sel must be set
                            surr_sel=1;
                            llr_sym{isurr,2}=var(loglik_rat_sym(:,surr_sel),0,1);
                            llr_umi{isurr,2}=var(loglik_rat_umi(:,surr_sel),0,1);
                            llr_sym_hfixed{isurr,2}=reshape(var(loglik_rat_sym_hfixed(:,surr_sel,:),0,1),[1 1 nhfix]);
                            llr_umi_hfixed{isurr,2}=reshape(var(loglik_rat_umi_hfixed(:,surr_sel,:),0,1),[1 1 nhfix]);
                        end
                        %
                        for imv=1:2% mean and variance
                            r.su.(ipg_strings{ipg}).sym{imv,ithr_type}(ithr,isurr)=llr_sym{isurr,imv};
                            r.su.(ipg_strings{ipg}).umi{imv,ithr_type}(ithr,isurr)=llr_umi{isurr,imv};
                            r.su.(ipg_strings{ipg}).sym_hfixed{imv,ithr_type}(ithr,isurr,:)=llr_sym_hfixed{isurr,imv};
                            r.su.(ipg_strings{ipg}).umi_hfixed{imv,ithr_type}(ithr,isurr,:)=llr_umi_hfixed{isurr,imv};
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
                    for isurr=1:nsurr+nconform
                        for imv=1:2% mean and variance
                            r.su.(ipg_strings{ipg}).sym{imv,ithr_type}(ithr,isurr)=llr_sym{isurr,imv};
                            r.su.(ipg_strings{ipg}).umi{imv,ithr_type}(ithr,isurr)=llr_umi{isurr,imv};
                            r.su.(ipg_strings{ipg}).sym_hfixed{imv,ithr_type}(ithr,isurr,:)=llr_sym_hfixed{isurr,imv};
                            r.su.(ipg_strings{ipg}).umi_hfixed{imv,ithr_type}(ithr,isurr,:)=llr_umi_hfixed{isurr,imv};
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
if ~exist('plot_opts') %allow for setting plot_opts.frac_keep_list
    plot_opts=struct;
end
plot_opts.ipg_min=ipg_min;
plot_opts.data_fullname=data_fullname;
plot_opts.llr_field='su';
plot_opts.nconform=nconform;
plot_opts.nsurr=nsurr;
if ~isempty(sel_desc)
    plot_opts.sel_desc=sel_desc;
end
if (if_plot)
    psg_umi_triplike_plot(r,plot_opts);
end
if (if_plota) | (if_auto)
    [plot_opts_used,figh,s]=psg_umi_triplike_plota(r,plot_opts);
    if (if_auto)
        if exist(auto.db_file,'file')
            db=getfield(load(auto.db_file),'db');
        else
            db=struct;
        end
        data_shortname=data_fullname;
        data_shortname=strrep(data_shortname,'.mat','');
        data_shortname=strrep(data_shortname,'/',filesep);
        data_shortname=strrep(data_shortname,'\',filesep);
        data_shortname=cat(2,filesep,data_shortname);
        data_shortname=data_shortname(1+max(find(data_shortname==filesep)):end);
        if isempty(sel_desc)
            data_fieldname=data_shortname;
        else
            data_fieldname=cat(2,data_shortname,'_',sel_desc);
        end
        db.(data_fieldname).r=r;
        db.(data_fieldname).s=s;
        db.(data_fieldname).select.sel_string=sel_string;
        db.(data_fieldname).select.sel_desc=sel_desc;
        save(auto.db_file,'db');
        disp(sprintf('saved results from %s in %s',data_fieldname,auto.db_file));
    end
end

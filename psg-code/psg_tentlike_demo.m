% psg_tentlike_demo.m: analyze rank-choice probability data of tent configurations 
% determine max-likelihood Dirichlet prior and compatibility necessary inequalities for addtree
% 
% For background, see .../jv/ey07977/psg_umi_notes.doc.
% 
% structured like psg_umi_triplike_demo
% in contrast to the counts in psg_tent_stats, here the tripods includes
% judgments in which not all members of the triangle are in the same trial.
%
% Tents:
% first 3 components a tripod (R(z;b,c), R(z;c,a), R(z;a,b))
% last 3 components are a triplet (R(a;b,c), R(b;c,a), R(c;a,b))
% tripods are determined from triplet choices (via psg_triplet_choices), so they will include
% data from all three judgments from same trial, as well as data from separate trials.
%
% results are in the structure r:
%  r.dirichlet: fitting of choice probability distribution
%  r.ad: log likelihood ratios for addtree necesary condition
%
% private calculation (ipg=1): a and h calculated only from the tents that meet threshold
% global calculation (ipg=2): a and h calculated from all tents
%
% nfit_min is number of tents required for the threshold to be used; set to 3 if not provided
%
% test files: psg_psg_tentlike_test_4stims_addtree_*.mat, four stimuli
%    *= bal and balweak are unambiguously addtree
%    *= full is close to being addtree, there are choice probs close to 1/2 that need to be flipped
%    *= yes has missing data and is inconsistent with transitivity, so cannot be made addtree
%
% if_fixa=0: (default): a is fitted allowing h to vary, and also fitted with h drawn from h_fixlist
% if_fixa=1:            a is fitted allowing h to vary, but assigned a fixed value for all of h in h_fixlist
%
% 20Mar23: add auto mode (if_auto=1; must set data_fullname and optionally set the structure auto)
% 20Mar23: add if_fixa
% 21Mar23: move ipg loop to outside and speed up by checking whether new threshold requires recalculation
% 04Apr23: automate saving results in a database
% 05Apr23: allow for external control of plot_opts.frac_keep_list
% 04May23: start adding 'flip_one' conform surrogate, change defaults, code cleanup (ipg_string)
% 11May23: add psg_choicedata_makeeven
% 15Jun23: add standard error of the mean for orig data and conform surrogate
% 23Jun23: add psg_select_choicedata.  Note that the field name in db includes the selection descriptor, if present
% 27Jul23: add to prompt that when no metadata are present (no reordering of stimuli), then ncloser is assumed for columm 4 of data
% 31Jul23: modify so that when no metadata are present, ncloser is properly determined based on responses_colnames
% 01Aug23: all reads go through psg_read_choicedata
%
% See also:  PSG_UMI_TRIPLIKE_DEMO, PSG_TENT_STATS, PSG_TRIPLET_CHOICES, 
% LOGLIK_BETA, LOGLIK_BETA_DEMO2, PSG_READ_CHOICEDATA, PSG_CHOICEDATA_MAKEEVEN, PSG_UMI_TRIPLIKE_PLOTA, NCHOOSEK2SEQ_3VR,
% PSG_INEQ_LOGIC, PSG_PERMUTES_LOGIC, PSG_INEQ_APPLY, PSG_UMI_TRIP_TENT_RUN, PSG_CONFORM, PSG_SELECT_CHOICEDATA.
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
%
ineq_logic_types={'exclude_addtree_trans','exclude_trans_tent'};
nineq=length(ineq_logic_types);
%
ncomps=6; %six rank choice probabilities to be compared
%
partitions=cell(0);
for ineq=1:nineq
    disp(sprintf('creating inequality logic for %s',ineq_logic_types{ineq}));
    partitions{ineq}=psg_ineq_logic(ncomps,ineq_logic_types{ineq},setfield([],'if_log',1));
end
disp('creating permutation logic for surrogates')
permutes=psg_permutes_logic(ncomps,'flip_each');
disp(sprintf(' size is %3.0f x %3.0f',size(permutes)));
nflips=size(permutes,2);
%
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
    if_fixa=getinp('1 to use fixed value for a','d',[0 1],0);
    if if_fixa
        a_fixval=getinp('value','f',[0 Inf]);
    end
    if_makeeven=getinp('1 to reduce all triads to an even number of trials','d',[0 1],0);
end
%
opts_loglik=struct;
opts_loglik.qvec=0.5; %discrete part always is at 0.5
%
if ~exist('opts_tent')
    opts_tent=struct;
end
%
if ~exist('nfit_min') nfit_min=3;end %minimum number of data points (triads for dirichlet, or six comparisons for tent) required to attempt a fit
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
    if (auto.if_reorder>1)
        [data,sa,opts_read_used]=psg_read_choicedata(data_fullname,setup_fullname);
    else
        [data,sa,opts_read_used]=psg_read_choicedata(data_fullname,setup_fullname,setfield([],'nometa',1));
    end
    nstims=length(unique(data(:,[1:3])));
    if_conform=auto.if_conform;
else
    %
    if_del=getinp('1 to delete large variables','d',[0 1],1);
    if_plot=getinp('1 for detailed plots','d',[0 1],0);
    if_plota=getinp('1 for summary (asymptotic) plots','d',[0 1]);
    switch getinp('0 to generate random unstructured rank choice probabilities, 1 to read data (-1 to skip reordering of stimuli)','d',[-1 1],1)
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
nt=3; %number of points in a triangle
triplets=nchoosek([1:nstims],nt); %triplets: unordered subsets of 3
ntriplets=nchoosek(nstims,nt); %ntriplets: number of unordered subsets of 3
%
% ncloser_triplets: [ntriplets,3]: N(d(a,b)<d(a,c)), N(d(b,c)<d(b,a)), N(d(c,a)<d(c,b))
% ntrials_triplets: [ntriplets,3]: total trials in above
[ncloser_triplets,ntrials_triplets,abc_list]=psg_triplet_choices(nstims,data); %extract triplets and sort
%
disp(sprintf('number of trials   after sorting: %6.0f',sum(ntrials_triplets(:)))) 
disp(sprintf('number of triads   after sorting: %6.0f',sum(ntrials_triplets(:)>0)));
disp(sprintf('number of triplets after sorting: %6.0f',sum(sum(ntrials_triplets,2)>0)));
if (sum(ntrials_triplets(:))~=ntrials_found)
    warning('mismatch with number of trials in data')
end
disp(sprintf('number of trials per triad range from %6.0f to %6.0f',min(ntrials_triplets(:)),max(ntrials_triplets(:))));
%
%now create tents from the triplets
%
ntriplets_exclude=nchoosek(nstims-1,nt); %number of triplets that exclude a given stimulus
ntents=nstims*ntriplets_exclude;
tents=zeros(ntents,nt+1);
%set up table of all possible tents: z,a,b,c, where a,b,c are distinct from z, and a<b<c
for istim=1:nstims
    tents([1:ntriplets_exclude]+(istim-1)*ntriplets_exclude,:)=[repmat(istim,ntriplets_exclude,1) nchoosek(setdiff([1:nstims],istim),nt)];
end
ncloser=zeros(ntents,nt*2); 
ntrials=zeros(ntents,nt*2);
disp('creating tents from triplets');
%abc_lists=tents(:,1+[1:nt]);
%
% for ncloser, ntrials, there is one row for every tent.  Rows correspond to rows of tents, which specify z,a,b,c
% ncloser: [ntriplets,6]: N(d(z,b)<d(z,c)), N(d(z,c)<d(z,a)), N(d(z,a)< d(z,b), N(d(a,b)<d(a,c)), N(d(b,c)<d(b,a)), N(d(c,a)<d(c,b))
% ntrials: [ntriplets,6]: total trials in above
%
% last 3 columns of close, ntrials are from triplets
triplet_rows=nchoosek2seq_3vr(tents(:,1+[1:nt]),nstims); %find rows in triplet table that match cols 2-4 
ncloser(:,nt+[1:nt])=ncloser_triplets(triplet_rows,:);
ntrials(:,nt+[1:nt])=ntrials_triplets(triplet_rows,:);
%
% first 3 columns of nclose, ntrials are the tripod
%note that col 1 and subsequent cols of tent are not necessarily in ascending order
%so we can't simply extract from the triplet table!!!!!
%
% a, b, and c are always i norder a<b<c in triplet table and tent table
%this loop takes care of odd permutation if z is in middle of the selected pair of a,b,c
%then we adjust second column of tent table (N(d(z,c)<d(z,a)) because a is always > c
for it=1:nt %it=1 for a, 2 for b, 3 for c
    tcol=[1 1+setdiff([1:nt],it)]; %[1 3 4], [1 2 4], or [1 2 3]
    tent_toks=tents(:,tcol);
    z_lo=find(tent_toks(:,1)<tent_toks(:,2));
    z_mi=find(tent_toks(:,1)>tent_toks(:,2) & tent_toks(:,1)<tent_toks(:,3));
    z_hi=find(tent_toks(:,1)>tent_toks(:,3));
    disp(sprintf(' tcol length(lo,mi,hi): %2.0f %2.0f %2.0f   %8.0f %8.0f %8.0f tot %10.0f',...
        tcol,length(z_lo),length(z_mi),length(z_hi),length(union(union(z_lo,z_mi),z_hi))));
    %z_lo: triad already has z first, i.e., z, tcol(2), tcol(3), so we need
    %first column of triplet table, which is r(z;tcol(2),tcol(3))
    tentz_lo_rows=nchoosek2seq_3vr(tent_toks(z_lo,:),nstims);
    ncloser(z_lo,it)=ncloser_triplets(tentz_lo_rows,1);
    ntrials(z_lo,it)=ntrials_triplets(tentz_lo_rows,1);
    %z_mi: triad has z in the middle, i.e., tcol(2), z, tcol(3), so we need
    %second column of triplet table, which is r(z;tcol(3),tcol(2)), we t invert the choices since the order is an odd permutation in triplets
    tentz_mi_rows=nchoosek2seq_3vr(tent_toks(z_mi,[2 1 3]),nstims);
    ncloser(z_mi,it)=ntrials_triplets(tentz_mi_rows,2)-ncloser_triplets(tentz_mi_rows,2);
    ntrials(z_mi,it)=ntrials_triplets(tentz_mi_rows,2);
    %z_hi: triad has z lset, i.e., tcol(2), tcol(3), z so we need
    %third column of triplet table, which is r(z; tcol(2), tcol(3))
    tentz_hi_rows=nchoosek2seq_3vr(tent_toks(z_hi,[2 3 1]),nstims);
    ncloser(z_hi,it)=ncloser_triplets(tentz_hi_rows,3);
    ntrials(z_hi,it)=ntrials_triplets(tentz_hi_rows,3);
end
ncloser(:,2)=ntrials(:,2)-ncloser(:,2); %invert choices since second column of tents is N(d(z,c)<d(z,a)) and a>c
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
disp('finding Dirichlet and edge parameters for each threshold');
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
%
%now analyze for consistency with addtree
%
r.adt.thr_types=thr_types;
r.adt.columns_tallies={'thr','ntents','ntrials'};
r.adt.columns_a={'a','loglike_per_trial'}; %values of a and h determined from the selected trials
r.adt.columns_ah={'a','h','loglike_per_trial'}; %values of a and h determined from the selected trials
r.adt.columns_adt={'llr_addtree_trans_vs_trans_tent'}; %not excluded by addtree_trans vs not excluded by anything trans
r.adt.columns_adt_hfixed=r.adt.columns_adt; %like adt, but with fixed values of h
r.adt.global.a=r.dirichlet.a(1,1,:); % values with h fixed
%compute using global a and h from unthresholded Dirichlet and save in r.adt.global.ah
if r.dirichlet.ah(1,2)>=0 %use full fit if h>=0
    r.adt.global.ah=r.dirichlet.ah(1,1:2);
else %otherwise use best fit with h=0
    r.adt.global.ah=[r.dirichlet.a(1,1,1),0];
end
%compute later using private a and h, to go in r.adt.private.[a|ah]{ithr_type}
r.adt.private.a=cell(1,nthr_types); 
r.adt.private.ah=cell(1,nthr_types);
%
r.adt.global_private_d1={'mean of sum','variance of sum'};
r.adt.global_private_d2={'threshold type'};
r.adt.global.adt=cell(2,nthr_types); %d1: mean of sum, d2: var of sum
r.adt.global.adt_hfixed=cell(2,nthr_types);
%
r.adt.private.adt=cell(2,nthr_types);
r.adt.private.adt_hfixed=cell(2,nthr_types);
%
r.adt.llr_d1={'threshold value'};
r.adt.llr_d2={'orig data','flip all','flip any'};
r.adt.llr_d3={'hfixed'};
nsurr=length(r.adt.llr_d2); %three kinds of surrogates: native, flip all, flip any
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
    logic_type.adt='exclude_addtree_trans';
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
    r.adt.llr_d2{end+1}=strrep(opts_conform.method,'_',' ');
else
    nconform=0;
end
r.nsurr=nsurr;
r.nconform=nconform;
llr_adt=cell(nsurr+nconform,2); %summed log likelihood ratio across trials, and summed variance of total 
llr_adt_hfixed=cell(nsurr+nconform,2);
llr_adt=cell(nsurr,2); %summed log likelihood ratio across trials, and summed variance of total 
llr_adt_hfixed=cell(nsurr,2);
surr_list={1,[1 nflips],[1:nflips]};
obs_all=[reshape(ncloser',[ncomps 1 ntents]),reshape(ntrials',[ncomps 1 ntents])];
if (if_fast~=0)
    liks_all=zeros(nineq,nflips,ntents);
    liks_hfixed_all=zeros(nineq,nflips,ntents,nhfix);
    %if_fast=1: calculate probabilities for all triplets
    ah=r.adt.global.ah;
    ah_fixed=[squeeze(r.dirichlet.a(1,1,:)),h_fixlist(:)];
    obs_all=[reshape(ncloser',[ncomps 1 ntents]),reshape(ntrials',[ncomps 1 ntents])];
    params.a=ah(1);
    params.h=ah(2);
    liks_all=psg_ineq_apply(params,obs_all,partitions,permutes);
    for ihfix=1:nhfix
        params.a=ah_fixed(ihfix,1);
        params.h=ah_fixed(ihfix,2);
        liks_hfixed_all(:,:,:,ihfix)=psg_ineq_apply(params,obs_all,partitions,permutes);
    end
    disp(sprintf(' global calculations done.'));
end
ipg_strings={'private','global'};
for ipg=ipg_min:2 %private and global
    disp(sprintf('%10s calculations',ipg_strings{ipg}));
    for ithr_type=1:nthr_types %three kinds of thresholds: min, max, average
        if_ok=1;
        thr=0; %threshold
        ithr=1; %threshold pointer
        disp(sprintf('analyzing for addtree likelihood ratio for threshold type %s',thr_types{ithr_type}));
        nuse_prev=-1; %will allow for reuse if increasing the threshold doesn't change the number of triplets/tents used
        while (if_ok)
            switch thr_types{ithr_type}
                case 'min'
                    tents_use=find(min(ntrials,[],2)>=thr);
                    thr_val=thr;
                case 'max'
                    tents_use=find(max(ntrials,[],2)>=thr);
                    thr_val=thr;
                case 'avg'
                    tents_use=find(sum(ntrials,2)>=thr);
                    thr_val=thr/ncomps; %average not total
            end
            if (length(tents_use)>=nfit_min)
                ntents_use=length(tents_use);
                if ntents_use~=nuse_prev
                    did_or_skipped='did'; %have to calculate
                    nuse_prev=ntents_use;
                    ntrials_use=sum(sum(ntrials(tents_use,:)));
                    r.adt.tallies{ithr_type}(ithr,:)=[thr_val ntents_use ntrials_use]; %threshold, number of triplets, number of trials
                    %compute private best-fitting a and h
                    data_use=[reshape(ncloser(tents_use,:),ncomps*ntents_use,1) reshape(ntrials(tents_use,:),ncomps*ntents_use,1)];
                    %fit with assuming fixed values of h
                    for ihfix=1:nhfix
                        if (if_fixa==0)
                            [fit_a,nll_a,exitflag_a]=fminbnd(@(x) -loglik_beta(x,data_use,setfield(opts_loglik,'hvec',h_fixlist(ihfix))),...
                                a_limits(1),a_limits(2)); %optimize assuming discrete part
                        else
                            fit_a=a_fixval;
                            nll_a=-loglik_beta(fit_a,data_use,setfield(opts_loglik,'hvec',h_fixlist(ihfix)));
                        end
                        r.adt.private.a{ithr_type}(ithr,:,ihfix)=[fit_a,-nll_a/ntrials_use];
                    end
                    ah_init=[r.adt.private.a{ithr_type}(ithr,1,1);h_init]; %optimize with discrete part, using a_only fit as starting point
                    [fit_ah,nll_ah,exitflag_ah,output_ah]=fminsearch(@(x) -loglik_beta(x(1),data_use,setfield(opts_loglik,'hvec',x(2))),ah_init);
                    if fit_ah(2)>=0
                        r.adt.private.ah{ithr_type}(ithr,:)=[fit_ah(:)',-nll_ah/ntrials_use];
                    else
                        r.adt.private.ah{ithr_type}(ithr,:)=[fit_a,0,-nll_a/ntrials_use];
                    end
                    liks=zeros(nineq,nflips,ntents_use);
                    liks_hfixed=zeros(nineq,nflips,ntents_use,nhfix);
                    %
                    %fast global option:calculate probabilities for all triplets and later select
                    %
                    if if_fast~=0 & ipg==2
                        liks=liks_all(:,:,tents_use); %d1: addtree_trans vs trans_tent, d2: flips, d3: tent
                        liks_hfixed=liks_hfixed_all(:,:,tents_use,:); %d1: addtree_trans vs trans_tent, d2: flips, d3: tent, d4:h_fixed
                    else %if_fast==0
                        if (ipg==1) %private
                            ah=r.adt.private.ah{ithr_type}(ithr,:); %a and h both fitted
                            ah_fixed=[squeeze(r.adt.private.a{ithr_type}(ithr,1,:)),h_fixlist(:)]; %a fitted, h fixed
                        else %global
                            ah=r.adt.global.ah;
                            ah_fixed=[squeeze(r.dirichlet.a(1,1,:)),h_fixlist(:)];
                        end
                        params.a=ah(1);
                        params.h=ah(2);
                        liks=psg_ineq_apply(params,obs_all(:,:,tents_use),partitions,permutes);
                        for ihfix=1:nhfix
                            params.a=ah_fixed(ihfix,1);
                            params.h=ah_fixed(ihfix,2);
                            liks_hfixed(:,:,:,ihfix)=psg_ineq_apply(params,obs_all(:,:,tents_use),partitions,permutes);
                        end
                    end %if_fast
                    %do statistics
                    %likelihood of addtree and trans,/likelihood(trans), i.e., exclude_addtree_trans/exclude_trans_tent';
                    %dimensions reordered to match those of loglik_rat_sym|umi[|_hfixed] of psg_umi_triplike_demo
                    loglikrats=transpose(log(reshape(liks(1,:,:)./liks(2,:,:),[nflips ntents_use]))); %after transpose: d1: tents, d2: flips
                    loglikrats_hfixed=permute(log(reshape(liks_hfixed(1,:,:,:)./liks_hfixed(2,:,:,:),[nflips ntents_use nhfix])),[2 1 3]); %d1:tents, d2:flips, d3:h
                    for isurr=1:nsurr+nconform %for each kind of surrogate
                        if (isurr<=nsurr)
                            surr_sel=surr_list{isurr}; %for isurr=1, this is just the original data (1)
                            llr_adt{isurr,1}=sum(mean(loglikrats(:,surr_sel),2),1);
                            llr_adt_hfixed{isurr,1}=reshape(sum(mean(loglikrats_hfixed(:,surr_sel,:),2),1),[1 1 nhfix]);
                            if (isurr>1)
                                %each tent contributes independently to the variance
                                %variance for each tent is normalized by N not N-1, since we have all the values
                                llr_adt{isurr,2}=sum(var(loglikrats(:,surr_sel),1,2),1);
                                llr_adt_hfixed{isurr,2}=reshape(sum(var(loglikrats_hfixed(:,surr_sel,:),1,2),1),[1 1 nhfix]);
                            else %isurr=1: original data. Here, goal is for psg_umi_triplike_plota to compute standard error of the mean
                                %which is sqrt(var)/ntents_use, but
                                %psg_umi_triplike_plota will find square root and then divide by ntents_use
                                %so here we just compute var, normalized by N-1 since it is a sample
                                %here, surr_sel=1
                                llr_adt{isurr,2}=var(loglikrats(:,surr_sel),0,1);
                                llr_adt_hfixed{isurr,2}=reshape(var(loglikrats_hfixed(:,surr_sel,:),0,1),[1 1 nhfix]);
                            end
                        else %do conform
                            %select the appropriate flip for each triplet
                            loglik_rat_adt_conform=zeros(ntents_use,1);
                            loglik_rat_adt_hfixed_conform=zeros(ntents_use,nhfix);
                            for itptr=1:ntents_use
                                it=tents_use(itptr);
                                %logic differs from psg_umi_triplike_demo since we use loglikrats[_hfixed]
                                loglik_rat_adt_conform(itptr)=loglikrats(itptr,which_flip_conform.adt(it));
                                loglik_rat_adt_hfixed_conform(itptr,:)=reshape(loglikrats_hfixed(itptr,which_flip_conform.adt(it),:),[1 nhfix]);
                            end
                            llr_adt{isurr,1}=sum(loglik_rat_adt_conform);
                            llr_adt_hfixed{isurr,1}=reshape(sum(loglik_rat_adt_hfixed_conform),[1 1 nhfix]);
                            %variances are calculated as for original data, but surr_sel must be set
                            surr_sel=1;
                            llr_adt{isurr,2}=var(loglikrats(:,surr_sel),0,1);
                            llr_adt_hfixed{isurr,2}=reshape(var(loglikrats_hfixed(:,surr_sel,:),0,1),[1 1 nhfix]);
                        end
                        for imv=1:2% mean and variance
                            r.adt.(ipg_strings{ipg}).adt{imv,ithr_type}(ithr,isurr)=llr_adt{isurr,imv};
                            r.adt.(ipg_strings{ipg}).adt_hfixed{imv,ithr_type}(ithr,isurr,:)=llr_adt_hfixed{isurr,imv};
                        end %imv
                     end %isurr
                else
                    did_or_skipped='skp';
                    r.adt.tallies{ithr_type}(ithr,:)=r.adt.tallies{ithr_type}(ithr-1,:);
                    r.adt.tallies{ithr_type}(ithr,1)=thr_val; %threshold is new
                    if (ipg==1)
                        r.adt.private.a{ithr_type}(ithr,:,:)=r.adt.private.a{ithr_type}(ithr-1,:,:);
                        r.adt.private.ah{ithr_type}(ithr,:)=r.adt.private.ah{ithr_type}(ithr-1,:);
                    end
                    for isurr=1:nsurr+nconform %for each kind of surrogate
                        for imv=1:2% mean and variance
                            r.adt.(ipg_strings{ipg}).adt{imv,ithr_type}(ithr,isurr)=llr_adt{isurr,imv};
                            r.adt.(ipg_strings{ipg}).adt_hfixed{imv,ithr_type}(ithr,isurr,:)=llr_adt_hfixed{isurr,imv};
                        end %imv
                     end %isurr
                end %nuse_prev
                disp(sprintf('%s ipg %3.0f ithr_type %3.0f ithr %3.0f thr %3.0f ntents_use %6.0f size(liks) %4.0f %4.0f %6.0f size(liks_hfixed) %4.0f %4.0f %6.0f %4.0f',...
                    did_or_skipped,ipg,ithr_type,ithr,thr,ntents_use,size(liks),size(liks_hfixed)));
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
plot_opts.llr_field='adt';
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

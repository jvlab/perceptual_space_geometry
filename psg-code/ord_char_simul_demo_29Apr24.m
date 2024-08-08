%ord_char_simul_demo: simulate analysis of rank order for toy models
%
% this is in conjunction with ord_char_sim_*.doc manuscript, and  
% psg_umi_triplike*.m (for symmetry and ultrametric analysis)
% psg_tentlike*.m (for addtree analysis)
%
% There is a difference in the way statisticss are accumulated here vs. in psg_[umi_trip|tent]like_demo:
% because here we have an exhaustive sampling of all triplets and triads
% 
% llr_[sym|umi|adt] is the sum of the log likelihood ratios across all triplets or tents
%    and for surrogates, surrogates are averaged within each triplet or tent
% llr_vm_[sym|umi|adt] is the variance of the log likelihood ratios across all triplets or tents
%    and for surrogates, surrogates are averaged within each triplet or
%    tent; this can be used for computing the effects of sampling a random
%    subset of triplets or tents
% llr_sv_[sym|umi|adt] is the sum of the variances within surrogates, so it
%    reflects the the distribution of surrogate values (and 0 for the real data)
% llr_vv_[sym|umi|adt] is the variance of the variances within surrogates, so it reflects
%    the sensitivity of llr_sv to sampling a random subset of triplets or tents
%
%  All of these are computed with a denominator of N, not N-1, since the
%  samplings here are complete.
%
%  In contrast, in psg_[umi_trip|tent]like_demo], the variances for the
%  real data correspond to llr_vm above, computed with N-1, and the variances for the
%  surrogate data correspond to llr_sv above
%
% With Poisson-distributed trial counts, triplets and tents with no trials are not considered
% 
%   Uses graph toolbox.
%
% See also: ORD_CHAR_SIMUL_DIST, ORD_CHAR_SIMUL_PLOT, PSG_TRIPLET_CHOICES, PSG_TENT_CHOICES,
% PSG_UMI_TRIPLIKE, PSG_TENTLIKE, PSG_INEQ_LOGIC, PSG_PERMUTES_LOGIC.
%
if ~exist('opts_triplike')
    opts_triplike=struct;
    opts_triplike.if_fast=[];
end
%
gsets={'tree','star','ladder','grid'}; %each set of geometries has a common graph and same number of points
%
if getinp('1 to use debug params','d',[0 1])
    if ~exist('nlevels') nlevels=3; end
    if ~exist('narms') narms=4; end
    if ~exist('npts_per_arm') npts_per_arm=2; end
    if ~exist('nrungs') nrungs=4; end
    if ~exist('ngridsize') ngridsize=3; end
    if ~exist('h_fixlist') h_fixlist=[0 0.001 0.1]; end %fixed nonzero values of h for fitting
    if ~exist('a_fixlist') a_fixlist=[0.125 1]; end %differs from psg_[umi_trip|tent]like_demo, want to step through many specific values
    if ~exist('rules') rules=[0.125 0; 0 0.125; 0 1]; end %rules(:,1) is button error, rules(:,2) is sigma
    if ~exist('trials_list') trials_list=[1 2 8 16]; end
end
if ~exist('nlevels') nlevels=4; end
if ~exist('narms') narms=4; end
if ~exist('npts_per_arm') npts_per_arm=6; end
if ~exist('nrungs') nrungs=8; end
if ~exist('ngridsize') ngridsize=4; end
%
if ~exist('h_fixlist') h_fixlist=[0 0.001 0.01 0.0625, 0.125 0.25 0.5]; end %fixed nonzero values of h for fitting
if ~exist('a_fixlist') a_fixlist=[0.125 0.25 0.5 1 2]; end %differs from psg_[umi_trip|tent]like_demo, want to step through many specific values
if ~exist('rules') rules=[0.25 0; 0.125 0; 0.0625 0; 0.03125 0; 0 0.0625; 0 0.125; 0 0.25; 0 0.5; 0 1;0 2]; end %rules(:,1) is button error, rules(:,2) is sigma
if ~exist('trials_list') trials_list=[1 2 4 8 16 32 64]; end
%
if ~exist('if_normdist') if_normdist=1; end %set to normalize distances
if ~exist('nbins_cp') nbins_cp=21; end %bins for choice probabilities
bincenters_cp=[0.5:nbins_cp-0.5]*(1/(nbins_cp-1));
%
%for fitting dirichlet params, adapted from psg_[umi_trip|tent]like_demo
if ~exist('a_limits') a_limits=[2^-10 2^3]; end
if ~exist('h_limits') h_limits=[0 1]; end
if ~exist('h_init') h_init=0; end %initial value for h
h_fixlist=unique([0 h_fixlist(:)']);
nhfix=length(h_fixlist);
nafix=length(a_fixlist);
%
opts_loglik=struct;
opts_loglik.qvec=0.5; %discrete part always is at 0.5
%
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
if (if_frozen~=0) %also need to repeat this below prior to each simulation of choice probabilities
    rng('default');
    if (if_frozen<0)
        rand(1,abs(if_frozen));
    end
else
    rng('shuffle');
end
if_poisson=getinp('1 to use Poisson distribution for trial-repeat counts','d',[0 1],0);
%
nrules=size(rules,1);
ntrials_list=length(trials_list);
%
%choose a set of geometries
%
for ig=1:length(gsets)
    disp(sprintf('%1.0f->%s',ig,gsets{ig}));
end
ig=getinp('choice','d',[1 length(gsets)],1);
gset=gsets{ig};
%
S=struct;
S.coords=zeros(0,2);
S.edges=zeros(0,2);
S.ranks=zeros(0,1);
S.edgelengths=zeros(0,1);
switch gset
    case 'tree' %a binary tree with nleels levels
        nlevels=getinp('number of levels in the tree','d',[1 6],nlevels); %number of levels in a binary tree toy model
        lab_string=sprintf('%s, %2.0f levels',gset,nlevels);
        dist_types_nr=cell(0);
        %create graph and coords of points in a binary tree with nlevels
        for ilev=1:nlevels
            npts=2^(nlevels-ilev);
            spacing=2^(ilev-1);
            posits=spacing*([0:npts-1]-(npts-1)/2);
            if (ilev>1)
                pointnos=size(S.coords,1)+[1:npts];
                points_prev=size(S.coords,1)-npts*2+[1:npts*2];
                edges_prev=size(S.edges,1);
                S.edges=[S.edges;[pointnos(:),points_prev(1:2:end)']];
                S.edges=[S.edges;[pointnos(:),points_prev(2:2:end)']];
            end
            S.coords=[S.coords;[posits(:),repmat(ilev,npts,1)]];
            S.ranks=[S.ranks;repmat(ilev,npts,1)];
        end
    case 'star'
        narms=getinp('number of arms','d',[1 8],narms); %number of arms
        npts_per_arm=getinp('number of points on each arm','d',[1 12],npts_per_arm); %points on each arm
        vertex_angs=getinp('angle (deg)','f',[0 360],[0:narms-1]*360/narms);
        vertex_angs=vertex_angs(1+mod([0:narms-1],length(vertex_angs)));
        lab_string_vtx=sprintf('%5.2f ',vertex_angs);
        lab_string=sprintf('%s, %2.0f arms, %2.0f pts each, angles %s deg',gset,narms,npts_per_arm,lab_string_vtx);
        angs=vertex_angs*(pi/180);
        dist_types_nr={'graph_weighted','line','ring_as_graph'};
        npts=1+narms*npts_per_arm;
        S.coords=zeros(npts,2);
        S.ranks=npts_per_arm+2-[1 repmat(2:npts_per_arm+1,1,narms)]'; %point at origin has highest rank as a tree
        S.edges=zeros(npts-1,2);
        %create graph and coords of points for a star with one point at center, and narms with npts_per_arms on each arm
        for ilev=2:npts_per_arm+1
            S.coords(1+narms*(ilev-2)+[1:narms],:)=(ilev-1)*[cos(angs'),sin(angs')];
            S.edges(narms*(ilev-2)+[1:narms],:)=[1+narms*(ilev-2)+[1:narms]',max(1,1+narms*(ilev-3)+[1:narms]')];
        end
    case {'ladder','grid'}
        if strcmp(gset,'ladder')
            nx=getinp('number of rungs','d',[2 32],nrungs);
            ny=2;
            lab_string=sprintf('%s, %2.0f rungs',gset,nrungs);
        else
            nx=getinp('grid size','d',[2 6],ngridsize);
            ny=nx;
            lab_string=sprintf('%s, %2.0f x %2.0f',gset,nx,nx);
        end
        dist_types_nr={'ultra','graph_weighted','line','ring_as_graph'};
        npts=nx*ny;
        S.ranks=ones(npts,1); %ranks are trivial
        S.edges=zeros(0,2);
        [coords_x,coords_y]=meshgrid([0:nx-1],[0:ny-1]);
        S.coords=[coords_x(:),coords_y(:)];
        for ix=1:nx %edges along x axis
            S.edges=[S.edges;[1:ny-1]'+(ix-1)*ny,[2:ny]'+(ix-1)*ny];
        end
        for iy=1:ny %edges along y axis
            S.edges=[S.edges;iy+ny*[0:nx-2]',iy+ny*[1:nx-1]'];
        end
end
%common calculations for all graphs
S.nedges=size(S.edges,1);
S.npts=size(S.coords,1);
% S.edgelengths=sqrt(...
%     (S.coords(S.edges(:,1),1)-S.coords(S.edges(:,2),1)).^2+...
%     (S.coords(S.edges(:,1),2)-S.coords(S.edges(:,2),2)).^2); %should have same length within rank
S.edgelengths=sqrt(diff(reshape(S.coords(S.edges,1),S.nedges,2),[],2).^2+...
    diff(reshape(S.coords(S.edges,2),S.nedges,2),[],2).^2); %should have same length within rank
S.graph=graph();
S.graph=addnode(S.graph,S.npts);
S.graph=addedge(S.graph,S.edges(:,1),S.edges(:,2),S.edgelengths);
%
%show the domain and connect edges
%
figure;
set(gcf,'Position',[100,100,1200,800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name',lab_string);
subplot(1,3,1);
plot(S.coords(:,1)-mean(S.coords(:,1)),S.coords(:,2)-mean(S.coords(:,2)),'k.','MarkerSize',10,'Color','k');
hold on;
for iedge=1:size(S.edges,1)
    plot(S.coords(S.edges(iedge,:),1)-mean(S.coords(:,1)),S.coords(S.edges(iedge,:),2)-mean(S.coords(:,2)),'k','LineWidth',1);
end
axis equal;
set(gca,'XLim',(max(abs(S.coords(:,1)-mean(S.coords(:,1))))+0.5)*[-1 1]);
set(gca,'YLim',(max(abs(S.coords(:,2)-mean(S.coords(:,2))))+0.5)*[-1 1]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
title('points as a graph');
%
subplot(1,3,2);
plot([1:S.npts],zeros(1,S.npts),'k.','MarkerSize',10,'Color','k');
hold on;
plot([1 S.npts],[0 0],'k','LineWidth',1);
axis equal;
set(gca,'XLim',[0.5 S.npts+0.5]);
set(gca,'YLim',[-1.5 1.5]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
title('points as a line');
%
subplot(1,3,3);
rad=S.npts/(2*pi);
angle=[0:S.npts-1]*2*pi/S.npts;
plot(rad*cos(angle),rad*sin(angle),'k.','MarkerSize',10,'Color','k');
hold on;
plot(rad*cos(angle([1:end 1])),rad*sin(angle([1:end 1])),'k','LineWidth',1);
axis equal;
set(gca,'XLim',(rad+1)*[-1 1]);
set(gca,'YLim',(rad+1)*[-1 1]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
title('points as a ring');
%
axes('Position',[0.01,0.02,0.01,0.01]); %for text
text(0,0,lab_string,'Interpreter','none','FontSize',10);
axis off;
%
%compute and show all pairwise distances
%
dist_types_all=ord_char_simul_dist; %all distance types are possible
for itype=1:length(dist_types_all)
    if ~isempty(strmatch(dist_types_all{itype},dist_types_nr,'exact'))
        disp(sprintf('distance type %2.0f->%20s [redundant or not recommended]',itype,dist_types_all{itype}));
    else
        disp(sprintf('distance type %2.0f->%20s',itype,dist_types_all{itype}));
    end
end
%
%choose from library of distances
%
dist_list=getinp('choice(s)','d',[1 length(dist_types_all)],[1:length(dist_types_all)]);
dist_types=cell(1,length(dist_list));
ndist_types=length(dist_list);
for idist=1:ndist_types
    dist_types{idist}=dist_types_all{dist_list(idist)};
end
dists=struct;
%
for itype=1:ndist_types
    dist_type=dist_types{itype};
    dists.(dist_type)=zeros(S.npts);
    for ip1=1:S.npts
        for ip2=1:S.npts
            dists.(dist_type)(ip1,ip2)=ord_char_simul_dist(dist_type,[ip1,ip2],S);
        end
    end
    %normalize distances to have a mean of 1 if requested
    if if_normdist
       dist_mean=sum(dists.(dist_type)(:))/S.npts/(S.npts-1); %don't count self-distances
       dists.(dist_type)=dists.(dist_type)/dist_mean;
    end
end
%show distances
figure;
set(gcf,'Position',[100,100,800,800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name','distances');
[nr,nc]=nicesubp(ndist_types);
for itype=1:ndist_types
    dist_type=dist_types{itype};
    subplot(nr,nc,itype);
    imagesc(dists.(dist_type),[0 max(1,max(dists.(dist_type)(:)))]);
    set(gca,'XTick',[1:S.npts]);
    set(gca,'YTick',get(gca,'XTick'));
    title(sprintf('%s, if_normdist=%2.0f',dist_type,if_normdist),'Interpreter','none');
    axis tight;
    axis equal;
    colorbar;
end
axes('Position',[0.01,0.02,0.01,0.01]); %for text
text(0,0,lab_string,'Interpreter','none','FontSize',10);
axis off;
%
%setups for symmetry and umi
%
ncomps_su=3; %number of comparisons for symmetry and umi
flipconfigs_su=int2nary([0:2^ncomps_su-1]',2);  %rows are [0 0 0;1 0 0;0 1 0;1 1 0; 0 0 1;1 0 1;0 1 1;1 1 1];
nflips_su=size(flipconfigs_su,1); 
surr_list_su={1,[1 nflips_su],[1:nflips_su]};
%
%setups for addtree
ineq_logic_types_adt={'exclude_addtree_trans','exclude_trans_tent'};
nineq_adt=length(ineq_logic_types_adt);
%
ncomps_adt=6; %six rank choice probabilities to be compared
%
partitions_adt=cell(0);
for ineq=1:nineq_adt
    disp(sprintf('creating inequality logic for %s',ineq_logic_types_adt{ineq}));
    partitions_adt{ineq}=psg_ineq_logic(ncomps_adt,ineq_logic_types_adt{ineq},setfield([],'if_log',1));
end
disp('creating permutation logic for surrogates')
permutes_adt=psg_permutes_logic(ncomps_adt,'flip_each');
disp(sprintf(' size is %3.0f x %3.0f',size(permutes_adt)));
nflips_adt=size(permutes_adt,2);
%no need for flipconfigs_adt, this functoin is performed by permutes_adt
surr_list_adt={1,[1 nflips_adt],1:nflips_adt};
%
rmeta=struct;
rmeta.llr_d1={'ntrials_ptr'};
rmeta.llr_d2={'orig data','flip all','flip any'};
rmeta.llr_d3={'a index'};
rmeta.llr_d4={'h index'};
rmeta.llr={'llr: log likelihood ratios summed across trials and triplets or tents'};
rmeta.llr_vm={'llr_vm: variances of mean log likelihood ratios across surrogates, within triplets or tents'};
rmeta.llr_sv={'llr_vm: sum of variances of log likelihood ratios across surrogates, within triplets or tents'};
rmeta.llr_vv={'llr_vm: variances of variances of log likelihood ratios across surrogates, within triplets or tents'};
obs_orig_su=zeros(ncomps_su,2);
%
nsurr=length(rmeta.llr_d2); %three surrogates (1: original data, 2 flip all, 3: flip any)
%
%step through all decision rules and geometries
%
results=cell(ndist_types,nrules);
ntriplets=nchoosek(S.npts,3); %exhaustively sample the triplets
ntriads=3*ntriplets;
ntents=4*nchoosek(S.npts,4); %exhaustively sample tents
for irule=1:nrules
    disp(' ')
    button=rules(irule,1);
    sigma=rules(irule,2);
    rule_string=sprintf('button error %5.3f sigma %5.3f',button,sigma);
    rule_string_nice=cat(2,'\',sprintf('epsilon_b_u_t_t_o_n=%5.3f ',button),'\',sprintf('sigma=%4.2f',sigma));
    for itype=1:ndist_types
        dist_type=dist_types{itype};
        disp(sprintf(' decision rule: %s    distance type: %s',rule_string,dist_type));
        tic;
        %compute true probabilities for every possible triad
        cp=zeros(S.npts,S.npts,S.npts);
        cp_table=zeros(ntriads,4); %ir, ix, iy, cp
        triad_ptr=0;
        for ir=1:S.npts
            if (if_frozen~=0) %reset random number generator so that simulation does not depend on which other conditions were run
                rng('default');
                if (if_frozen<0)
                    rand(1,abs(if_frozen));
                end
            end
            cp(ir,ir,:)=NaN;
            cp(ir,:,ir)=NaN;
            for ix=2:S.npts
                for iy=1:ix-1
                    if (ix~=ir) & (iy~=ir)
                        dx=ord_char_simul_dist(dist_type,[ir,ix],S);
                        dy=ord_char_simul_dist(dist_type,[ir,iy],S);
                        %compute choice probability bias (1 if x always judged closer to r than y, -1 if x always judged further) (without button error)
                        if sigma==0
                            cpb=sign(dy-dx);
                        else
                            cpb=erf((dy-dx)/(2*sigma)); %as in Waraich JN 2023, 2 sigma in denom
                        end
                        %convert choice prob bias to choice prob, based on button error
                        cp_xy=0.5*(1+cpb*(1-2*button));
                        cp(ir,ix,iy)=cp_xy;
                        cp(ir,iy,ix)=1-cp_xy;
                        %tabulate
                        triad_ptr=triad_ptr+1;
                        cp_table(triad_ptr,:)=[ir ix iy cp(ir,ix,iy)];
                    end
                end %iy
            end %ix
        end %ir
        r=rmeta;
        r.geometry_set=lab_string;
        r.cp_table=cp_table;
        r.button=button;
        r.sigma=sigma;
        r.dist_type=dist_type;
        r.if_normdist=if_normdist;
        r.rule_string=rule_string;
        r.rule_string_nice=rule_string_nice;
        r.trials_list=trials_list;
        r.trials_if_poisson=if_poisson;
        r.a_fixlist=a_fixlist;
        r.h_fixlist=h_fixlist;
        r.dirichlet=struct();
        r.dirichlet_dims={'a_ptr (0 for fit)','h_ptr (0 for fit)','itrial_ptr'};
        r.dirichlet.a=NaN(nafix+1,nhfix+1,ntrials_list);
        r.dirichlet.h=NaN(nafix+1,nhfix+1,ntrials_list);
        r.dirichlet.nll=NaN(nafix+1,nhfix+1,ntrials_list); %negative log likelihoods
        r.C=cell(1,ntrials_list);
        %sums of log likelihood raios
        r.llr_sym=NaN(ntrials_list,3,nafix+1,nhfix+1); %d1={'ntrials_ptr'}; d2={'orig data','flip all','flip any'};d3={'a index'}; rmeta.llr_d4={'h index'};
        r.llr_umi=NaN(ntrials_list,3,nafix+1,nhfix+1); %d1={'ntrials_ptr'}; d2={'orig data','flip all','flip any'};d3={'a index'}; rmeta.llr_d4={'h index'};
        r.llr_adt=NaN(ntrials_list,3,nafix+1,nhfix+1); %d1={'ntrials_ptr'}; d2={'orig data','flip all','flip any'};d3={'a index'}; rmeta.llr_d4={'h index'};
        % variances of surrogate means 
        r.llr_vm_sym=r.llr_sym;
        r.llr_vm_umi=r.llr_umi; 
        r.llr_vm_adt=r.llr_adt;
        % sums of surrogate variances
        r.llr_sv_sym=r.llr_sym;
        r.llr_sv_umi=r.llr_umi;
        r.llr_sv_adt=r.llr_adt;
        % variances of surrogate variances
        r.llr_vv_sym=r.llr_sym;
        r.llr_vv_umi=r.llr_umi;
        r.llr_vv_adt=r.llr_adt;
        %
        r.triplet_counts=zeros(1,ntrials_list);
        r.tent_counts=zeros(1,ntrials_list);
        %
        for itrial_ptr=1:ntrials_list
            ntrials=trials_list(itrial_ptr);
            if (if_poisson)
                ntrials_each=poissrnd(repmat(ntrials,ntriads,1));
            else
                ntrials_each=repmat(ntrials,ntriads,1);
            end
            %
            %simulate all trials
            %
            C=binornd(ntrials_each,cp_table(:,4));
            r.C{itrial_ptr}=C;
            %create triplet list
            triad_data=[cp_table(:,1:3),C,ntrials_each];
            triad_data=triad_data(ntrials_each>0,:); %remove triads that have not been shown
            [ncloser_triplets,ntrials_triplets,abc_list]=psg_triplet_choices(S.npts,triad_data);
            %create tent list
            [ncloser_tents,ntrials_tents]=psg_tent_choices(S.npts,triad_data,ncloser_triplets,ntrials_triplets);
            %remove all triplets and tents that have no trials
            triplets_shown=find(sum(ntrials_triplets,2)>0);
            ncloser_triplets=ncloser_triplets(triplets_shown,:);
            ntrials_triplets=ntrials_triplets(triplets_shown,:);
            tents_shown=find(sum(ntrials_tents,2)>0);
            ncloser_tents=ncloser_tents(tents_shown,:);
            ntrials_tents=ntrials_tents(tents_shown,:);
            %if not poisson, check counts
            if ~if_poisson
                if ~all(ntrials_triplets(:)==ntrials)
                    disp(sprintf('warning: mismatch in triplet trial count, ntrials=%4.0f',ntrials));
                end
                if size(ncloser_triplets,1)~=ntriplets
                    disp(sprintf('warning: ntriplets mismatch, expecting %4.0f found %4.0f',ntriplets,size(ncloser_triplets,1)));
                end
                if ~all(ntrials_tents(:)==ntrials)
                    disp(sprintf('warning: mismatch in tent    trial count, ntrials=%4.0f',ntrials));
                end
                if size(ncloser_tents,1)~=ntents
                    disp(sprintf('warning: ntents mismatch, expecting %4.0f found %4.0f',ntents,size(ncloser_tents,1)));
                end
            end
            %the number of triplets or tents may not be ntriplets or ntents, if some triads are never included
            r.triplet_counts(1,itrial_ptr)=size(ncloser_triplets,1);
            r.tent_counts(1,itrial_ptr)=size(ncloser_tents,1);
            %
            data_fitprior=[ncloser_triplets(:),ntrials_triplets(:)]; %use all data
            for ia=0:nafix
                for ih=0:nhfix
                    have_ah=0;
                    if (ia>0)
                        if (ih>0)
                            have_ah=1;
                            r.dirichlet.a(ia+1,ih+1,itrial_ptr)=a_fixlist(ia);
                            r.dirichlet.h(ia+1,ih+1,itrial_ptr)=h_fixlist(ih);
                        elseif ntrials>1 %fit h, keep a fixed
                             [fit_h,nll_h,exitflag_h]=fminbnd(@(x) -loglik_beta(a_fixlist(ia),data_fitprior,...
                                 setfield(opts_loglik,'hvec',x)),h_limits(1),h_limits(2)); %optimize assuming a=a_fixlist(ia)
                             have_ah=1;
                             r.dirichlet.a(ia+1,ih+1,itrial_ptr)=a_fixlist(ia);
                             r.dirichlet.h(ia+1,1,itrial_ptr)=fit_h;
                        end
                    elseif ntrials>1
                        if (ih>0) %fit a, keep h fixed
                            [fit_a,nll_a,exitflag_a]=fminbnd(@(x) -loglik_beta(x,data_fitprior,...
                                setfield(opts_loglik,'hvec',h_fixlist(ih))),a_limits(1),a_limits(2)); %optimize assuming discrete part = h_fixlist(ih)
                            have_ah=1;
                            r.dirichlet.a(1,ih+1,itrial_ptr)=fit_a;
                            r.dirichlet.h(1,ih+1,itrial_ptr)=h_fixlist(ih);
                        elseif ntrials>2 %fit a, fit h
                            %find best a for h=0
                            [a_init,nll_a,exitflag_a]=fminbnd(@(x) -loglik_beta(x,data_fitprior,...
                                setfield(opts_loglik,'hvec',0)),a_limits(1),a_limits(2)); %optimize assuming discrete part=0
                                ah_init=[a_init;h_init];
                                [fit_ah,nll_ah,exitflag_ah,output_ah]=fminsearch(@(x) -loglik_beta(x(1),data_fitprior,setfield(opts_loglik,'hvec',x(2))),ah_init);
                                have_ah=1;
                                if fit_ah(2)>0
                                    r.dirichlet.a(1,1,itrial_ptr)=fit_ah(1);
                                    r.dirichlet.h(1,1,itrial_ptr)=fit_ah(2);
                                else %cannot have h<0
                                    r.dirichlet.a(1,1,itrial_ptr)=a_init;
                                    r.dirichlet.h(1,1,itrial_ptr)=0;
                                end
                        end
                    end
                    if have_ah
                        r.dirichlet.nll(ia+1,ih+1,itrial_ptr)=-loglik_beta(...
                           r.dirichlet.a(ia+1,ih+1,itrial_ptr),data_fitprior,...
                           setfield(opts_loglik,'hvec',r.dirichlet.h(ia+1,ih+1,itrial_ptr)));
                    end
                end %ih
            end %ia
            obs_adt=[reshape(ncloser_tents',[ncomps_adt 1 size(ncloser_tents,1)]),reshape(ntrials_tents',[ncomps_adt 1 size(ncloser_tents,1)])]; %setup for tents
            for ia=0:nafix
                for ih=0:nhfix
                    llrs_sym=NaN(length(triplets_shown),nflips_su);
                    llrs_umi=NaN(length(triplets_shown),nflips_su);
                    if ~isnan(r.dirichlet.nll(ia+1,ih+1,itrial_ptr))
                        params_ah=struct;
                        params_ah.a=r.dirichlet.a(ia+1,ih+1,itrial_ptr);
                        params_ah.h=r.dirichlet.h(ia+1,ih+1,itrial_ptr);
                        %
                        %symmetry and umi calculation, following logic of psg_umi_triplike_demo
                        %
                        for itriplet=1:length(triplets_shown)
                            obs_orig_su(:,1)=ncloser_triplets(itriplet,:)';
                            obs_orig_su(:,2)=ntrials_triplets(itriplet,:)';
                            obs_orig_flip_su=obs_orig_su(:,2)-obs_orig_su(:,1);
                            for iflip=1:nflips_su
                                obs_su=obs_orig_su;
                                whichflip=find(flipconfigs_su(iflip,:)==1);
                                obs_su(whichflip,1)=obs_orig_flip_su(whichflip);
                                likrat=psg_umi_triplike(params_ah,obs_su,opts_triplike); %calculate the likelihood ratios
                                llrs_sym(itriplet,iflip)=log(likrat.sym);
                                llrs_umi(itriplet,iflip)=log(likrat.umi_trans);
                            end
                        end %itriplet
                        %
                        %addtree
                        %
                        liks_adt=psg_ineq_apply(params_ah,obs_adt,partitions_adt,permutes_adt); % size is [ni nflips nsets]
                        llrs_adt=transpose(log(reshape(liks_adt(1,:,:)./liks_adt(2,:,:),[nflips_adt size(ncloser_tents,1)]))); %after transpose: d1: tents, d2: flips
                        %
                        %accumulate statistics
                        %  note that size of d3 in llrs_[sym|umi] and llrs_adt differ
                        %  and length of entries in surr_list_su and surr_list_adt also differ
                        %
                        for isurr=1:nsurr
                            r.llr_sym(itrial_ptr,isurr,ia+1,ih+1)=sum(mean(llrs_sym(:,surr_list_su{isurr}),2),1);
                            r.llr_umi(itrial_ptr,isurr,ia+1,ih+1)=sum(mean(llrs_umi(:,surr_list_su{isurr}),2),1);
                            r.llr_adt(itrial_ptr,isurr,ia+1,ih+1)=sum(mean(llrs_adt(:,surr_list_adt{isurr}),2),1);
                            %variances of means
                            r.llr_vm_sym(itrial_ptr,isurr,ia+1,ih+1)=var(mean(llrs_sym(:,surr_list_su{isurr}),2),1,1);
                            r.llr_vm_umi(itrial_ptr,isurr,ia+1,ih+1)=var(mean(llrs_umi(:,surr_list_su{isurr}),2),1,1);
                            r.llr_vm_adt(itrial_ptr,isurr,ia+1,ih+1)=var(mean(llrs_adt(:,surr_list_adt{isurr}),2),1,1);
                            %sum of variances
                            r.llr_sv_sym(itrial_ptr,isurr,ia+1,ih+1)=sum(var(llrs_sym(:,surr_list_su{isurr}),1,2),1);
                            r.llr_sv_umi(itrial_ptr,isurr,ia+1,ih+1)=sum(var(llrs_umi(:,surr_list_su{isurr}),1,2),1);
                            r.llr_sv_adt(itrial_ptr,isurr,ia+1,ih+1)=sum(var(llrs_adt(:,surr_list_adt{isurr}),1,2),1);
                            %variances of variances
                            r.llr_vv_sym(itrial_ptr,isurr,ia+1,ih+1)=var(var(llrs_sym(:,surr_list_su{isurr}),1,2),1,1);
                            r.llr_vv_umi(itrial_ptr,isurr,ia+1,ih+1)=var(var(llrs_umi(:,surr_list_su{isurr}),1,2),1,1);
                            r.llr_vv_adt(itrial_ptr,isurr,ia+1,ih+1)=var(var(llrs_adt(:,surr_list_adt{isurr}),1,2),1,1);
                        end %isurr
                    end %ia, ih available
                end %ih
            end %ia
        end %itrial_ptr
        results{itype,irule}=r;
        disp(sprintf('   time: %8.4f s',toc));
    end %itype
end %irule
%show underlying and empiric choice probabilities
for itrial_ptr=0:ntrials_list
    figure; %choice probabilities
    set(gcf,'Position',[100 100 1200 800]);
    set(gcf,'Numbertitle','off');
    if (itrial_ptr>0)
        name_string=sprintf('empiric choice probabilities, %4.2f trials, if_poisson=%1.0f',trials_list(itrial_ptr),if_poisson);
    else
        name_string='choice probabilities';
    end
    set(gcf,'Name',name_string);
    for irule=1:nrules
        button=rules(irule,1);
        sigma=rules(irule,2);
        for itype=1:ndist_types
            dist_type=dist_types{itype};
            subplot(ndist_types,nrules,irule+(itype-1)*nrules);
            if (itrial_ptr>0)
                ntrials=trials_list(itrial_ptr);
                cp_list=results{itype,irule}.C{itrial_ptr}/ntrials;
                bincenters_cp=[0:ntrials]/ntrials;
            else
                cp_list=results{itype,irule}.cp_table(:,4);
                bincenters_cp=[0.5:nbins_cp-0.5]*(1/(nbins_cp-1));
            end
            hist([cp_list;1-cp_list],bincenters_cp);
            set(gca,'XLim',[0 1]);
            set(gca,'XTick',[0 0.5 1]);
            if (itype==1)
                title(results{itype,irule}.rule_string_nice)
            end
            if (irule==1)
                ylabel(results{itype,irule}.dist_type,'Interpreter','none');
            end
        end %itype
    end %irule
    axes('Position',[0.01,0.02,0.01,0.01]); %for text
    text(0,0,cat(2,name_string,'; ',lab_string),'Interpreter','none','FontSize',10);
    axis off;
end
%ord_char_simul_plot: plot results of simulate analysis of rank order for toy models
%
% run after load results from the output workspace of ord_char_simul_demo
%
% lots of kludges to recalculate column labels, a_sel, h_sel if col_mode=2 (a,h vary across column, but also indexing
%      interacts with page of plot (sym, umi, adt(
% 
% 01May24: allow for manual setting of ordinate scales with ylims_override
% 09Aug24: allow for setting of ticks with ytick_override
% 10Aug24: allow for boxplots of surrogates even if normalizing by trials
% 11Aug24: add a mode (col_mode==2) in which columns range over prior 
%
% See also: ORD_CHAR_SIMUL_DEMO, PSG_UMI_TRIPLIKE_PLOTA, PSG_LIKE_ANALTABLE.
%
nrats=3; %sym, umi, adt: ratios to plot
nratsp=nrats+1; % also plot umi raw
%
if ~exist('col_mode')
    col_mode=1;
end
col_mode=getinp('column mode: 1 for decision rule, 2 for prior','d',[1 2],col_mode);
%plot details
if ~exist('ymarg') ymarg=[-0.05 0.05]; end
if ~exist('xmarg') xmarg=[-0.5 0.5]; end
if ~exist('a_posit') a_posit=0.3; end
if ~exist('h_posit') h_posit=0.1; end
if ~exist('fontsize_ahlabel') fontsize_ahlabel=7;end
if ~exist('box_halfwidth') box_halfwidth=0.06; end %box half-width for "flip all" surrogates
if ~exist('ylims_override') ylims_override=NaN(nratsp,2); end %ylims_override(suar,:) are low and hi ranges for sym, ultra, addtree, ultra_raw
if ~exist('ytick_override') ytick_override=cell(nratsp,1); end
% 
%results{itype,irule}=r is how data are saved in ord_char_simul_demo
%
geometry_set_def='tree'; %in case missing from r
%
nrules=size(results,2);
ndist_types=size(results,1);

if ~isfield(results{1,1},'geometry_set')
    geometry_set=geometry_set_def;
else
    geometry_set=results{1,1}.geometry_set;
end
disp(sprintf('geometry set is %s',geometry_set));
for itype=1:ndist_types
    disp(sprintf(' distance type %1.0f-> %s',itype,results{itype,1}.dist_type));
end
type_ptrs=getinp('distance type(s) to plot','d',[1 ndist_types],[1:ndist_types]);
%
for irule=1:nrules
    disp(sprintf(' rule %1.0f-> %s',irule,results{1,irule}.rule_string));
end
switch col_mode
    case 1
        rule_ptrs=getinp('decision rule(s) to plot','d',[1 nrules],[1:nrules]);
    case 2
        rule_ptrs=getinp('one decision rule to plot','d',[1 nrules],ceil(nrules/2));
end
%
npts=length(unique(results{1,1}.cp_table(:,1)));
ntriplets=nchoosek(npts,3); %exhaustively sample the triplets
ntriads=3*ntriplets;
ntents=4*nchoosek(npts,4); %exhaustively sample tents
disp(sprintf(' npts %6.0f      ntriads %6.0f       ntriplets %6.0f      ntents %6.0f',npts,ntriads,ntriplets,ntents));
trials_list_all=results{1,1}.trials_list;
disp('trials list (all)');
disp(trials_list_all);
ntrials_list_all=length(trials_list_all);
trials_list_ptrs=getinp('pointers to above list','d',[1 ntrials_list_all],[1:ntrials_list_all]);
ntrials_list_sel=length(trials_list_ptrs);
trials_list_sel=trials_list_all(trials_list_ptrs);
disp('trials list (selected)')
disp(trials_list_sel);
%
xval_mode=getinp('1 for x-axis is log2(trials), 2 for each simulation separately','d',[1 2],...
    2-double(length(unique(trials_list_sel))==length(trials_list_sel)));
switch xval_mode
    case 1
        xvals=log2(trials_list_sel);
        xvals_label=trials_list_sel;
    case 2
        xvals=[1:length(trials_list_sel)];
        xvals_label=trials_list_sel;
end
xrange=[xvals(1)+xmarg(1),xvals(end)+xmarg(2)];
%
if ~exist('if_norm_trials') if_norm_trials=0; end
if_norm_trials=getinp('1 to normalize by number of trials, -1 to normalize by trials and show surrogates','d',[-1 1],if_norm_trials);
if if_norm_trials<=0
    if ~exist('if_show_flipall') if_show_flipall=[0 1 1 1]; end %show flip-all surrogates for sym, ultra, addtree, raw-ultra
    if ~exist('if_show_flipany') if_show_flipany=[1 1 1 1]; end %show flip-any surrogates for sym, ultra, addtree, raw-ultra
else %don't show surrogates, so axis is expanded
    if ~exist('if_show_flipall') if_show_flipall=[0 0 0 0]; end
    if ~exist('if_show_flipany') if_show_flipany=[0 0 0 0]; end
end
%
a_fixlist=results{1,1}.a_fixlist;
disp('a_fixlist');
disp(a_fixlist);
h_fixlist=results{1,1}.h_fixlist;
disp('h_fixlist');
disp(h_fixlist);
switch col_mode
    case 1
        ah_pairs=1;
            a_sel_list=getinp('choice for a (0: fit a)','d',[0 length(a_fixlist)],0);
            h_selaug_list=getinp('choice for h (0: fit h, -1: use h=0 for sym and adt, lowest h>0 for umi)','d',[-1 length(h_fixlist)],-1);
      case 2
        ah_pairs=getinp('number of (a,h) pairs for columns','d',[1 10]);
            a_sel_list=getinp('list of choices for a (0: fit a)','d',[0 length(a_fixlist)],zeros(1,ah_pairs));
            h_selaug_list=getinp('list of choices for h (0: fit h, -1: use h=0 for sym and adt, lowest h>0 for umi)','d',[-1 length(h_fixlist)],repmat(-1,1,ah_pairs));
end
%use initial values for a, h; will be recalculated later if col_mode=2
a_sel=a_sel_list(1);
h_selaug=h_selaug_list(1);
if (a_sel==0) a_string='fit'; else a_string=sprintf('%5.3f',a_fixlist(a_sel)); end
if (h_selaug==0)
    h_sel_sa=h_selaug;
    h_sel_umi=h_sel_sa;
    h_string_sa='fit';
    h_string_umi=h_string_sa;
elseif (h_selaug>0)
    h_sel_sa=h_selaug;
    h_sel_umi=h_sel_sa;
    h_string_sa=sprintf('%6.4f',h_fixlist(h_sel_sa));
    h_string_umi=h_string_sa;
else
    h_sel_sa=min(find(h_fixlist==0));
    h_sel_umi=min(find(h_fixlist==min(h_fixlist(h_fixlist>0))));
    h_string_sa=sprintf('%6.4f',h_fixlist(h_sel_sa));
    h_string_umi=sprintf('%6.4f',h_fixlist(h_sel_umi));
end
if (col_mode==1) ncols=length(rule_ptrs); end
if (col_mode==2) ncols=ah_pairs; end
%
%fill in for backward compatibility
for itype=1:ndist_types
    for irule=1:nrules
        if ~isfield(results{itype,irule},'geometry_set')
            results{itype,irule}.geometry_set=geometry_set_def;
        end
        if ~isfield(results{itype,irule},'trials_if_poisson')
            results{itype,irule}.trials_if_poisson=0;
        end
        if ~isfield(results{itype,irule},'triplet_counts')
            results{itype,irule}.triplet_counts=repmat(ntriplets,1,ntrials_list_all);
        end
        if ~isfield(results{itype,irule},'tent_counts')
            results{itype,irule}.tent_counts=repmat(ntents,1,ntrials_list_all);
        end
    end
end
if_poisson=results{1,1}.trials_if_poisson;
%log2(trials_list_sel
% for each index (Isym, Iumi, Iadt, raw Iumi), show how llr depends on
% number of trials per triad, for each geometry type and each decision rule
nsurrs=length(results{1,1}.llr_d2);
%
for suar=1:nratsp %sym, umi, adt, umi raw
    if_sublogh=0;
    switch suar
        case 1
            sua_brief='sym';
            sua_string='symmetry';
            nsets_name='triplet_counts';
            apriori_sub=log(3/4);
            apriori_plot=apriori_sub;
            h_string=h_string_sa;
            h_sel=h_sel_sa;
        case 2
            if_sublogh=1;
            sua_brief='umi';
            sua_string='ultrametric';
            nsets_name='triplet_counts';
            apriori_sub=0;
            apriori_plot=apriori_sub;
            h_string=h_string_umi;
            h_sel=h_sel_umi;
        case 3
            sua_brief='adt';
            sua_string='addtree';
            nsets_name='tent_counts';
            apriori_sub=log(2/3);
            apriori_plot=apriori_sub;
            h_string=h_string_sa;
            h_sel=h_sel_sa;
        case 4
            sua_brief='umi';
            sua_string='ultrametric (raw)';
            nsets_name='triplet_counts';
            apriori_sub=0;
            apriori_plot=NaN;
            h_string=h_string_umi;
            h_sel=h_sel_umi;
    end
    switch col_mode
        case 1
            name_string=sprintf('%s %s: a=%s h=%s; normalize by trials: %1.0f; poisson counts: %1.0f',...
                geometry_set,sua_string,a_string,h_string,if_norm_trials,if_poisson);
        case 2
            name_string=sprintf('%s %s: decision rule: %s; normalize by trials: %1.0f; poisson counts: %1.0f',...
                geometry_set,sua_string,results{1,rule_ptrs(1)}.rule_string,if_norm_trials,if_poisson);
    end
    figure;
    set(gcf,'Position',[100 100 1200 800]);
    set(gcf,'Numbertitle','off');
    set(gcf,'Name',name_string);
    yrange=[Inf,-Inf];
    ha=cell(nrules,ndist_types);
    for icol=1:ncols
        switch col_mode
            case 1
                irule_ptr=icol;
                a_sel=a_sel_list(1);
                h_selaug_list=h_selaug_list(1);
            case 2
                irule_ptr=1;
                a_sel=a_sel_list(1+mod(icol-1,ah_pairs));
                h_selaug=h_selaug_list(1+mod(icol-1,ah_pairs));
                if (a_sel==0) a_string='fit'; else a_string=sprintf('%5.3f',a_fixlist(a_sel)); end
                if (h_selaug==0)
                    h_sel_sa=h_selaug;
                    h_sel_umi=h_sel_sa;
                    h_string_sa='fit';
                    h_string_umi=h_string_sa;
                elseif (h_selaug>0)
                    h_sel_sa=h_selaug;
                    h_sel_umi=h_sel_sa;
                    h_string_sa=sprintf('%6.4f',h_fixlist(h_sel_sa));
                    h_string_umi=h_string_sa;
                else
                    h_sel_sa=min(find(h_fixlist==0));
                    h_sel_umi=min(find(h_fixlist==min(h_fixlist(h_fixlist>0))));
                    h_string_sa=sprintf('%6.4f',h_fixlist(h_sel_sa));
                    h_string_umi=sprintf('%6.4f',h_fixlist(h_sel_umi));
                end
                switch suar
                    case 1
                        h_string=h_string_sa;
                        h_sel=h_sel_sa;
                    case 2
                        h_string=h_string_umi;
                        h_sel=h_sel_umi;
                    case 3
                        h_string=h_string_sa;
                        h_sel=h_sel_sa;
                    case 4
                        h_string=h_string_umi;
                        h_sel=h_sel_umi;
                end
        end
        irule=rule_ptrs(irule_ptr);
        %
        for itype_ptr=1:length(type_ptrs)
            itype=type_ptrs(itype_ptr);
            nsets=results{itype,irule}.(nsets_name); %allow for different number of sets for each trial count
            %
            ha{icol,itype}=subplot(length(type_ptrs),ncols,icol+(itype_ptr-1)*ncols);
            %
            rd=results{itype,irule}.(cat(2,'llr_',sua_brief)); %log likelihood for sym, umi, or adt
            rd_vm=results{itype,irule}.(cat(2,'llr_vm_',sua_brief)); %variance of the within-[triplet|tent] mean
            rd_sv=results{itype,irule}.(cat(2,'llr_sv_',sua_brief)); %sum of the within-[triplet|tent] variances
            rd_vv=results{itype,irule}.(cat(2,'llr_vv_',sua_brief)); %variance of the within-[triplet|tent] variances
            %
            rdp=rd(:,:,a_sel+1,h_sel+1)./repmat(nsets(:),1,nsurrs); %d2 is real data or surrogates           
            if (if_sublogh)
                sublogh=log(results{itype,irule}.dirichlet.h(a_sel+1,h_sel+1,:));
                sublogh=sublogh(:);
            else
                sublogh=zeros(ntrials_list_all,1);
            end
            rdp=rdp-repmat(sublogh,1,size(rdp,2)); %subtract offset if Iumi is plotted (suar=2)
            %error bars for surrogates: from psg_umi_triplike_plota:
            %   vars=ruse{2,ithr_type}(:,:,ihfix); %d1: threshold level, d2: surrogate type
            %   eb_stds=sqrt(vars)./repmat(nsets,1,nsurr+nconform);
            ebs_stds=sqrt(rd_sv(:,:,a_sel+1,h_sel+1))./repmat(nsets(:),1,nsurrs);
            %
            if (if_norm_trials~=0)
                rdp=apriori_sub+(rdp-apriori_sub)./repmat(trials_list_all(:),1,nsurrs); %normalize by number of trials per triad
                ebs_stds=ebs_stds./repmat(trials_list_all(:),1,nsurrs);
            end
            %
            %from psg_like_analtable:
            %hs=psg_like_plot(abscissa_val+opts.box_halfwidth*[-1 1 1 -1 -1],data.flip_all+data.flip_all_sd*[1 1 -1 -1 1],'k',data.h,d23);
            hp=plot(xvals,rdp(trials_list_ptrs,1),'ks');
            set(hp,'MarkerFaceColor','k');
            hold on;
            for is=2:3 %is=2: flip all, is=3: flip all
                bwm=is-1; %double the box width if (if_show_flipall(suar)) %show flip-all surrogates?
                if (is==2) if_show_flip=if_show_flipall(suar); end
                if (is==3) if_show_flip=if_show_flipany(suar); end
                if (if_show_flip)
                    for itrial=1:ntrials_list_sel
                        plot(xvals(itrial)+bwm*box_halfwidth*[-1 1 1 -1 -1],rdp(trials_list_ptrs(itrial),is)+ebs_stds(trials_list_ptrs(itrial),is)*[1 1 -1 -1 1],'k'); %flip all
                        plot(repmat(xvals(itrial),1,2),rdp(trials_list_ptrs(itrial),[1 is]),'k');
                    end
                end
            end
            %
            plot(xrange,repmat(apriori_plot,1,2),'k:'); % a priori line
            %
            yr=get(gca,'YLim');
            yrange=[min(yrange(1),yr(1)),max(yrange(2),yr(2))];
            %
            set(gca,'XLim',xrange);
            set(gca,'XTick',xvals);
            set(gca,'XTickLabel',xvals_label);
            if (itype_ptr==1)
                if (col_mode==1) col_header=results{itype,irule}.rule_string_nice; end
                if (col_mode==2) col_header=cat(2,'a=',a_string,' h=',h_string); end
                title(col_header);
            end
            if (itype_ptr==length(type_ptrs))
                xlabel('trials per triad');
            end
            if (icol==1)
                ylabel(cat(2,'sim: ',results{itype,irule}.dist_type),'Interpreter','none');
            end
        end %itype_ptr
    end %irule_ptr
    ylims=yrange+ymarg;
    %uniformize axis scales and add labels for values of a and h if fitted
    if ~any(isnan(ylims_override(suar,:))) %override on ylims for plotting
        ylims=ylims_override(suar,:);
    end
    for icol=1:ncols
        switch col_mode
            case 1
                irule_ptr=icol;
            case 2 %need to recalculate a,h selections
                irule_ptr=1;
                a_sel=a_sel_list(1+mod(icol-1,ah_pairs));
                h_selaug=h_selaug_list(1+mod(icol-1,ah_pairs));
                if (h_selaug==0)
                    h_sel_sa=h_selaug;
                    h_sel_umi=h_sel_sa;
                    h_string_sa='fit';
                    h_string_umi=h_string_sa;
                elseif (h_selaug>0)
                    h_sel_sa=h_selaug;
                    h_sel_umi=h_sel_sa;
                    h_string_sa=sprintf('%6.4f',h_fixlist(h_sel_sa));
                    h_string_umi=h_string_sa;
                else
                    h_sel_sa=min(find(h_fixlist==0));
                    h_sel_umi=min(find(h_fixlist==min(h_fixlist(h_fixlist>0))));
                    h_string_sa=sprintf('%6.4f',h_fixlist(h_sel_sa));
                    h_string_umi=sprintf('%6.4f',h_fixlist(h_sel_umi));
                end
                switch suar
                    case 1
                        h_string=h_string_sa;
                        h_sel=h_sel_sa;
                    case 2
                        h_string=h_string_umi;
                        h_sel=h_sel_umi;
                    case 3
                        h_string=h_string_sa;
                        h_sel=h_sel_sa;
                    case 4
                        h_string=h_string_umi;
                        h_sel=h_sel_umi;
                end

        end
        irule=rule_ptrs(irule_ptr);
        for itype_ptr=1:length(type_ptrs)
            itype=type_ptrs(itype_ptr);
            axes(ha{icol,itype});
            set(gca,'YLim',ylims);
            if ~isempty(ytick_override{suar})
                set(gca,'YTick',ytick_override{suar});
            end
            %add values of a and h if fit, one below each point
            if a_sel==0
                a_vals=results{itype,irule}.dirichlet.a(a_sel+1,h_sel+1,:);
                for itrial=1:ntrials_list_sel
                    if ~isnan(a_vals(trials_list_ptrs(itrial)))
                        text(xvals(itrial),ylims(1)+a_posit*diff(ylims),sprintf('%4.2f',a_vals(trials_list_ptrs(itrial))),'FontSize',fontsize_ahlabel);
                    end
                end
                text(xvals(1)+xmarg(1)*0.75,ylims(1)+a_posit*diff(ylims),'a=','FontSize',fontsize_ahlabel);
            end
            if h_sel==0
                h_vals=results{itype,irule}.dirichlet.h(a_sel+1,h_sel+1,:);
                for itrial=1:ntrials_list_sel
                    if ~isnan(h_vals(trials_list_ptrs(itrial)))
                        text(xvals(itrial),ylims(1)+h_posit*diff(ylims),sprintf('%4.2f',h_vals(trials_list_ptrs(itrial))),'FontSize',fontsize_ahlabel);
                    end
                end
                text(xvals(1)+xmarg(1)*0.75,ylims(1)+h_posit*diff(ylims),'h=','FontSize',fontsize_ahlabel);
            end
        end
    end
    %
    axes('Position',[0.01,0.02,0.01,0.01]); %for text
    text(0,0,name_string,'Interpreter','none','FontSize',10);
    axis off;
end

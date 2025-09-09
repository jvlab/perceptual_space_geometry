%psg_isomap_demo: embed psg coordinates via isomap
%
% Demonstrate psg_isomap, which requires graph toolikit
%
% looks at eigenvalues of MDS solution of standard and isomap embeddings
% plots:
% * eigenvalues (prior to square root) for mds of standard and isomap distances
% * fraction of non-hyperbolic power: sum of positive eigs/sum of abs(eigs) for mds of standard and isomap distances
% * participation ratio: sum(abs(eigs))^2/(sum of eigs^2) for mds of standardd and isomap distances
% * edge ratio: number of edges per node/nn for isomap graph
% as a function of minimum number of nearest neighbors (nn) in the isomap graph
%
%  Should be applicable to btc datasets and also the 5-domain Waraich datasets
%
% Requires graph toolkit.
%
% 05Sep25: Modularize psg_isomap
%
% See also: PSG_GET_COORDSETS, COOTODSQ, PSG_ISOMAP.
%
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_qpred') opts_qpred=struct(); end %for psg_qformpred
if ~exist('opts_isomap') opts_isomap=struct(); end %for psg_isomap
%
opts_qpred=filldefault(opts_qpred,'qform_datafile_def','../stim/btc_allraysfixedb_avg_100surrs_madj.mat');
opts_qpred=filldefault(opts_qpred,'qform_modeltype',12); %if_symm=1 (symmetrize around origin), if_axes=1 (symmetrize bc, de, tuvw); ifaug=1 (augmented coords)
%
if ~exist('neigs_show_max') neigs_show_max=10; end
%
opts_read=filldefault(opts_read,'if_spray',0); %default to not define rays by single points
[sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred);
nsets=length(sets);
%
% choose datasets and dimensions
%
disp('datasets available:');
for iset=1:nsets
    disp(sprintf(' set %2.0f: dim range [%3.0f %3.0f] label: %s',iset,min(sets{iset}.dim_list),max(sets{iset}.dim_list),sets{iset}.label));
end
% for each dataset, get dimension list to analyze, min/max number of nearest neighbors
nstims=zeros(1,nsets);
dim_list=sets{1}.dim_list;
for iset=1:nsets
    nstims(iset)=sets{1}.nstims;
    dim_list=intersect(dim_list,sets{iset}.dim_list);
end
dim_list_anal=getinp('dims to analyze','d',[1 max(dim_list)],[2:max(dim_list)]);
nn_list=getinp('nth nearest neighbor list','d',[1 min(nstims)-1],[1:4]);
nbr_mindist=zeros(max(dim_list),max(nn_list),nsets);
nbr_mtx_nonzero=zeros(max(dim_list),max(nn_list),nsets);
eivals_iso=cell(max(dim_list),max(nn_list),nsets);
coords_iso=cell(max(dim_list),max(nn_list),nsets);
eivals_std=cell(max(dim_list),1,nsets);
coords_std=cell(max(dim_list),1,nsets);
%
% compute distances based on networks with connectivity via restricted number of nearest neighbors
%
for iset=1:nsets
    for dim_ptr=1:length(dim_list_anal);
        dim_anal=dim_list_anal(dim_ptr);
        coords=ds{iset}{dim_anal};
        dists=sqrt(cootodsq(coords));
        for nn_ptr=1:length(nn_list)
            nn=nn_list(nn_ptr);
            opts_isomap.min_nn=nn;
            [eivals_iso{dim_anal,nn,iset},coords_iso{dim_anal,nn,iset},opts_isomap_used]=psg_isomap(dists,opts_isomap);
            nbr_mindist(dim_anal,nn,iset)=opts_isomap_used.nbr_mindist;
            nbr_mtx=opts_isomap_used.nbr_mtx;
            nbr_mtx_nonzero(dim_anal,nn,iset)=(sum(nbr_mtx(:)>0)/2);
        end %nn_ptr
        opts_isomap.min_nn=0;
        [eivals_std{dim_anal,1,iset},coords_std{dim_anal,1,iset}]=psg_isomap(dists,opts_isomap);
    end %dim_ptr
end %iset
%
%plot
%
labels_short=cell(1,nsets+2);
for iset=1:nsets
    last_sep=max([0,max(find(sets{iset}.label=='/')),max(find(sets{iset}.label=='\'))]);
    labels_short{iset}=sets{iset}.label(last_sep+1:end);
end
labels_short{iset+1}='mean';
labels_short{iset+2}='median';
nrows=1+length(nn_list);
ncols=2;
neigs_show=min(min(nstims),neigs_show_max);
rng('default');
colors=rand(3,nsets)'*0.75;
pwr_ratio=zeros(max(dim_list),1+max(nn_list),nsets); %fraction of power that is in the positive eigs
part_ratio=zeros(max(dim_list),1+max(nn_list),nsets); %participation ratio
edge_ratio=zeros(max(dim_list),1+max(nn_list),nsets); %fraction of edges that are occupied in the isomap graph
for dim_ptr=1:length(dim_list_anal)
    dim_anal=dim_list_anal(dim_ptr);
    fig_label=sprintf('isomap analysis, dim=%2.0f',dim_anal);
    figure;
    set(gcf,'Position',[50 50 1000 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',fig_label);
    %plot normalized eigenvalues, for each nn value
    for nn_ptr=0:length(nn_list)
        subplot(nrows,ncols,1+(nn_ptr*ncols));
        if nn_ptr==0
            eigs_plot=eivals_std(dim_anal,1,:);
            lab='std';
        else
            nn=nn_list(nn_ptr);
            eigs_plot=eivals_iso(dim_anal,nn,:);
            lab=sprintf('iso, nn=%2.0f',nn);
        end
        neigs_show_all=0;
        for iset=1:nsets
            eigs_show=eigs_plot{1,1,iset};
            neigs_show=min(length(eigs_show),neigs_show_max);
            hp=plot([1:neigs_show],eigs_show(1:neigs_show)/max(abs(eigs_show)),'k','LineWidth',1); %normalize eigenvalues by max
            set(hp,'Color',colors(iset,:));
            hold on;
            neigs_show_all=max(neigs_show_all,neigs_show);
        end
        plot([1 neigs_show_all],[0 0],'k:');
        xlabel('eiv');
        set(gca,'XLim',[1 neigs_show_all]);        
        set(gca,'YLim',[-0.25 1]);
        set(gca,'XTick',[1:neigs_show_all]);
        title(lab);
    end
    %compute and plot fraction of power in pos eigs, as function of nn
    subplot(2,2,2);
    pwr_ratio(dim_anal,1,:)=1;
    for iset=1:nsets
        for nn_ptr=1:length(nn_list)
            nn=nn_list(nn_ptr);
            eivs=eivals_iso{dim_anal,nn,iset};
            pwr_ratio(dim_anal,1+nn,iset)=sum(eivs(eivs>0))/sum(abs(eivs));
        end
        hp=plot([0 nn_list],pwr_ratio(dim_anal,1+[0 nn_list],iset),'k','LineWidth',1);
        set(hp,'Color',colors(iset,:));
        hold on;
    end
    plot([0 nn_list],mean(pwr_ratio(dim_anal,1+[0 nn_list],:),3),'k','LineWidth',2);
    plot([0 nn_list],median(pwr_ratio(dim_anal,1+[0 nn_list],:),3),'k:','LineWidth',2);
    legend(labels_short,'FontSize',6,'Interpreter','none','Location','SouthEast');
    xlabel('number of nrst nbrs');
    ylabel('Euclidean power frac');
    set(gca,'XLim',[0 max(nn_list)]);
    set(gca,'XTick',[0 nn_list]);
    set(gca,'XTickLabel',strvcat('std',num2str(nn_list(:))));
    set(gca,'YLim',[0.5 1]); %at least half of the power must be in the positive eigs
    title('frac of power in pos eivs');
    %compute and plot participation ratio as functon of nn
    subplot(4,2,6);
    for iset=1:nsets
        for nn_ptr=1:length(nn_list)
            nn=nn_list(nn_ptr);
            eivs=eivals_iso{dim_anal,nn,iset};
            part_ratio(dim_anal,1+nn,iset)=(sum(abs(eivs)).^2)/sum(eivs.^2);
        end
        eivs=eivals_std{dim_anal,1,iset};
        part_ratio(dim_anal,1,iset)=(sum(abs(eivs)).^2)/sum(eivs.^2);
        hp=plot([0 nn_list],part_ratio(dim_anal,1+[0 nn_list],iset),'k','LineWidth',1);
        set(hp,'Color',colors(iset,:));
        hold on;
    end
    plot([0 nn_list],mean(part_ratio(dim_anal,1+[0 nn_list],:),3),'k','LineWidth',2);
    plot([0 nn_list],median(part_ratio(dim_anal,1+[0 nn_list],:),3),'k:','LineWidth',2);
    xlabel('number of nrst nbrs');
    ylabel('dimension');
    set(gca,'XLim',[0 max(nn_list)]);
    set(gca,'XTick',[0 nn_list]);
    set(gca,'XTickLabel',strvcat('std',num2str(nn_list(:))));
    set(gca,'YLim',[0 max(get(gca,'YLim'))]);
    title('participation ratio');
    %compute and plot edge fraction as function of nn
    subplot(4,2,8);
    for iset=1:nsets
        for nn_ptr=1:length(nn_list)
            nn=nn_list(nn_ptr);
            edge_ratio(dim_anal,1+nn,iset)=2*nbr_mtx_nonzero(dim_anal,nn,iset)/nstims(iset)/nn;
        end
        edge_ratio(dim_anal,1,iset)=1;
        hp=plot([0 nn_list],log2(edge_ratio(dim_anal,1+[0 nn_list],iset)),'k','LineWidth',1);
        set(hp,'Color',colors(iset,:));
        hold on;
    end
    plot([0 nn_list],log2(mean(edge_ratio(dim_anal,1+[0 nn_list],:),3)),'k','LineWidth',2);
    plot([0 nn_list],log2(median(edge_ratio(dim_anal,1+[0 nn_list],:),3)),'k:','LineWidth',2);
    xlabel('number of nrst nbrs');
    ylabel('log2(edge density)');
    set(gca,'XLim',[0 max(nn_list)]);
    set(gca,'XTick',[0 nn_list]);
    set(gca,'XTickLabel',strvcat('std',num2str(nn_list(:))));
    set(gca,'YLim',[0 ceil(log2(max(nstims)))]);
    title('edge ratio');
    %
    axes('Position',[0.01,0.04,0.01,0.01]); %for text
    text(0,0,fig_label,'Interpreter','none','FontSize',8);
    axis off;
end

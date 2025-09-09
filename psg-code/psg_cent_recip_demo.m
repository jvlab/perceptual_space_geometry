%psg_cent_recip_demo: demonwstrate omputation of centrality and reciprocity indices
% 
% centrality and reciprocity indices, from
% Tversky, A. & Hutchinson, J. W. Nearest neighbor analysis of
% psychological spaces (1988) Psychological Review, 93(1), 3–22. 
%
%  Should be applicable to btc datasets and also the 5-domain Waraich datasets
%
% See also: PSG_GET_COORDSETS, COOTODSQ, PSG_CENT_RECIP.
%    
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_qpred') opts_qpred=struct(); end %for psg_qformpred
if ~exist('opts_cr') opts_cr=struct(); end % for psg_cent_recip
%
opts_qpred=filldefault(opts_qpred,'qform_datafile_def','../stim/btc_allraysfixedb_avg_100surrs_madj.mat');
opts_qpred=filldefault(opts_qpred,'qform_modeltype',12); %if_symm=1 (symmetrize around origin), if_axes=1 (symmetrize bc, de, tuvw); ifaug=1 (augmented coords)
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
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
if (if_frozen~=0)
    rng('default');
    if (if_frozen<0)
        rand(1,abs(if_frozen));
    end
else
    rng('shuffle');
end
%
% compute centrality and reciprocity indices for each dimension
%
ind_cent=zeros(max(dim_list),nsets);
ind_recip=zeros(max(dim_list),nsets);
ind_cent_unitized=zeros(max(dim_list),nsets); %ind_cent transformed to [0 1]
ind_recip_unitized=zeros(max(dim_list),nsets); %ind_recip transformed to [0 1]
opts_rc_used=cell(max(dim_list),nsets);
for iset=1:nsets
    for dim_ptr=1:length(dim_list_anal);
        dim_anal=dim_list_anal(dim_ptr);
        coords=ds{iset}{dim_anal};
        dists=sqrt(cootodsq(coords));
        [ind_cent(dim_anal,iset),ind_recip(dim_ana,iset),opts_rc_used{dim_anal,iset}]=psg_cent_recip(dists,opts_rc);
        ind_cent_unitized(dim_anal,iset)=(ind_cent(dim_anal,iset)-1)./(opts.rc_used.ind_cent_max-1); %0: no central elements, 1: a single central element (star)
        ind_recip_unitized(dim_anal,iset)=1-(ind_recip(dim_anal,iset)-1)./(opts.rc_used.ind_recip_max-1); %0: non-reciprrocal, 1: reciprocal
    end %dim_ptr
end %iset
% disp(sprintf('total ties: %3.0f',nties_tot));
% %
% %plot
% %
% labels_short=cell(1,nsets+2);
% for iset=1:nsets
%     last_sep=max([0,max(find(sets{iset}.label=='/')),max(find(sets{iset}.label=='\'))]);
%     labels_short{iset}=sets{iset}.label(last_sep+1:end);
% end
% labels_short{iset+1}='mean';
% labels_short{iset+2}='median';
% ncols=2;
% rng('default');
% colors=rand(3,nsets)'*0.75;
% %
% fig_label='centrality and reciprocity';
% figure;
% set(gcf,'Position',[50 50 1000 800]);
% set(gcf,'NumberTitle','off');
% set(gcf,'Name',fig_label);
% for dim_ptr=1:length(dim_list_anal)
%     dim_anal=dim_list_anal(dim_ptr);
%     % %plot normalized eigenvalues, for each nn value
%     % for nn_ptr=0:length(nn_list)
%     %     subplot(nrows,ncols,1+(nn_ptr*ncols));
%     %     if nn_ptr==0
%     %         eigs_plot=eivals_std(dim_anal,1,:);
%     %         lab='std';
%     %     else
%     %         nn=nn_list(nn_ptr);
%     %         eigs_plot=eivals_iso(dim_anal,nn,:);
%     %         lab=sprintf('iso, nn=%2.0f',nn);
%     %     end
%     %     neigs_show_all=0;
%     %     for iset=1:nsets
%     %         eigs_show=eigs_plot{1,1,iset};
%     %         neigs_show=min(length(eigs_show),neigs_show_max);
%     %         hp=plot([1:neigs_show],eigs_show(1:neigs_show)/max(abs(eigs_show)),'k','LineWidth',1); %normalize eigenvalues by max
%     %         set(hp,'Color',colors(iset,:));
%     %         hold on;
%     %         neigs_show_all=max(neigs_show_all,neigs_show);
%     %     end
%     %     plot([1 neigs_show_all],[0 0],'k:');
%     %     xlabel('eiv');
%     %     set(gca,'XLim',[1 neigs_show_all]);        
%     %     set(gca,'YLim',[-0.25 1]);
%     %     set(gca,'XTick',[1:neigs_show_all]);
%     %     title(lab);
%     % end
%     %compute and plot fraction of power in pos eigs, as function of nn
%     %
%     axes('Position',[0.01,0.04,0.01,0.01]); %for text
%     text(0,0,fig_label,'Interpreter','none','FontSize',8);
%     axis off;
end

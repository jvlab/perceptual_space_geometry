%psg_procrustes_regr_demo: apply Procrustes and regression routines
%
% results contains key outputs
%
%   See also:  PROCRUSTES, REGRESS, PERSP_XFORM_FIND, PSG_GET_COORDSETS, PSG_PROCRUSTES_REGR_TEST,PSG_PARSE_FILENAME.
%
if ~exist('min_dim') min_dim=1; end
if ~exist('max_dim') max_dim=3; end
if ~exist('if_regr_const') if_regr_const=1; end
if ~exist('if_center') if_center=1; end
if ~exist('nshuff_procrustes') nshuff_procrustes=1000; end
if ~exist('nshuff_resids') nshuff_resids=1000; end
if ~exist('nshuff_quantile') nshuff_quantile=0.05; end
if ~exist('opts_read') opts_read=struct; end
if ~exist('opts_rays') opts_rays=struct; end
if ~exist('opts_qpred') opts_qpred=struct; end
if ~exist('opts_persp') opts_persp=struct; end
opts_persp=filldefault(opts_persp,'method','fmin');
opts_persp=filldefault(opts_persp,'if_cycle',1); %in case method='oneshot'
%
if ~exist('if_pairplots') if_pairplots=1; end
if ~exist('color_data') color_data=[0 0 0]; end
if ~exist('color_procrustes') color_procrustes=[0 0.5 0]; end
if ~exist('color_check') color_check=[1 0 0]; end
if ~exist('color_affine') color_affine=[1 0 0.7]; end
if ~exist('color_projective') color_projective=[0 0 0.7]; end
if ~exist('line_width') line_width=1;end
%
fit_types={'affine','projective'};
%
if_builtin=getinp('use built-in datasets','d',[0 1]);
if ~exist('file_names') file_names={'./psg_data/bgca3pt_coords_MC_sess01_10.mat','./psg_data/bgca3pt_coords_BL_sess01_10.mat'}; end
if if_builtin
    nsets=2;
else
    disp('enter at least two datasets');
    [sets,ds,sas,rayss,opts_read_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred);
    nsets=length(sets);
end
%
min_dim=getinp('minimum dimension to use','d',[1 7],min_dim);
max_dim=getinp('maximum dimension to use','d',[min_dim 7],max(min_dim,max_dim));
results=struct;
subj_ids=cell(1,nsets);
paradigm_name=[];
if (if_builtin)
    results.names=file_names;
    for iset=1:nsets
        subj_ids{iset}=getfield(psg_parse_filename(file_names{iset}),'subj_id');
        if isempty(paradigm_name)
            paradigm_name=getfield(psg_parse_filename(file_names{iset}),'paradigm_name');
        end
    end
else
    for iset=1:nsets
        results.names{1,iset}=sets{iset}.label_long;
        if strcmp(sets{iset}.type,'data')
            subj_ids{iset}=getfield(psg_parse_filename(sets{iset}.label_long),'subj_id');
            if isempty(paradigm_name)
                paradigm_name=getfield(psg_parse_filename(sets{iset}.label_long),'paradigm_name');
            end
        else
            subj_ids{iset}='qform';
        end
    end
end
results.d_dim_names={'ref_ptr','adj_ptr','idim','shuffles'};
results.procrustes=struct();
results.affine=struct();
results.projective=struct();
%
for iset=1:nsets-1
    for jset=iset+1:nsets
        for ijrev=1:2
            switch ijrev
                case 1
                    ref_ptr=iset;
                    adj_ptr=jset;
                case 2
                    ref_ptr=jset;
                    adj_ptr=iset;
            end %ijrev
            if (if_builtin)
                ref_file=file_names{ref_ptr};
                adj_file=file_names{adj_ptr};
            else
                ref_file=sets{ref_ptr}.label_long;
                adj_file=sets{adj_ptr}.label_long;
            end
            disp('   ');
            disp(sprintf('reference dataset: %s',ref_file));
            disp(sprintf('dataset to be adjusted: %s',adj_file));
            for idim=min_dim:max_dim
                if if_builtin
                    ref=getfield(load(ref_file),sprintf('dim%1.0f',idim));
                    adj=getfield(load(adj_file),sprintf('dim%1.0f',idim));
                else
                    ref=ds{ref_ptr}{idim};
                    adj=ds{adj_ptr}{idim};
                end
                %
                tstring=sprintf(' ref %s dim %1.0f adj %s dim %1.0f regr const %1.0f center %1.0f',...
                    ref_file,idim,adj_file,idim,if_regr_const,if_center);
                disp(' ');
                disp(tstring);
                %
                rank_ref=rank(ref);
                rank_adj=rank(adj);
                if (idim>rank_ref) | (idim>rank_adj)
                    disp(sprintf(' dim %3.0f skipped.  rank(ref)=%3.0f rank(adj)=%3.0f',idim,rank_ref,rank_adj));
                    results.procrustes.d(ref_ptr,adj_ptr,idim)=0;
                    results.procrustes.d_shuff(ref_ptr,adj_ptr,idim,:)=zeros([1 1 1 nshuff_procrustes]);
                    for iafpe=1:2
                        results.(fit_types{iafpe}).d(ref_ptr,adj_ptr,idim)=0;
                        results.(fit_types{iafpe}).d_shuff(ref_ptr,adj_ptr,idim,:)=zeros([1,1,1,nshuff_resids]);
                    end
                else
                    %
                    npts=size(ref,1);
                    if if_center
                        ref=ref-repmat(mean(ref,1),npts,1);
                        adj=adj-repmat(mean(adj,1),npts,1);
                    end
                    %
                    [d,adj_proc,transform]=procrustes(ref,adj);
                    %disp('transform');
                    %disp(transform);
                    adj_check=transform.b*adj*transform.T+transform.c;
                    %
                    ref_pr=ref-adj_proc; %residuals after Procrustes
                    %
                    %recalculate normalized measure of deviation (d) from Procrustes
                    %
                    d_den=sum(sum((ref-repmat(mean(ref,1),npts,1)).^2,1));
                    d_num=sum(sum((ref_pr).^2,1));
                    d_check=d_num/d_den;
                    disp(sprintf('d (procrustes): %12.7f   d_check: %12.7f   diff: %12.7f',d,d_check,d-d_check))
                    results.procrustes.d(ref_ptr,adj_ptr,idim)=d;
                    %
                    % compare d with shuffle surrogates
                    %
                    if (nshuff_procrustes>0) & iset==1 & jset==2 & ijrev==1
                        perms_proc=zeros(nshuff_procrustes,npts); %set up a consistent set of shuffles across dimensions and dataset pairs
                        for ishuff=1:nshuff_procrustes
                            perms_proc(ishuff,:)=randperm(npts);
                        end
                    end
                    if nshuff_procrustes>0
                        d_shuff_proc=zeros(1,nshuff_procrustes);
                        for ishuff=1:nshuff_procrustes
                            d_shuff_proc(ishuff)=procrustes(ref,adj(perms_proc(ishuff,:),:));
                        end
                        disp(sprintf('d is >= d_shuff for %5.0f of %5.0f shuffles, d=%12.7f, min(shuffle)=%12.7f',...
                            sum(double(d>=d_shuff_proc)),nshuff_procrustes,d,min(d_shuff_proc)));
                        results.procrustes.d_shuff(ref_ptr,adj_ptr,idim,:)=reshape(d_shuff_proc,[1 1 1 nshuff_procrustes]);
                    end
                    %
                    % fit an affine or perspective (projective) transformation to original data and to residuals
                    %
                    adj_xform=cell(1,2);
                    d_afpe=zeros(2,2); %iorre,iafpe
                    for iafpe=1:2  % 1: affine, 2: perspective
                        for iorre=1:2 % 1: original, 2: residuals
                            switch iorre
                                case 1
                                    ref_tofit=ref;
                                    ref_type='from scratch';
                                case 2
                                    ref_tofit=ref_pr;
                                    ref_type='Procrustes residual';
                            end
                            switch iafpe
                                case 1
                                    c=zeros(size(adj,2),1); %affine transform is special case of projective transform 
                                    [a,b]=persp_fit(c,adj,ref_tofit);
                                case 2
                                    [persp,y_fit]=persp_xform_find(adj,ref_tofit,opts_persp);
                                    a=persp(1:end-1,1:end-1);
                                    b=persp(end,1:end-1);
                                    c=persp(1:end-1,end);
                            end
                            adj_fit=persp_apply(a,b,c,adj);
                            if (iorre==1)
                                adj_xform{iafpe}=adj_fit;
                            end
                            d_afpe(iorre,iafpe)=sum((ref_tofit(:)-adj_fit(:)).^2)/d_den;
                            %
                            if (iafpe==2)
                                persp_yfit_maxdev=max(abs(adj_fit(:)-y_fit(:)));
                                disp(sprintf('fitting %20s with %12s transformation: d=%12.7f, max dev for persp_yfit: %12.7f',...
                                    ref_type,fit_types{iafpe},d_afpe(iorre,iafpe),persp_yfit_maxdev));
                            else
                                disp(sprintf('fitting %20s with %12s transformation: d=%12.7f',ref_type,fit_types{iafpe},d_afpe(iorre,iafpe)));
                            end
                        end %orig or resids
                    end %affine or perspective
                    for iafpe=1:2
                        results.(fit_types{iafpe}).d(ref_ptr,adj_ptr,idim)=d_afpe(1,iafpe);
                    end
                    %
                    % fit an affine or perspective (projective) transformation to original data with shuffled residuals
                    %
                    if (nshuff_resids>0) & iset==1 & jset==2 & ijrev==1
                        perms_resids=zeros(nshuff_resids,npts); %set up a consistent set of shuffles across dimensions and dataset pairs
                        for ishuff=1:nshuff_resids
                            perms_resids(ishuff,:)=randperm(npts);
                        end
                    end
                    if nshuff_resids>0
                        for ishuff=1:nshuff_resids
                            if (ishuff==3)
                                perms_resids(ishuff,:)=[1:npts];
                            end
                            ref_tofit=adj_proc+ref_pr(perms_resids(ishuff,:),:); %reconstitute with shuffled residuals; ref_pr=ref-adj_proc; %residuals after Procrustes
                            for iafpe=1:2  % 1: affine, 2: perspective
                                switch iafpe
                                    case 1
                                        c=zeros(size(adj,2),1); %affine transform is special case of projective transform 
                                        [a,b]=persp_fit(c,adj,ref_tofit);
                                    case 2
                                        [persp,y_fit]=persp_xform_find(adj,ref_tofit,opts_persp);
                                        a=persp(1:end-1,1:end-1);
                                        b=persp(end,1:end-1);
                                        c=persp(1:end-1,end);
                                end
                                adj_fit=persp_apply(a,b,c,adj);
                                d_den_shuff=sum(sum((ref_tofit-repmat(mean(ref_tofit,1),npts,1)).^2,1));
                                d_shuff_afpe(iafpe,ishuff)=sum((ref_tofit(:)-adj_fit(:)).^2)/d_den_shuff;
                            end %iafpe
                        end %ishuff
                        for iafpe=1:2
                            disp(sprintf('fitting from scratch with %12s but shuffled residuals: d>=shuffled value in %3.0f of %5.0f shuffles',...
                                fit_types{iafpe},sum(double(d_afpe(1,iafpe)>=d_shuff_afpe(iafpe,:))),nshuff_resids));
                            results.(fit_types{iafpe}).d_shuff(ref_ptr,adj_ptr,idim,:)=reshape(d_shuff_afpe(iafpe,:),[1 1 1 nshuff_resids]);
                        end
                    end
                    %
                    % plots
                    %
                    if ((idim==2) | (idim==3)) & if_pairplots
                        figure;
                        set(gcf,'Position',[50 50 1000 800]);
                        set(gcf,'NumberTitle','off');
                        set(gcf,'Name',sprintf('ref subj %s, adj subj %s, dim %2.0f',subj_ids{ref_ptr},subj_ids{adj_ptr},idim));
                        switch idim
                            case 2
                                hr=plot(ref(:,1),ref(:,2),'k.','MarkerSize',8);
                                hold on;
                                hac=plot(adj_check(:,1),adj_check(:,2),'x','MarkerSize',6);
                                hap=plot(adj_proc(:,1),adj_proc(:,2),'.','MarkerSize',8);
                                har1=plot(adj_xform{1}(:,1),adj_xform{1}(:,2),'.','MarkerSize',8);
                                har2=plot(adj_xform{2}(:,1),adj_xform{2}(:,2),'.','MarkerSize',8);
                                for k=1:npts
                                    hpc=plot([ref(k,1) adj_proc(k,1)],[ref(k,2) adj_proc(k,2)]);
                                    hp1=plot([ref(k,1) adj_xform{1}(k,1)],[ref(k,2) adj_xform{1}(k,2)]);
                                    hp2=plot([ref(k,1) adj_xform{2}(k,1)],[ref(k,2) adj_xform{2}(k,2)]);
                                    set(hpc,'Color',color_check);
                                    set(hp1,'Color',color_affine);
                                    set(hp2,'Color',color_projective);
                                end
                            case 3
                                hr=plot3(ref(:,1),ref(:,2),ref(:,3),'k.','MarkerSize',8);
                                hold on;
                                hac=plot3(adj_check(:,1),adj_check(:,2),adj_check(:,3),'x','MarkerSize',6);
                                hap=plot3(adj_proc(:,1),adj_proc(:,2),adj_proc(:,3),'.','MarkerSize',8);
                                har1=plot3(adj_xform{1}(:,1),adj_xform{1}(:,2),adj_xform{1}(:,3),'.','MarkerSize',8);
                                har2=plot3(adj_xform{2}(:,1),adj_xform{2}(:,2),adj_xform{2}(:,3),'.','MarkerSize',8);
                                for k=1:npts
                                    hpc=plot3([ref(k,1) adj_proc(k,1)],[ref(k,2) adj_proc(k,2)],[ref(k,3) adj_proc(k,3)]);
                                    hp1=plot3([ref(k,1) adj_xform{1}(k,1)],[ref(k,2) adj_xform{1}(k,2)],[ref(k,3) adj_xform{1}(k,3)]);
                                    hp2=plot3([ref(k,1) adj_xform{2}(k,1)],[ref(k,2) adj_xform{2}(k,2)],[ref(k,3) adj_xform{2}(k,3)]);
                                    set(hpc,'Color',color_check);
                                    set(hp1,'Color',color_affine);
                                    set(hp2,'Color',color_projective);
                                end
                                view(3);
                                zlabel('dim 3');
                        end
                        set(hr,'Color',color_data);
                        set(hap,'Color',color_procrustes);
                        set(hpc,'Color',color_check);
                        set(har1,'Color',color_affine);
                        set(har2,'Color',color_projective');
                        xlabel('dim 1');
                        ylabel('dim 2');
                        axis equal;
                        legend([hr;hap;hac;har1;har2],{'ref','adjusted (Procrustes)','check (Procrustes)','affine','projective'});
                        box on;
                        grid on;
                        title(sprintf('d (procrustes, regression, projective) %7.4f %7.4f %7.4f',d,d_afpe(1,:)));
                        %
                        axes('Position',[0.01,0.01,0.01,0.01]); %for text
                        text(0,0,tstring,'Interpreter','none','FontSize',10);
                        axis off;
                    end
                end %insufficient rank
            end %idim
        end %ijrev
    end %jset
end %iset
%
%summary plot
%
figure;
set(gcf,'Position',[100 100 1200 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name',cat(2,'Procrustes and regression: ',paradigm_name));
for ref_ptr=1:nsets
    for adj_ptr=1:nsets
        if (adj_ptr~=ref_ptr)
            subplot(nsets,nsets,adj_ptr+(ref_ptr-1)*nsets);
            hpp=plot([min_dim:max_dim],squeeze(results.procrustes.d(ref_ptr,adj_ptr,[min_dim:max_dim])),'LineWidth',line_width);
            set(hpp,'Color',color_procrustes);
            hold on;
            hpar1=plot([min_dim:max_dim],squeeze(results.affine.d(ref_ptr,adj_ptr,[min_dim:max_dim])),'LineWidth',line_width);
            set(hpar1,'Color',color_affine);
            hold on;
            hpar2=plot([min_dim:max_dim],squeeze(results.projective.d(ref_ptr,adj_ptr,[min_dim:max_dim])),'LineWidth',line_width);
            set(hpar2,'Color',color_projective);
            hold on;
            if nshuff_procrustes>0
                hppq=plot([min_dim:max_dim],squeeze(quantile(results.procrustes.d_shuff(ref_ptr,adj_ptr,:,:),nshuff_quantile,4)),':','LineWidth',line_width);
                set(hppq,'Color',color_procrustes);
            end
            if nshuff_resids>0
                hpar1q=plot([min_dim:max_dim],squeeze(quantile(results.affine.d_shuff(ref_ptr,adj_ptr,:,:),nshuff_quantile,4)),':','LineWidth',line_width);
                set(hpar1q,'Color',color_affine);
                hpar2q=plot([min_dim:max_dim],squeeze(quantile(results.projective.d_shuff(ref_ptr,adj_ptr,:,:),nshuff_quantile,4)),':','LineWidth',line_width);
                set(hpar2q,'Color',color_projective);
            end
            if (ref_ptr==1) & (adj_ptr==2)
                legend([hpp;hpar1;hpar2],{'Procrustes','affine','projective'},'FontSize',7,'Location','NorthEast');
            end
            set(gca,'XLim',[0.5 max_dim+0.5]);
            set(gca,'YLim',[0 1]);
            set(gca,'XTick',[1:max_dim]);
            set(gca,'XTickLabel',[1:max_dim]);
            %set(gca,'YLim',[0 min(1,1.1*max(results.procrustes.d(:)))]);
            set(gca,'YLim',[0 1]);
            set(gca,'YTick',[0:.25:1]);
            xlabel(cat(2,'dims     adj: ',subj_ids{adj_ptr}),'FontSize',9);
            ylabel(cat(2,'ref: ',subj_ids{ref_ptr}),'FontSize',9);
            title('d (unexplained var)','FontSize',9);
        end
    end
end
axes('Position',[0.01,0.01,0.01,0.01]); %for text
text(0,0,paradigm_name,'Interpreter','none','FontSize',10);
axis off;

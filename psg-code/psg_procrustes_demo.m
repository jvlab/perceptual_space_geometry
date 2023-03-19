%psg_procrustes_demo: demonstrate visualization of psg coordinates
%
% Uses Procrustes method to compare one or more sets of coordinates; can be experimental or predicted
%
% Compares each coordinate set with itself (as more dimensions are added to model)
%   Note that for coordinates from a quadratic form model, since each dimension is orthogonal to lower
%   ones, this analysis produces a trivial rotation
% Compares one coordinate set with another via classic (pairwise) Procrustes analysis on all datasets, across all choices of
%    model dimension in each, with and without scaling; producing results_procrustes.
%    Also plots goodness of fit and scaling as function of dimension and interactively
%    plots two or more datasets together
%
% Main results are in results_procrustes{iset,jset,if_scaling}.
%
%  See also: PSG_GET_COORDSETS, PSG_FINDRAYS, PSG_QFORMPRED, PSG_PLOTCOORDS,PSG_VISUALIZE_DEMO, PROCRUSTES,
%    PSG_COLORS_LEGACY.
%
if ~exist('opts_plot') opts_plot=struct(); end %for psg_plotcoords
if ~exist('opts_vis') opts_vis=struct(); end %for psg_visualize
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_qpred') opts_qpred=struct(); end %for psg_qformpred
opts_qpred=filldefault(opts_qpred,'qform_datafile_def','../stim/btc_allraysfixedb_avg_100surrs_madj.mat');
opts_qpred=filldefault(opts_qpred,'qform_modeltype',12); %if_symm=1 (symmetrize around origin), if_axes=1 (symmetrize bc, de, tuvw); ifaug=1 (augmented coords)
if ~exist('yplot_range_log') yplot_range_log=10.^[-2 2]; end %to prevent log problems
if ~exist('dfmt') dfmt='%s, [md: %1.0f]'; end %label formatting for Procrustes combined plots
%
opts_plot=psg_colors_legacy(opts_plot);
%
[sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred);
nsets=length(sets);
nstims=size(ds{1}{1},1);
dim_list_all=[];
for iset=1:nsets
    dim_list_all=union(dim_list_all,sets{iset}.dim_list);
end
%
% below here, code differs c/w psg_consensus_demo.
%
%computations for pairwise mode: one set with itself and pairwise analyses
%
results_procrustes=cell(nsets,nsets,2); %d1: reference set, d2: adjustable set, d3: allowing scaling or not
%note that results_procrustes{iset,jset,1}.d(i,j)=results_procrustes{jset,iset,1}.d(j,i)
for iset=1:nsets
    for jset=1:nsets
        disp(' ');
        disp(sprintf(' analyzing  with reference  set %2.0f: dim range [%3.0f %3.0f] label: %s',iset,min(sets{iset}.dim_list),max(sets{iset}.dim_list),sets{iset}.label));
        disp(sprintf(' analyzing  with adjustable set %2.0f: dim range [%3.0f %3.0f] label: %s',jset,min(sets{jset}.dim_list),max(sets{jset}.dim_list),sets{jset}.label));
        for if_scaling=1:2 %1: allow scaling, 2: not
            switch if_scaling
                case 1
                    allow_scaling=true;
                case 2
                    allow_scaling=false;
            end
            for ref_ptr=1:length(sets{iset}.dim_list)
                ref_dim=sets{iset}.dim_list(ref_ptr);
                for adj_ptr=1:length(sets{jset}.dim_list)
                    adj_dim=sets{jset}.dim_list(adj_ptr);
                    %
                    ref_coords=ds{iset}{ref_dim};
                    adj_coords=ds{jset}{adj_dim};
                    if ref_dim<adj_dim
                        ref_coords=[ref_coords,zeros(nstims,adj_dim-ref_dim)];
                    end
                    [d,z,transform]=procrustes(ref_coords,adj_coords,'Scaling',allow_scaling);
                    results_procrustes{iset,jset,if_scaling}.ref_label=sets{iset}.label;
                    results_procrustes{iset,jset,if_scaling}.adj_label=sets{jset}.label;
                    results_procrustes{iset,jset,if_scaling}.d(ref_dim,adj_dim)=d;
                    results_procrustes{iset,jset,if_scaling}.adj_transformed{ref_dim,adj_dim}=z;
                    results_procrustes{iset,jset,if_scaling}.scaling(ref_dim,adj_dim)=transform.b;
                    results_procrustes{iset,jset,if_scaling}.orthog{ref_dim,adj_dim}=transform.T;
                    results_procrustes{iset,jset,if_scaling}.translation{ref_dim,adj_dim}=transform.c(1,:);
                    results_procrustes{iset,jset,if_scaling}.transform='adj_transformed=scaling*adj_coords*orthog+repmat(translation,nstims,1) is closest fit to ref_coords';
                end %adj_ptr
            end %ref_ptr
        end %if_scaling
    end %jset 
end %iset
%
%simple plots
%
vplot_list={'d','scaling'};
ylim_list={[0 2],10.^[-1 1]};
ytickvals_list={[0 .5 1 2],10.^[-1:.5:1]};
if_logplot_list=[0 1];
%
%for each set, show effect of adding one more dimension
%
figure;
set(gcf,'Position',[50 50 1200 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name','Procrustes analysis within set');
for ifig=1:2
    vplot=vplot_list{ifig};
    ylims=ylim_list{ifig};
    ytickvals=ytickvals_list{ifig};
    if_logplot=if_logplot_list(ifig);
    ncol=max(nsets,4);
    for idir=1:2 %idir=1: ref is one dimension lower; idir=2: ref is one dimension higher
        for iset=1:nsets
            dims_both=intersect(sets{iset}.dim_list,1+sets{iset}.dim_list); %make sure we have idim and idim+1
            subplot(4,ncol,iset+((ifig-1)*2+idir-1)*ncol);
            plotvals=zeros(length(dims_both),2); %d2 is for scaling true vs scaling false
            for if_scaling=1:2
                rp=results_procrustes{iset,iset,if_scaling};
                for idimptr=1:length(dims_both)
                    idim=dims_both(idimptr);
                    if (idir==1)
                        plotvals(idimptr,if_scaling)=rp.(vplot)(idim-1,idim);
                        dims_plot=dims_both-1;
                    else
                        plotvals(idimptr,if_scaling)=rp.(vplot)(idim,idim-1);
                        dims_plot=dims_both;
                    end %idir
                end
            end
            if (if_logplot)
                plotvals=max(min(plotvals,yplot_range_log(2)),yplot_range_log(1));
                semilogy(dims_plot,plotvals);
            else
                plot(dims_plot,plotvals);
            end
            hold on;
            set(gca,'XLim',[-0.5 0.5]+[min(dim_list_all) max(dim_list_all)]);
            plot(get(gca,'XLim'),[1 1],'k:');
            set(gca,'XTick',dims_plot);
            set(gca,'YLim',ylims);
            set(gca,'YTick',ytickvals);
            xlabel('dim');
            ylabel(vplot);
            title(sprintf('set %1.0f',iset));
            legend({'scale','noscale'},'FontSize',7)
        end %iset
    end %idir
end
%
%compare one set with another
%
for ifig=1:2
    vplot=vplot_list{ifig};
    ylims=ylim_list{ifig};
    ytickvals=ytickvals_list{ifig};
    if_logplot=if_logplot_list(ifig);
    fig_title=sprintf('Procrustes analysis between sets: %s',vplot);
    %
    figure;
    set(gcf,'Position',[50 50 1200 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',fig_title);
    for iset=1:nsets
        for jset=1:nsets
            dims_both=intersect(sets{iset}.dim_list,sets{jset}.dim_list);
            subplot(nsets,nsets,jset+(iset-1)*nsets);
            plotvals=zeros(length(dims_both),2); %d2 is for scaling true vs scaling false
            for if_scaling=1:2
                rp=results_procrustes{iset,jset,if_scaling};
                for idimptr=1:length(dims_both)
                    idim=dims_both(idimptr);
                    plotvals(idimptr,if_scaling)=rp.(vplot)(idim,idim);
                end
            end
            if (if_logplot)
                plotvals=max(min(plotvals,yplot_range_log(2)),yplot_range_log(1));
                semilogy(dims_both,plotvals);
            else
                plot(dims_both,plotvals);
            end
            hold on;
            set(gca,'XLim',[-0.5 0.5]+[min(dim_list_all) max(dim_list_all)]);
            plot(get(gca,'XLim'),[1 1],'k:');
            set(gca,'XTick',dims_both);
            set(gca,'YLim',ylims);
            set(gca,'YTick',ytickvals);
            xlabel('dim');
            ylabel(vplot);
            title(sprintf('ref %1.0f adj %1.0f',iset,jset));
            legend({'scale','noscale'},'FontSize',7)
        end %jset
    end %iset
    axes('Position',[0.01,0.02,0.01,0.01]); %for text
    text(0,0,fig_title,'Interpreter','none');
    axis off;
end
%
%if requested, plot a reference dataset and a second set aligned to it
%
if_ok=0;
while (if_ok==0)
     disp('datasets available:');
    for iset=1:nsets
        disp(sprintf(' set %2.0f: dim range [%3.0f %3.0f] label: %s',iset,min(sets{iset}.dim_list),max(sets{iset}.dim_list),sets{iset}.label));
    end
    if_ok=1-getinp('1 to plot a reference dataset and a second set aligned to it','d',[0 1]);
    if (if_ok==0)
        iset=getinp('reference dataset','d',[1 nsets]);
        jset=getinp('dataset to align','d',[1 nsets]);
        dims_avail=intersect(sets{iset}.dim_list,sets{jset}.dim_list);
        if ~isempty(dims_avail)
            idim=getinp('reference dataset model dimension','d',[min(sets{iset}.dim_list),max(sets{iset}.dim_list)]);
            jdim=getinp('aligned dataset model dimension','d',[min(sets{jset}.dim_list),max(sets{jset}.dim_list)]);
             %set up library of what can be plotted
            lib=cell(1,4);
            %reference dataset
            lib{1}.name='reference dataset';
            lib{1}.d=ds{iset}{idim};
            lib{1}.sa=sas{iset};
            lib{1}.rays=rayss{iset};
            lib{1}.opts_vis=setfield(opts_vis,'file_string',sprintf(dfmt,sets{iset}.label,idim));
            %dataset to be adjusted
            lib{2}.name='dataset to be adjusted';
            lib{2}.d=ds{jset}{jdim};
            lib{2}.sa=sas{jset};
            lib{2}.rays=rayss{jset};
            lib{2}.opts_vis=setfield(opts_vis,'file_string',sprintf(dfmt,sets{jset}.label,jdim));
            %dataset after adjustment, with scaling
            lib{3}.name='dataset adjusted with scaling';
            lib{3}.d=results_procrustes{iset,jset,1}.adj_transformed{idim,jdim};
            lib{3}.sa=sas{jset};
            lib{3}.rays=rayss{jset};
            lib{3}.opts_vis=setfield(opts_vis,'file_string',sprintf(dfmt,cat(2,sets{jset}.label,sprintf('->Proc(scaled, %4.2f)',results_procrustes{iset,jset,1}.scaling(idim,jdim))),jdim));
            %dataset after adjustment, no scaling
            lib{4}.name='dataset adjusted without scaling';
            lib{4}.d=results_procrustes{iset,jset,2}.adj_transformed{idim,jdim};
            lib{4}.sa=sas{jset};
            lib{4}.rays=rayss{jset};
            lib{4}.opts_vis=setfield(opts_vis,'file_string',sprintf(dfmt,cat(2,sets{jset}.label,'->Proc (no scaling)'),jdim));
            %
            %ensure that the adjusted dataset has the same number of dimensions
            %
            for il=2:4
                if (jdim>idim)
                    lib{il}.d=lib{il}.d(:,1:idim);
                end
                if (jdim<idim)
                    lib{il}.d=[lib{il}.d zeros(size(lib{il}.d,1),idim-jdim)];
                end
            end
            if idim<=2
                plotformats=[idim 2];
            elseif idim==3
                plotformats=[idim 2;idim 3];
            elseif idim==4
                plotformats=[idim 2;idim 3;idim 4];
            else
                plotformats=[idim 3];
            end
            if_ok_plot=0;
            while if_ok_plot==0
                disp('plot choices');
                for il=1:length(lib)
                    disp(sprintf('%2.0f->%35s (%s)',il,lib{il}.name,lib{il}.opts_vis.file_string));
                end
                lib_list=getinp('choice(s), 0 to end and choose other datasets or model dimensions','d',[0 length(lib)]);
                if min(lib_list)==0
                    if_ok_plot=1;
                else
                    dm=cell(0);
                    sam=cell(0);
                    raysm=cell(0);
                    opts_vism=cell(0);
                    for il=1:length(lib_list)
                        dm{il}{1}=lib{lib_list(il)}.d;
                        sam{il}=lib{lib_list(il)}.sa;
                        raysm{il}=lib{lib_list(il)}.rays;
                        opts_vism{il}=lib{lib_list(il)}.opts_vis;
                    end
                    opts_plotm=struct;
                    if isfield(opts_plot,'colors')
                        opts_plotm.colors=opts_plot.colors;
                    end
                    opts_multm=struct;
                    [opts_vism_used,opts_plotm_used,opts_multm_used]=psg_visualize(plotformats,dm,sam,raysm,opts_vism,opts_plotm,opts_multm);
                end
            end %if_ok_plot
        end %available dimensions
    end %ok to plot
end %if_ok

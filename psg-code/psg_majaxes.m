function [results,opts_used]=psg_majaxes(d_ref,sa_ref,d_adj,sa_adj,results_geo,opts)
% [results,opts_used]=psg_majaxes(d_ref,sa_ref,d_adj,sa_adj,results_geo,opts) analyzes
%     an affine geometric model to determine major axes and plots
%     projections onto these axes
%
%  d_ref, sa_ref: data and sa (labeling) structure for reference dataset, typically from psg_read_coorddata
%  d_adj, sa_adj: data and sa (labeling) structure for adjustable dataset, typically from psg_read_coorddata
%       sa_ref,sa_adj: only relevant field is typenames; typically should match but this is not checked
%  results_geo: results field from psg_geomodels_run
%  opts: options (can be empty)
%     opts.if_log: 1 to log
%     opts.model_class_list: list of model classes, defaults to {'affine','pwaffine'};
%        should also run for {'affine','procrustes'} but un-interesting
%        since eigenvalues should all be 1, and eigenvectors are degenerate
%     opts.tol_neg: tolerance for negative eigenvalues
%     opts.tol_match: tolerance for matching eigenvalues
%     opts.plot_pairs: size [nplots 2], each row is [ref dim, adj dim] to plot, defaults to zeros(0,2)
%     opts.plot_colormap: color map for plots, defaults to 'jet'
%     opts.plot_submean: 1 to subtract means from plots
%     opts.plot_flipsign: 1 to allow sign of projections to be flipped to best match across piecewise transforms
%
%   consistency of results_geo (i.e., same models for each dimension) is not checked
%
%  results: analysis results
%  opts_used: options used
%    opts_used.figh: figure handles, if opts.plot_pairs is non-empty 
%      figh is figh(iplot,im_ptr): one plot for each geometric model analyzed)
%
%   See also: PSG_GEOMODELS_RUN, PSG_GEOMODELS_DEFINE.
%
if nargin<=5 opts=struct; end
%
opts=filldefault(opts,'if_log',0);
opts=filldefault(opts,'model_class_list',{'affine','pwaffine'});
opts=filldefault(opts,'tol_neg',10^-6);
opts=filldefault(opts,'tol_match',10^-4);
opts=filldefault(opts,'plot_pairs',zeros(0,2));
opts=filldefault(opts,'plot_colormap','jet');
opts=filldefault(opts,'plot_submean',1);
opts=filldefault(opts,'plot_flipsign',1);
%
model_types_def=psg_geomodels_define();
opts.model_types_def=model_types_def;
opts.warnings=[];
opts.figh=[];
%
results=cell(size(results_geo));
%
%scan the geometry results
nds_ref=size(results_geo,1);
nds_adj=size(results_geo,2);
d_ref_list=[];
d_adj_list=[];
for id_ref=1:nds_ref
    for id_adj=1:nds_adj
        res=results_geo{id_ref,id_adj};
        if ~isempty(res)
            d_ref_list=[d_ref_list,res.ref_dim];
            d_adj_list=[d_adj_list,res.adj_dim];
            model_types=res.model_types_def.model_types;
            ref_file=res.ref_file;
            adj_file=res.adj_file;
            results{id_ref,id_adj}.ref_dim=res.ref_dim;
            results{id_ref,id_adj}.adj_dim=res.adj_dim;
            results{id_ref,id_adj}.ref_file=ref_file;
            results{id_ref,id_adj}.adj_file=adj_file;
        end
    end %id_adj
end %id_ref
d_ref_list=unique(d_ref_list);
d_adj_list=unique(d_adj_list);
if opts.if_log
    disp(sprintf('reference dataset: %s',ref_file));
    disp('dimensions available:');
    disp(d_ref_list);
    disp('stimuli');
    disp(sa_ref.typenames);
    disp(sprintf('adjusted dataset: %s',adj_file));
    disp('dimensions available:');
    disp(d_adj_list);
    disp('stimuli');
    disp(sa_adj.typenames);
    disp('model types used for fitting geometric transformation from adjusted to reference');
    disp(model_types);
    disp('model classes to be anayzed');
    disp(opts.model_class_list);
end
%calculations
adj_ref_labels={'adj','ref'};
nar=length(adj_ref_labels);
%
im_ptr=0;
for im=1:length(model_types)
    model_type=model_types{im};
    mclass=strmatch(model_types_def.(model_type).class,opts.model_class_list,'exact');
    if mclass>0 
        im_ptr=im_ptr+1;
        if opts.if_log
            disp(sprintf('analyzing model %25s (class: %s)',model_type,opts.model_class_list{mclass}))
        end   
        for id_ref=1:nds_ref
            for id_adj=1:nds_adj
                if ~isempty(res)
                    ref_dim=results{id_ref,id_adj}.ref_dim;
                    adj_dim=results{id_ref,id_adj}.adj_dim;
                    results{id_ref,id_adj}.model_types{im_ptr}=model_type;
                    transform_struct=results_geo{id_ref,id_adj}.transforms{im};
                    results{id_ref,id_adj}.ref.typenames=sa_ref.typenames;
                    results{id_ref,id_adj}.adj.typenames=sa_adj.typenames;
                    b=transform_struct.b;
                    Tlist=transform_struct.T;
                    npw=size(Tlist,3);
                    for ipw=1:npw
                        T=Tlist(:,:,ipw); %size of T is adj_dim x max(adj_dim,ref_dim),
                        %and if adj_dim> ref_dim, then columns ref_dim+1:end should be zero
                        num_eigs=adj_dim; % this is min(size(T)), number of eigenvalues computed
                        num_eigs_nz=min(ref_dim,adj_dim); %expected number of nonzero eigenvalues
                        for iar=1:nar
                            lab=adj_ref_labels{iar};
                            switch iar
                                case 1
                                    % lab='adj'; %compute eigenvals and eigenvecs of T*Ttranspose
                                    A=T*transpose(T); %think of T as a stretching matrix M * a rotation matrix R
                                    d_coords=d_adj{adj_dim}; %coordinates in dataset that is adjusted
                                    %A is square, size is adj_dim
                                case 2
                                    % lab='ref'; %compute eigenvals and eigenvecs of Ttranspose*T
                                    A=transpose(T)*T; %think of T as a rotation matrix R * a stretching matrix M
                                    d_coords=d_ref{ref_dim}; %coordinates in reference dataset
                                    %A is square, size is max(adj_dim,ref_dim)
                            end
                            [eivecs,eivals,opts]=psg_majaxes_eigs(A,sprintf('%s [%s ipw %1.0f]',lab,model_type,ipw),ref_dim,adj_dim,opts);
                            results{id_ref,id_adj}.(lab).magnifs{im_ptr}(:,ipw)=b*sqrt(eivals); % since A=(MR)*transpose(MR)
                            results{id_ref,id_adj}.(lab).eivecs{im_ptr}(:,:,ipw)=eivecs;
                            results{id_ref,id_adj}.(lab).magnif_ratio{im_ptr}(:,ipw)=sqrt(eivals(1)/eivals(num_eigs_nz)); %ratio of highest to lowest magnification factor
                            prj_dim=size(d_coords,2); %number of dimensions to project onto, less than size(A) if ref_dim<adj_dim lab='ref'
                            %disp(sprintf('ref_dim %2.0f adj_dim %2.0f lab %s prj_dim %1.0f size(A) %2.0f %2.0f',ref_dim,adj_dim,lab,prj_dim,size(A)));
                            eivecs_project=eivecs(1:prj_dim,1:prj_dim); %if id_ref<id_adj, remaining eigenvecs are units
                            %rows of d_coords are the coordinates for each stimulus, either in adj set or ref set
                            results{id_ref,id_adj}.(lab).projections{im_ptr}(:,:,ipw)=d_coords*eivecs_project; %projections of original coordinates onto eigenvectors
                        end %adj or ref
                        %check that magnification factors agree (sqrt of eigenvals)
                        eig_diff=max(abs(results{id_ref,id_adj}.ref.magnifs{im_ptr}(1:num_eigs)-results{id_ref,id_adj}.adj.magnifs{im_ptr}(1:num_eigs)));
                        if eig_diff>opts.tol_match
                            warning_text=sprintf('magnifications disagree [%s ipw %1.0f] ref_dim %2.0f adj_dim %2.0f, disparity is %15.12f',model_type,ipw,ref_dim,adj_dim,eig_diff);
                            if opts.if_log
                                disp(warning_text);
                            end
                            opts.warnings=strvcat(opts.warnings,warning_text);
                        end
                    end %ipw: each transformation matrix
                end %non-empty
            end %id_adj
        end %id_ref
    end %model type is affine
end
models_analyzed=im_ptr;
%plots
nplots=size(opts.plot_pairs,1);
opts.figh=cell(nplots,models_analyzed);
for iplot=1:nplots
    ref_plot=opts.plot_pairs(iplot,1);
    adj_plot=opts.plot_pairs(iplot,2);
    id_ref=min(find(d_ref_list==ref_plot));
    id_adj=min(find(d_adj_list==adj_plot));
    if ~isempty(id_ref) & ~isempty(id_adj)
        res_plot=results{id_ref,id_adj};
        for im_ptr=1:models_analyzed
            fig_label=sprintf(' ref dim %2.0f  adj dim %2.0f  model %s',ref_plot,adj_plot,res_plot.model_types{im_ptr});
            opts.figh{iplot,im_ptr}=figure;
            set(gcf,'Position',[100 100 1200 800]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',fig_label);
            npw=size(res_plot.ref.eivecs{im_ptr},3);
            nrows=max(2,npw);
            if (opts.if_log)
                disp(sprintf('plotting %s',fig_label));
            end
            for iar=1:nar
                for ipw=1:npw
                    lab=adj_ref_labels{iar};
                    subplot(nrows,nar+1,(nar+1)*(ipw-1)+iar);
                    z=res_plot.(lab).projections{im_ptr}(:,:,ipw); %projections
                    nstims=size(z,1);
                    neivs=size(z,2);
                    if opts.plot_submean
                        z=z-repmat(mean(z,1),nstims,1);
                    end
                    flips=cell(1,neivs);
                    if opts.plot_flipsign
                        if (ipw==1)
                            z1=z;
                        else
                            %if more than one transform, flip signs for the later transforms so that there
                            %is the closest match to the first transform
                            zdots=z'*z1;
                            for ieiv=1:neivs
                                meiv=find(abs(zdots(ieiv,:))==max(abs(zdots(ieiv,:))));
                                sgn=sign(zdots(ieiv,meiv));
                                z(:,ieiv)=z(:,ieiv)*sgn;
                                if sgn<0
                                    flips{ieiv}='-';
                                else
                                    flips{ieiv}='+';                              
                                end
                            end
                        end
                    end
                    imagesc(z',max(abs(z(:)))*[-1 1]);
                    set(gca,'XTick',[1:nstims]);
                    set(gca,'XTickLabel',res_plot.(lab).typenames);
                    set(gca,'YTick',[1:neivs]); %label each eigenvector with magnif factor
                    ytl=cell(1,neivs);
                    for ieiv=1:neivs
                        ytl{ieiv}=sprintf('%1s eiv %1.0f: x %5.2f',flips{ieiv},ieiv,res_plot.(lab).magnifs{im_ptr}(ieiv));
                    end
                    set(gca,'YTickLabel',ytl);
                    colormap(opts.plot_colormap);
                    title(sprintf('%s, transform %1.0f',lab,ipw));
                end %ipw
            end %adj or ref
            subplot(1,nar+1,nar+1);
            hc=colorbar;
            set(hc,'TickLabelsMode','manual');
            set(hc,'Ticks',[0:.1:1]);
            set(hc,'TickLabels',[-1:.2:1]);
            axis off;     
            %
            axes('Position',[0.01,0.02,0.01,0.01]); %for text
            text(0,0,fig_label,'Interpreter','none','FontSize',8);
            axis off;
            axes('Position',[0.01,0.04,0.01,0.01]); %for text
            text(0,0,cat(2,'ref: ',ref_file),'Interpreter','none','FontSize',8);
            axis off;
            axes('Position',[0.01,0.06,0.01,0.01]); %for text
            text(0,0,cat(2,'adj: ',adj_file),'Interpreter','none','FontSize',8);
            axis off;
        end %im_ptr
    end %dims not avail
end %plot
%
opts_used=opts;
return

function [eivecs,eivals,opts_new]=psg_majaxes_eigs(A,label,ref_dim,adj_dim,opts)
%compute, sort, and check eigenvalues and eigenvecs
[eivecs,eivals]=eig(A);
eivals=real(diag(eivals)); %A is self-adjoint
if any(eivals<-opts.tol_neg)
    warning_text=sprintf('negative eigenvalue in %s set to zero for ref_dim %2.0f adj_dim %2.0f: %15.12f',label,ref_dim,adj_dim,min(eivals));
    if opts.if_log
        disp(warning_text);
    end
    opts.warnings=strvcat(opts.warnings,warning_text);
end
opts_new=opts;
eivals=max(eivals,0);
[eivals,sort_inds]=sort(eivals,'descend'); %obtain eigenvalues in descending order
eivecs=real(eivecs(:,sort_inds));
return

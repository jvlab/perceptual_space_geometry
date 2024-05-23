%psg_curv_loglik_plot: read a file (typically csv) of log likelihoods as
%functoin of curvature and dimension, and plot
%
% data files from Suniyya Waraich
%
%   See also:  PSG_COLORS_LIKE.
%
if ~exist('domain_name') domain_name='bc6pt9'; end
if ~exist('filename_base') filename_base='psg_data/*_curved_model_log_likelihoods_debiased.csv'; end
if ~exist('dim_colors') dim_colors={'k','b','m','r','g'}; end
%
plotvar_names_alt={'LogLikelihood','BiasEstimate','CorrectedLLRelativeToBestModel','CorrectedLLs'};
if ~exist('curv_name') curv_name='CurvatureOfSpace'; end %alternative could be Lambda_Mu
%
domain_name=getinp('experiment name','s',[],domain_name);
filename=strrep(filename_base,'*',domain_name);
filename=getinp('log likelihood file to read','s',[],filename);
T=readtable(filename);
%
sigmas=unique(T{:,'Sigma'});
if ~all(sigmas==1)
    warning('some sigmas are not 1');
    sigmas
end
colors_like=psg_colors_like;
subj_symbs_res=colors_like.subj_symbs_res; %reserved symbols
subj_fills_res=colors_like.subj_fills_res; %reserved fills
subj_symbs_unres='ovsdph<>h'; %other available symbols
%
if ~exist('plotvar_name') plotvar_name='CorrectedLLs';end
for ipv=1:length(plotvar_names_alt)
    disp(sprintf('%1.0f ->%s',ipv,plotvar_names_alt{ipv}));
end
ipv=getinp('choice','d',[1 length(plotvar_names_alt)],strmatch(plotvar_name,plotvar_names_alt,'exact'));
plotvar_name=plotvar_names_alt{ipv};
%
% choose subjects and define symbols
%
subj_list_all=unique(T{:,'Subject'});
for isubj=1:length(subj_list_all)
    disp(sprintf('subject %1.0f->%s',isubj,subj_list_all{isubj}));
end
subj_list=getinp('choice(s)','d',[1 length(subj_list_all)],1:length(subj_list_all));
%remove the reserved symbols to create available symbols for other subjects
for isubj=1:length(subj_list)
    subj_id=subj_list_all{subj_list(isubj)};
    if isfield(subj_symbs_res,subj_id)
        subj_symbs_unres=strrep(subj_symbs_unres,subj_symbs_res.(subj_id),'');
    end
end
%
dim_list_all=unique(T{:,'Dimension'})';
dim_list=getinp('dimension(s)','d',[min(dim_list_all) max(dim_list_all)],dim_list_all);
dim_list=intersect(dim_list,dim_list_all);
%
curv_list=unique(T{:,curv_name})';
curv_range=getinp('curvature range','f',[min(curv_list) max(curv_list)],[min(curv_list) max(curv_list)]);
%
tstring=cat(2,'curvature analysis ',domain_name);
figure;
set(gcf,'Name','curvature')
set(gcf,'Position',[100 100 1000 800]);
set(gcf,'NumberTitle','off');
n_unres=0;
hl=cell(0);
ht=[];
for isubj=1:length(subj_list)
    subj_id=subj_list_all{subj_list(isubj)};
    if isfield(subj_symbs_res,subj_id)
        subj_symb=subj_symbs_res.(subj_id);
        subj_fill=subj_fills_res.(subj_id);
    else
        n_unres=n_unres+1;
        subj_symb=subj_symbs_unres(1+mod(n_unres-1,length(subj_symbs_unres)));
        subj_fill=1;
    end
    T_subj=T(strmatch(subj_id,T.Subject,'exact'),:);
    for idim=1:length(dim_list)
        dim_plot=dim_list(idim);
        dim_color=dim_colors{1+mod(idim-1,length(dim_colors))};
        T_plot=T_subj(find(T_subj.Dimension==dim_plot),:);
        xyplot=[table2array(T_plot(:,curv_name)),table2array(T_plot(:,plotvar_name))];
        %sort and censor by curvature range
        curv_select=intersect(find(xyplot(:,1)>=curv_range(1)),find(xyplot(:,1)<=curv_range(2)));
        xyplot=sortrows(xyplot(curv_select,:));
        if ~isempty(xyplot)
            hp=plot(xyplot(:,1),xyplot(:,2),subj_symb);
            set(hp,'Color',dim_color);
            hold on;
            hline=plot(xyplot(:,1),xyplot(:,2),'k-');
            set(hline,'Color',dim_color);
            hold on;
            if (subj_fill)
                set(hp,'MarkerFaceColor',dim_color);
            end
            hl=[hl,hp];
            ht=strvcat(ht,sprintf('%s dim %1.0f',subj_id,dim_plot));
        end
    end
end
plot([0 0],get(gca,'YLim'),'k:'); %a vertical line at 0
legend(hl,ht,'Location','best');
xlabel(curv_name);
ylabel(plotvar_name);
title(tstring);


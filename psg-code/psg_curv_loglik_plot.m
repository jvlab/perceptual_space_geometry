%psg_curv_loglik_plot: read a file (typically csv) of log likelihoods as
%functoin of curvature and dimension, and plot
%
% data files from Suniyya Waraich
%
%   See also:  PSG_COLORS_LIKE.
%
if ~exist('domain_name') domain_name='bc6pt9'; end
if ~exist('filename_base') filename_base='psg_data/*_curved_model_log_likelihoods_debiased.csv'; end
if ~exist('dcolors') dcolors={'k','b','m','r','g'}; end
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
        symb=subj_symbs_res.(subj_id);
        fill=subj_fills_res.(subj_id);
    else
        n_unres=n_unres+1;
        symb=subj_symbs_unres(1+mod(n_unres-1,length(subj_symbs_unres)));
        fill=1;
    end
    T_subj=T(strmatch(subj_id,T.Subject,'exact'),:);
    for idim=1:length(dim_list)
        dim_plot=dim_list(idim);
        T_plot=T_subj(find(T_subj.Dimension==dim_plot),:);
        xvals=table2array(T_plot(:,curv_name));
        yvals=table2array(T_plot(:,plotvar_name));
        hp=plot(xvals,yvals,symb);
        hl=[hl,hp];
        ht=strvcat(ht,sprintf('%s dim %1.0f',subj_id,dim_plot));
        hold on;

    end
end
legend(hl,ht,'Location','best');
xlabel(curv_name);
ylabel(plotvar_name);
title(tstring);


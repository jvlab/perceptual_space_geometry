% mater_ratings_survey_rs: read materials ratings files and analyze them
% 
% creates a mat-file with consolidated 'ratings' data or reads a pre-existing one
% 
if ~exist('path_ratings') path_ratings='C:\Users\jdvicto\Documents\JV\ey7977\Materials\FabricAffordanceRatings'; end
if ~exist('string_ratings') string_ratings='_participant_ratings.csv'; end
%
%section to force btc defaults, even if rs_aux_defaults.mat has been created or modified
if ~exist('aux_force_filename') aux_force_filename='rs_aux_defaults_btc.mat'; end
auxs_force=struct;
opts_needed={'opts_read','opts_rays','opts_check','opts_align','opts_import','opts_qpred','opts_knit','opts_write'};
for k=1:length(opts_needed)
    auxs_force.(opts_needed{k})=rs_aux_force(opts_needed{k},[],aux_force_filename);
end
%
if getinp('1 to create a new ratings mat-file (otherwise, load a file with ''ratings'')','d',[0 1],~exist('ratings','var'));
    if_ok=0;
    while(if_ok==0)
        if_ok=1;
        dir_cmd=cat(2,'dir ',path_ratings,'\*',string_ratings,' /b');
        [dos_status,dos_out]=dos(dir_cmd);
        if dos_status==0
            filenames_ratings_all=cell(0);
            while ~isempty(dos_out)
                [filenames_ratings_all{end+1},dos_out]=strtok(dos_out);
            end
            filenames_ratings_all=filenames_ratings_all(1:end-1);
            n_aff_all=length(filenames_ratings_all);
            names_aff_all=cell(1,n_aff_all);
            for iaff=1:n_aff_all
                names_aff_all{iaff}=filenames_ratings_all{iaff}(1:end-length(string_ratings));
                disp(sprintf('affordance %2.0f->%s',iaff,names_aff_all{iaff}))
            end
            aff_sel=getinp('choice(s)','d',[1 n_aff_all],[1:n_aff_all]);
            n_aff=length(aff_sel);
            names_aff=names_aff_all(aff_sel);
            filenames_ratings=filenames_ratings_all(aff_sel);
            %
            tables=cell(1,n_aff);
            %
            %read csv files and check for consistency
            for iaff=1:n_aff
                msgs=[];
                %warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
                tables{iaff}=readtable(cat(2,path_ratings,filesep,filenames_ratings{iaff}));
                %warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');
                disp(sprintf(' ratings for %30s read from %50s',names_aff{iaff},filenames_ratings{iaff}));
                %expect that first variable name is 'image', and subject names are followed by 'average'
                vns=tables{iaff}.Properties.VariableNames;
                if iaff==1
                    av_col=strmatch('average',vns,'exact');
                    subj_names=vns(2:av_col-1);
                    if length(av_col)~=1
                         msgs=strvcat(msgs,'cannot find ''average'' column');
                    end
                    image_nos=table2array(tables{1}(:,'image'));
                end
                av_col_match=min(strmatch('average',vns,'exact'),Inf);
                if av_col_match~=av_col
                    msgs=strvcat(msgs,'''average'' column not where expected');
                end
                image_col=min(strmatch('image',vns,'exact'),Inf);
                if image_col~=1
                    msgs=strvcat(msgs,'''image'' column not where expected');
                end
                subj_names_match=vns(2:av_col-1);
                if ~strcmp(subj_names_match,subj_names)
                    msgs=strvcat(msgs,'subject names do not match');
                end
                image_nos_match=table2array(tables{iaff}(:,'image'));
                if ~all(image_nos_match==image_nos)
                    msgs=strvcat(msgs,'image numbers do not match');
                end
                if ~isempty(msgs)
                    disp(msgs);
                    if_ok=0;
                end
            end
        else
            disp(dos_out);
            if_ok=0;
        end
        if if_ok==0
            disp('cannot proceed, try selecting a subset of ratings files.')
        end
    end
    %subject conversions
    subj_names=strrep(subj_names,'_s1','');
    subj_names=strrep(subj_names,'GRACE','YX');
    %
    %convert ratings data into a structure
    ratings=struct;
    n_imgs=length(image_nos);
    n_subjs=length(subj_names);
    ratings.n_subjs=n_subjs;
    ratings.n_imgs=n_imgs;
    ratings.n_aff=n_aff;
    ratings.names_aff_all=names_aff_all;
    ratings.names_aff=names_aff;
    ratings.aff_sel=aff_sel;
    ratings.subj_names=subj_names;
    ratings.image_nos=image_nos;
    ratings.filenames_ratings_all=filenames_ratings_all;
    ratings.filenames_ratings=filenames_ratings;
    ratings.vals=zeros(n_imgs,n_aff,n_subjs);
    ratings.avgs=zeros(n_imgs,n_aff);
    for iaff=1:n_aff
        ratings.vals(:,iaff,:)=reshape(table2array(tables{iaff}(:,2:av_col-1)),[n_imgs,1,n_subjs]);
        ratings.avgs(:,iaff)=table2array(tables{iaff}(:,'average'));
    end
    ratings.avgs_recompute=reshape(mean(ratings.vals,3),[n_imgs,n_aff]);
    disp(sprintf(' comparing averages and recomputed values: max diff is %5.3f',max(abs(ratings.avgs(:)-ratings.avgs_recompute(:)))));
    if getinp('1 to save the ratings data as a mat-file','d',[0 1])
        ratings_matfile=getinp('file name (suggest affordance_ratings*.mat)','s',[]);
        save(ratings_matfile,'ratings');
    end
end
%from here down, only use fields of 'ratings'
    % ratings=struct;
    % n_imgs=length(image_nos);
    % n_subjs=length(subj_names);
    % ratings.names_aff_all=names_aff_all;
    % ratings.names_aff=names_aff;
    % ratings.aff_sel=aff_sel;
    % ratings.subj_names=subj_names;
    % ratings.image_nos=image_nos;
    % ratings.filenames_ratings_all=filenames_ratings_all;
    % ratings.filenames_ratings=filenames_ratings;
    % ratings.vals=zeros(n_imgs,n_aff,n_subjs);
    % ratings.avgs=zeros(n_imgs,n_aff);

%create representational space of affordances from each subject and average 
aux_import=struct();
opts_import=struct;
typenames=cell(ratings.n_imgs,1);
for k=1:ratings.n_imgs
    typenames{k}=cat(2,'img_',zpad(ratings.image_nos(k),3));
end
opts_import.typenames=typenames;
for isubj=1:ratings.n_subjs
    [recon_pcaxes,recon_coords,var_ex,var_tot,coord_maxdiff,opts_used]=psg_pcaoffset(ratings.vals(:,:,isubj));
    disp(sprintf(' subj %2.0f (%10s): frac of variance explained',isubj,ratings.subj_names{isubj}));
    %%%%%need to capture pca matrics from opts_used
% opts_used = 
%   struct with fields:
% 
%     if_log: 0
%     nd_max: Inf
%         qu: [80×10 double]
%         qs: [10×10 double]
%         qv: [10×10 double]
%     offset: [0 0 0 0 0 0 0 0 0 0]
% 
    disp(diff([0 var_ex])/var_ex(end));
    opts_import.subj_id=ratings.subj_names{isubj};
    opts_import.paradigm_type='affordance ratings';
    opts_import.label_long=cat(2,'affordance ratings, ',ratings.subj_names{isubj});
    data_ratings=rs_import_coordsets(recon_pcaxes,setfield(aux_import,'opts_import',opts_import));
    if isubj==1
        data_ratings_all=data_ratings;
    else
        data_ratings_all=rs_concat_coordsets(data_ratings_all,data_ratings);
    end
end
%do a consensus
%plot individual spaces and consensus
%

% %
% if ~exist('tol') tol=10^-5; end %for matching stat calcs with matlab
% if ~exist('dmax') dmax=5; end
% if ~exist('dims_plot') dims_plot=[2 3];end
% if ~exist('line_width') line_width=2; end
% if ~exist('marker_size') marker_size=12; end
% %
% dmax=getinp('maximum dimension to analyze','d',[2 10],dmax);
% dims_plot=getinp('dimensions to plot (0 for none)','d',[0 dmax],dims_plot);
% dims_plot=setdiff(dims_plot,[0 1]);
% if ~isempty(dims_plot)
%     if_plain=getinp('1 to also plot spaces without color-by-behavior','d',[0 1],0);
% end
% %
% if ~exist('opts_read') opts_read=struct(); end
% opts_read.input_type=1; %just data
% opts_read.if_warn=0; %may have different sets of stimuli or stimuli in different orders
% opts_read.if_auto=1; %no confirmation needed
% opts_read.type_class_def='hlid';
% opts_read.type_coords_def='zeros';
% opts_read.type_class_aux=opts_read.type_class_def;
% opts_read.ui_filter='*coords*consensus*sf0*dff*pr-all*';
% opts_read.if_log=0;
% %
% if ~exist('opts_check') opts_check=struct(); end
% opts_check.if_warn=0;
% %
% if ~exist('opts_disp') opts_disp=struct(); end
% opts_disp.coord_group_method='keeplow';
% opts_disp.if_legend=0;
% opts_disp.callout_amount=0.5;
% opts_disp.set_markersizes=marker_size;
% opts_disp.set_colors='k';
% %
% aux=struct;
% aux.nsets=0;
% aux.opts_check=opts_check;
% aux.opts_read=opts_read;
% %
% fullnames=[];
% [data_read,aux_read_out]=rs_get_coordsets(fullnames,aux);
% nsets=length(data_read.sets);
% disp(sprintf(' %3.0f sets read.',nsets));
% mean_beh=mean(data_beh,2);
% range_mean_beh=[min(mean_beh),max(mean_beh)];
% if ~exist('colors_beh') colors_beh=[1 0 0;0 1 1]; end
% %
% std_beh=std(data_beh,0,2);
% range_std_beh=[min(std_beh),max(std_beh)];
% if ~exist('sizes_beh') sizes_beh=[30;15]; end %render std dev by symbol size (smaller std is bigger symbol)
% %
% for iset=1:nsets
%     typenames=data_read.sas{iset}.typenames;
%     nstims=length(typenames);
%     ptrs_beh=zeros(1,nstims);
%     for k=1:nstims_beh
%         ptr_beh=strmatch(typenames_beh{k},typenames,'exact');
%         if length(ptr_beh)~=1
%             ptr_beh=0;
%         end
%         ptrs_beh(k)=ptr_beh;
%     end
%     filename=strrep(strrep(data_read.sets{iset}.label_long,'/',filesep),'\',filesep);
%     filename_end=max(find(cat(2,filesep,filename)==filesep));
%     filename_short=filename(filename_end:end);
%     filename_short=strrep(filename_short,'.mat','');
%     disp(sprintf('analyzing %s',filename_short));   
%     if any(ptrs_beh==0)
%         disp(sprintf('no coords for %1.0f behaviors found in %s',sum(ptrs_beh==0),filename_short));
%     end
%     y=mean_beh(find(ptrs_beh>0)); %regress agains behaviors that have coordinates
%     pvals=zeros(dmax,1);
%     frats=zeros(dmax,1);
%     rsquareds=zeros(dmax,1);
%     rsquareds_drop=zeros(dmax,1);
%     rmse=zeros(dmax,1);
%     rmse_drop=zeros(dmax,1);
%     b=cell(dmax,1);
%     b_intvl=cell(dmax,1);
%     b_drop=cell(dmax,1);
%     sig_string=cell(dmax,1);
%     for dim=1:dmax
%         coords=data_read.ds{iset}{dim};
%         x=coords(ptrs_beh(ptrs_beh>0),:); %regress against behaviors that have coords
%         n=size(x,1);
%         [b{dim},b_intvl{dim},r,r_intvl,stats]=regress(y,[ones(n,1),x]); %add a constant term
%         disp(sprintf(' dim %1.0f: regressors (constant and each pc), and 0.95 confidence limits',dim))
%         disp([b{dim},b_intvl{dim}]')
%         for k=1:dim
%             if sign(b_intvl{dim}(k+1,1))==sign(b_intvl{dim}(k+1,2))
%                 sig_string{dim}=cat(2,sig_string{dim},sprintf(' dim %2.0f ',k));
%             else
%                 sig_string{dim}=cat(2,sig_string{dim},'        ');
%             end
%         end
%         %stats: the R-square statistic, the F statistic, p value for the full model, and an estimate of the error variance.
% %       disp(sprintf('   p=%6.4f, F=%8.4f, R^2=%6.4f, from stats',stats(3),stats(2),stats(1)));
%         %recalculate stats, first principles
%         y_pred=[ones(n,1),x]*b{dim};
%         ss_model=sum((y_pred-mean(y)).^2);
%         ss_error=sum((y-y_pred).^2);
%         frats(dim)=(ss_model/dim)/(ss_error/(n-dim-1));
%         pvals(dim)=1-fcdf(frats(dim),dim,n-dim-1);
%         rsquareds(dim)=corr(y,y_pred).^2;
% %       disp(sprintf('   p=%6.4f, F=%8.4f, R^2=%6.4f, recalc',p,frat,Rsquared));
%         %
%         if abs(pvals(dim)-stats(3))>tol
%             disp(sprintf('mismatch of p: %7.3f (matlab) vs %7.3f (recalc)',pvals(dim),stats(3)));
%         end
%         if abs(frats(dim)-stats(2))>tol
%             disp(sprintf('mismatch of F-ratio: %7.3f (matlab) vs %7.3f (recalc)',frats(dim),stats(2)));
%         end
%         if abs(rsquareds(dim)-stats(1))>tol
%             disp(sprintf('mismatch of R-squared: %7.3f (matlab) vs %7.3f (recalc)',rsquareds(dim),stats(1)));
%         end
%         y_pred_drop=zeros(n,1);
%         for k=1:n
%             i_drop=setdiff([1:n],n);
%             x_drop=x(i_drop,:);
%             y_drop=y(i_drop,:); 
%             b_drop{dim}=regress(y_drop,[ones(n-1,1),x_drop]);
%             y_pred_drop(k)=[1 x(k,:)]*b_drop{dim};
%         end
%         ss_error_drop=sum((y-y_pred_drop).^2);
%         rmse(dim)=sqrt(ss_error/n);
%         rmse_drop(dim)=sqrt(ss_error_drop/n);
%         rsquareds_drop(dim)=corr(y,y_pred_drop).^2;
%     end %dim
%     disp('   dim     p    f-ratio      R^2    rmse     R^2_drop rmse_drop   signif regressors');
%     for dim=1:dmax
%         disp(sprintf('%5.0f  %7.3f %7.4f    %7.4f %7.4f    %7.4f %7.4f         %s',dim,pvals(dim),frats(dim),rsquareds(dim),rmse(dim),rsquareds_drop(dim),rmse_drop(dim),sig_string{dim}))
%     end
%     std_beh_have=std(data_beh(ptrs_beh>0,:),0,2);
%     rms_beh=sqrt(mean(std_beh_have.^2));
%     disp(sprintf(' rms dev for behavior: %7.4f',rms_beh));
%     %
%     %plot
%     %
%     if ~isempty(dims_plot)
%         for dim_plot=dims_plot
%             opts_disp.set_select=1;
%             opts_disp.dim_select=dim_plot;
%             %
%             if if_plain
%                 %simple plot, all points black
%                 data_read_oneset=rs_extract_coordsets(data_read,iset);
%                 aux_disp=rs_disp_coordsets(data_read_oneset,setfield(aux,'opts_disp',opts_disp));
%                 set(gcf,'Name',sprintf('dim %1.0f, %s',dim_plot,filename_short));
%                 %
%                 axes('Position',[0.01,0.04,0.01,0.01]); %for text
%                 text(0,0,filename_short,'Interpreter','none','FontSize',8);
%                 axis off;
%             end
%             %
%             %plot with custom colors for each point: make each point a different set
%             %callouts will have slightly different lengths, since they are normalized by rms within each set
%             %
%             ptrs_beh_have=ptrs_beh(ptrs_beh>0);
%             %
%             data_indiv=struct;
%             data_indiv.sets=cell(1,n);
%             data_indiv.sas=cell(1,n);
%             data_indiv.ds=cell(1,n);
%             for istim_beh=1:n
%                 istim=ptrs_beh_have(istim_beh); %which stimulus, in list of typenames
%                 data_indiv.sets{istim_beh}=data_read.sets{iset};
%                 data_indiv.sets{istim_beh}.nstims=1;
%                 data_indiv.sas{istim_beh}=data_read.sas{iset};
%                 data_indiv.sas{istim_beh}.nstims=1;
%                 data_indiv.sas{istim_beh}.typenames=data_read.sas{iset}.typenames(istim);
%                 data_indiv.sas{istim_beh}.btc_specoords=data_read.sas{iset}.btc_specoords(istim,:);
%                 data_indiv.ds{istim_beh}=data_read.ds{iset};
%                 istim=ptrs_beh_have(istim_beh);
%                 for dq=1:length(data_indiv.ds{istim_beh})
%                     data_indiv.ds{istim_beh}{dq}=data_indiv.ds{istim_beh}{dq}(istim,:);
%                 end
%                 % istim_beh
%                 % istim
%                 % typenames(istim)
%             end
%             opts_disp_indiv=opts_disp;
%             opts_disp_indiv.data_label_setsel_method='all';
%             opts_disp_indiv.set_select=[1:n]; %each set is a different stimulus
%             opts_disp_indiv.data_label_typenames_vary=1;
%             opts_disp_indiv.callout_center=mean(data_read.ds{iset}{dim_plot},1);
%             opts_disp_indiv.callout_colors='set_colors';
%             opts_disp_indiv.set_colors=cell(1,n);
%             %assign behavior value to each point based on data_indiv.sas{istim}.typenames          
%             for k=1:n
%                 frac_beh=(y(k)-range_mean_beh(1))/diff(range_mean_beh);
%                 opts_disp_indiv.set_colors{k}=[(1-frac_beh) frac_beh]*colors_beh;
%                 % [y(index_beh) opts_disp_indiv.set_colors{k}]
%                 istim=ptrs_beh_have(k);
%                 frac_std=(std_beh(istim)-range_std_beh(1))/diff(range_std_beh);
%                 size_std=[(1-frac_std) frac_std]*sizes_beh;
%                 opts_disp_indiv.set_markersizes(k)=round(size_std);
%                 % [k std_beh(istim) frac_std size_std]
%                 % typenames{k}
%             end
%             %
%             aux_disp_indiv=rs_disp_coordsets(data_indiv,setfield(aux,'opts_disp',opts_disp_indiv));
%             set(gcf,'Name',sprintf('dim %1.0f, %s',dim_plot,filename_short));
%             %
%             %plot vectors corresponding to what the regressors project on
%             %
%             coord_groups=aux_disp_indiv.opts_disp.coord_groups;
%             axis_handles=aux_disp_indiv.opts_disp.axis_handles;
%             for isub=1:length(axis_handles)
%                 axes(axis_handles{isub});
%                 hold on;
%                 cg=coord_groups(isub,:);
%                 vec=b{dim_plot}(1+cg); %first entry of b is constant term
%                 coords=data_read.ds{iset}{dim_plot}(:,cg);
%                 coord_mean=mean(coords,1);
%                 coord_scale=sqrt(mean(coords(:).^2)); %scale of coordinates
%                 vec_scale=coord_scale*vec./sqrt(sum(vec.^2)); %normalized
%                 for im=1:2
%                     im_sign=-3+2*im;
%                     if dim_plot>=3
%                         hp=plot3(coord_mean(1)+[0 im_sign*vec_scale(1)],coord_mean(2)+[0 im_sign*vec_scale(2)],coord_mean(3)+[0 im_sign*vec_scale(3)],'k');
%                     else
%                         hp=plot(coord_mean(1)+[0 im_sign*vec_scale(1)],coord_mean(2)+[0 im_sign*vec_scale(2)],'k');
%                     end
%                     set(hp,'LineWidth',line_width);
%                     set(hp,'Color',colors_beh(im,:));
%                 end
%             end
%             %
%             axes('Position',[0.01,0.04,0.01,0.01]); %for text
%             text(0,0,filename_short,'Interpreter','none','FontSize',8);
%             axis off;
%         end %dims_plot
%     end %dims_plot empty?
% end %iset

% mater_ratings_survey_rs: read materials ratings files and analyze them
% 
% creates a mat-file with consolidated 'ratings' data or reads a pre-existing one
%
% shows the representatoinal space for the affordances, based on average ratings, and based on Procrustes consensus of average ratings
% compares consensus with individual ratings
% allows for labelling of selective stimuli
%
% could use strategy of  hlid_vi_coords_behav_rs to use affordance ratings to color the points
% 
%  See also: RS_IMPORT_COORDSETS, RS_KNIT_COORDSETS, RS_DISP_COORDSETS, RS_CONCAT_COORDSETS, RS_EXTRACT_COORDSETS, MATER_RATINGS_STIMSELECT_RS.
%
if ~exist('path_ratings') path_ratings='C:\Users\jdvicto\Documents\JV\ey7977\Materials\FabricAffordanceRatings'; end
if ~exist('string_ratings') string_ratings='_participant_ratings.csv'; end
%
%these are the images that have affordance ratings
images_haveratings=[1,2,4,7,8,9,10,12,13,14,15,21,26,27,30,36,37,38,40,44,49,50,51,55,56,59,65,66,70,73,75,78,80,82,83,84,85,87,88,89,91,92,96,98,99,100,102,103,104,105,109,110,113,116,117,120,121,122,123,129,130];
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
    if getinp('1 to keep all images','d',[0 1],1)
        images_keep=image_nos;
    else
        disp('default is the list of fabrics with 334 x 334 pixel images');
        images_keep=getinp('image numbers to keep','d',[1 max(image_nos)],images_haveratings);
        images_keep=intersect(images_keep,image_nos);
    end
    image_ptrs=find(ismember(image_nos,images_keep));
    %
    %convert ratings data into a structure
    ratings=struct;
    n_imgs=length(image_nos(image_ptrs));
    n_subjs=length(subj_names);
    ratings.n_subjs=n_subjs;
    ratings.n_imgs=n_imgs;
    ratings.n_aff=n_aff;
    ratings.names_aff_all=names_aff_all;
    ratings.names_aff=names_aff;
    ratings.aff_sel=aff_sel;
    ratings.subj_names=subj_names;
    ratings.image_nos=image_nos(image_ptrs);
    ratings.filenames_ratings_all=filenames_ratings_all;
    ratings.filenames_ratings=filenames_ratings;
    ratings.vals=zeros(n_imgs,n_aff,n_subjs);
    ratings.avgs_orig=zeros(n_imgs,n_aff);
    for iaff=1:n_aff
        ratings.vals(:,iaff,:)=reshape(table2array(tables{iaff}(image_ptrs,2:av_col-1)),[n_imgs,1,n_subjs]);
        ratings.avgs_orig(:,iaff)=table2array(tables{iaff}(image_ptrs,'average'));
    end
    ratings.avgs=reshape(mean(ratings.vals,3),[n_imgs,n_aff]);
    disp(sprintf(' comparing averages and recomputed values: max diff is %5.3f',max(abs(ratings.avgs(:)-ratings.avgs_orig(:)))));
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
    % ratings.subj_names=subj_names;type_coords
    % ratings.image_nos=image_nos;
    % ratings.filenames_ratings_all=filenames_ratings_all;
    % ratings.filenames_ratings=filenames_ratings;
    % ratings.vals=zeros(n_imgs,n_aff,n_subjs);
    % ratings.avgs=zeros(n_imgs,n_aff);
if ~exist('nshuffs') nshuffs=100; end
dim_max=min(ratings.n_subjs,ratings.n_aff);
if ~exist('dim_max_consensus') dim_max_consensus=min(6,dim_max); end
dim_max_consensus=getinp('maximum dimension for consensus calculation','d',[1 dim_max],dim_max_consensus);
nshuffs=getinp('number of shuffles for consensus statistics (0 for none)','d',[0 1000],nshuffs);
%create representational space of affordances from each subject and average 
aux_import=struct();
opts_import=auxs_force.opts_import;
%
typenames=cell(ratings.n_imgs,1);
for k=1:ratings.n_imgs
    typenames{k}=cat(2,'img_',zpad(ratings.image_nos(k),3));
end
opts_import.typenames=typenames;
for isubj=1:ratings.n_subjs
    coords=ratings.vals(:,:,isubj);
    [recon_pcaxes,recon_coords,var_ex,var_tot,coord_maxdiff,pca_opts_used]=psg_pcaoffset(coords);
    disp(sprintf(' subj %2.0f (%10s): frac of variance explained',isubj,ratings.subj_names{isubj}));
    disp(diff([0 var_ex])/var_ex(end));
    opts_import.subj_id=ratings.subj_names{isubj};
    opts_import.paradigm_type='affordance ratings';
    opts_import.label=ratings.subj_names{isubj};
    opts_import.label_long=cat(2,'affordance ratings, ',ratings.subj_names{isubj});
    data_ratings=rs_import_coordsets(recon_pcaxes,setfield(aux_import,'opts_import',opts_import));
    data_ratings.sets{1}.pipeline.pca.pca_opts=pca_opts_used; %capture qu, qs, qv, and offset from pca; offset+(coords-offset)*qv is the reconstruction in the PC space. Since qv'*qv is the identity (up to rank), and coords-offset=qu*qs*qv',
    data_ratings.sets{1}.pipeline.pca.var_ex=var_ex; %cumulative variance explained
    if isubj==1
        data_ratings_all=data_ratings;
    else
        data_ratings_all=rs_concat_coordsets(data_ratings_all,data_ratings);
    end   
end
%average subject
coords_avg=mean(ratings.vals,3);
[recon_pcaxes,recon_coords,var_ex,var_tot,coord_maxdiff,pca_opts_used]=psg_pcaoffset(coords_avg);
disp(sprintf(' average across subjects: frac of variance explained'));
disp(diff([0 var_ex])/var_ex(end));
opts_import.subj_id='avg';
opts_import.label='avg';
opts_import.label_long=cat(2,'affordance ratings, ','avg');
data_ratings_combined=rs_import_coordsets(recon_pcaxes,setfield(aux_import,'opts_import',opts_import));
data_ratings_combined.sets{1}.pipeline.pca.pca_opts=pca_opts_used; %capture qu, qs, qv, and offset from pca; offset+(coords-offset)*qv is the reconstruction in the PC space. Since qv'*qv is the identity (up to rank), and coords-offset=qu*qs*qv',
data_ratings_combined.sets{1}.pipeline.pca.var_ex=var_ex; %cumulative variance explained
%
%do a consensus of the individual subjects with options for offset and scaling
%
aux_knit=struct();
aux_knit.opts_check.if_warn=0;
opts_knit=auxs_force.opts_knit;
opts_knit.if_pca=1; %rotate consensus to pca
opts_knit.dim_max_in=dim_max_consensus;
opts_knit.if_stats=double(nshuffs>0);
opts_knit.nshuffs=nshuffs;
ts_pca=cell(2,2); %the transformation to the consensus space, 
data_ratings_components=cell(2,2);
knit_string=cell(2,2);
%
set_labels=cell(1,5);
set_labels{1}='avg'; %first set is average, later sets are individual subjects
iset=0;
for if_offset=0:1
    for if_scale=0:1
        opts_knit.allow_offset=if_offset;
        opts_knit.allow_scale=if_scale;
        knit_string{1+if_offset,1+if_scale}=sprintf('if offset: %1.0f, if scale: %1.0f',if_offset,if_scale);
        iset=iset+1;
        set_labels{iset+1}=knit_string{1+if_offset,1+if_scale}; %first set is 'avg'        
        [data_ratings_knit,aux_knit_out]=rs_knit_coordsets(data_ratings_all,setfield(aux_knit,'opts_knit',opts_knit));
        ts_pca{1+if_offset,1+if_scale}=aux_knit_out.ts_pca;
        for k=1:length(aux_knit_out.components.sets)
            aux_knit_out.components.sets{k}=rmfield(aux_knit_out.components.sets{k},'pipeline');
        end
        data_ratings_components{1+if_offset,1+if_scale}=aux_knit_out.components;
        clear aux_knit_out; %this is large
        if nshuffs>0 %if plotted, annotate with options
            axes('Position',[0.25,0.04,0.01,0.01]); %for text
            text(0,0,knit_string{1+if_offset,1+if_scale},'Interpreter','none','FontSize',8);
            axis off;
        end
        data_ratings_knit.sets{1}.label=knit_string{1+if_offset,1+if_scale};
        data_ratings_knit.sets{1}.label_long=cat(2,' affordances, ',knit_string{1+if_offset,1+if_scale});
        data_ratings_combined=rs_concat_coordsets(data_ratings_combined,data_ratings_knit);
    end
end
%
%plot individual spaces and consensus
%
%compare consensuses created with different options, without alignment
aux_disp=struct();
opts_disp_basic=struct();
opts_disp_basic.data_label_interpreter='none';
opts_disp_basic.callout_amount=0.5;
%
opts_disp=opts_disp_basic;
opts_disp.set_labels=set_labels;
opts_disp.connect_sets_color_mode='split';
opts_disp.set_select=[2:5];
opts_disp.connect_sets_method='list';
opts_disp.connect_sets_list=[2 3;2 4;3 5;4 5;1 5];
aux_disp_out1=rs_disp_coordsets(data_ratings_combined,setfield(aux_disp,'opts_disp',opts_disp));
set(gcf,'Name',cat(2,'consensus options aligned: ',get(gcf,'Name')));
%
if_done=0;
while (if_done==0)
    offset_choice=getinp('offset choice','d',[0 1]);
    scale_choice=getinp('scale choice','d',[0 1]);
    index=2+2*offset_choice+scale_choice;
    %
    data_ratings_use=rs_extract_coordsets(data_ratings_combined,index);% the selected consensus optoin
    %
    %compare average with consensus (optional offset, scale), without alignment
    opts_disp.set_select=[1 index]; 
    opts_disp.connect_sets_method='all';
    aux_disp_out2=rs_disp_coordsets(data_ratings_combined,setfield(aux_disp,'opts_disp',opts_disp));
    set(gcf,'Name',cat(2,'avg vs consensus: ',get(gcf,'Name')));
    %and with alignment
    opts_knit_compare=opts_knit;
    opts_knit_compare.if_offset=0;
    opts_knit_compare.if_scale=0;
    opts_knit_compare.if_stats=0;
    opts_knit_compare.nshuffs=0;
    %with alignment (considering all of the consensus datasets) 
    [data_ratings_compare,aux_compare_out]=rs_knit_coordsets(data_ratings_combined,setfield(aux_knit,'opts_knit',opts_knit_compare));
    aux_disp_out3=rs_disp_coordsets(aux_compare_out.components,setfield(aux_disp,'opts_disp',opts_disp));  
    set(gcf,'Name',cat(2,'avg vs consensus, aligned: ',get(gcf,'Name')));
    %
    %individual subjects and consensus
    opts_disp_subjs=opts_disp_basic;
    opts_disp_subjs.set_labels=['consensus',ratings.subj_names];
    opts_disp_subjs.connect_sets_method='star';
    opts_disp_subjs.connect_sets_color_mode='last';
    %
    rng('default');
    colors_init=rand(ratings.n_subjs,3);
    colors_init=colors_init./repmat(sqrt(sum(colors_init.^2,2)),1,3);
    colors=cat(1,zeros(1,3),0.5+colors_init*0.5);
    opts_disp_subjs.set_colors=cell(1,ratings.n_subjs+1);
    for k=1:ratings.n_subjs+1
        opts_disp_subjs.set_colors{1,k}=colors(k,:);
    end
    aux_disp_out_subjs=rs_disp_coordsets(rs_concat_coordsets(data_ratings_use,data_ratings_components{1+offset_choice,1+scale_choice}),...
        setfield(aux_disp,'opts_disp',opts_disp_subjs));
    set(gcf,'Name',cat(2,'indiv subjects and consensus: ',get(gcf,'Name')));
    axes('Position',[0.01,0.04,0.01,0.01]); %for text
    text(0,0,knit_string{1+offset_choice,1+scale_choice},'Interpreter','none','FontSize',8);
    axis off;
    if_done=getinp('1 if done','d',[0 1]);
end
%
if_done=0;
image_nos_labelled=ratings.image_nos(:)';
while (if_done==0)
    %plot with alignment, but selective labeling, dims 2 and 3
    if getinp('replot with selective labeling','d',[0 1])
        opts_disp_sel=opts_disp;
        image_nos_labelled=getinp('image numbers to label','d',[min(ratings.image_nos) max(ratings.image_nos)],image_nos_labelled);
        image_nos_labelled=intersect(image_nos_labelled,ratings.image_nos);
        image_nos_labelled=image_nos_labelled(:)';
        image_nos_ptrs=zeros(1,length(image_nos_labelled));
        for k=1:length(image_nos_labelled)
            image_nos_ptrs(k)=find(ratings.image_nos==image_nos_labelled(k));
        end
        opts_disp_sel.data_label_method='list';
        opts_disp_sel.data_label_list=image_nos_ptrs;
        opts_disp_sel.callout_amount=getinp('callout length','f',[0 1],opts_disp_sel.callout_amount);
        %
        lab_string=sprintf('%3.0f points, %3.0f labelled',length(ratings.image_nos),length(image_nos_labelled));
        for dim_select=[2:3]
            opts_disp_sel.dim_select=dim_select;
            aux_disp_outsel=rs_disp_coordsets(aux_compare_out.components,setfield(aux_disp,'opts_disp',opts_disp_sel));
            set(gcf,'Name',cat(2,lab_string,', ',get(gcf,'Name')));
            axes('Position',[0.01,0.04,0.01,0.01]); %for text
            text(0,0,lab_string,'Interpreter','none','FontSize',8);
            axis off;
        end %dim_select
    else
        if_done=1;
    end
end

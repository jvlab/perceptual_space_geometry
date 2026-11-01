% mater_ratings_survey_rs: select stimuli based on materials ratings files
%
% streamlined version of mater_ratings_survey_rs, to select stimulus subsets based on similar coverage of the affordance space based on the available images
%  * individual subject data are only combined (Procrustes consensus) with offset and scaling allowed
%  * this is the only set these are aligned with the "average" subject's ratings from Grace
%  * random partitions of the images that are available as 334 x 334 are made into set A, set B, and their overlap
%  * summary stats of the above and sifting based on quantiles are done as an aid to selection
% 
%  See also: RS_IMPORT_COORDSETS, RS_KNIT_COORDSETS, RS_DISP_COORDSETS, RS_CONCAT_COORDSETS, RS_EXTRACT_COORDSETS, MATER_RATINGS_SURVEY_RS.
%
%section to force btc defaults, even if rs_aux_defaults.mat has been created or modified
if ~exist('aux_force_filename') aux_force_filename='rs_aux_defaults_btc.mat'; end
auxs_force=struct;
opts_needed={'opts_read','opts_rays','opts_check','opts_align','opts_import','opts_qpred','opts_knit','opts_write'};
for k=1:length(opts_needed)
    auxs_force.(opts_needed{k})=rs_aux_force(opts_needed{k},[],aux_force_filename);
end
%
if ~exist('path_ratings') path_ratings='C:\Users\jdvicto\Documents\JV\ey7977\Materials\FabricAffordanceRatings'; end
if ~exist('string_ratings') string_ratings='_participant_ratings.csv'; end
%
%affordance_ratings_all_10Oct26.mat:       all 80 fabrics, all affordances
%affordance_ratings_334_13Oct26.mat:       the 61 fabrics large enough to have 334 x 334 pixel images, all affordances
%affordance_ratings_334_6aff_14Oct26.mat:  the 61 fabrics large enough to have 334 x 334 pixel images, the 6 affordances used by Grace to construct the affordance space
%
if ~exist('ratings_fiile') ratings_file='affordance_ratings_334_13Oct26.mat'; end %61 images with ratings of affordances, pca based on all 10 affordances
load (ratings_file);
disp(sprintf(' loaded %s',ratings_file));
%
%these are the images that have affordance ratings
images_haveratings=[1,2,4,7,8,9,10,12,13,14,15,21,26,27,30,36,37,38,40,44,49,50,51,55,56,59,65,66,70,73,75,78,80,82,83,84,85,87,88,89,91,92,96,98,99,100,102,103,104,105,109,110,113,116,117,120,121,122,123,129,130];
%these are the images that have 334 x 334 patches, as of 10Oct26
images_have334=[1,2,4,15,26,27,30,36,37,38,40,44,49,50,51,55,56,59,65,66,70,73,75,78,80,82,83,84,85,87,88,91,92,96,99,102,103,104,109,120,123,129,130];
%
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
if ~exist('nshuffs') nshuffs=0; end
dim_max=min(ratings.n_subjs,ratings.n_aff);
if ~exist('dim_max_consensus') dim_max_consensus=min(3,dim_max); end
dim_max_consensus=getinp('maximum dimension for consensus calculation','d',[1 dim_max],dim_max_consensus);
nshuffs=getinp('number of shuffles for consensus statistics (0 for none)','d',[0 1000],nshuffs);
dims_show=getinp('dimensions to show','d',[2 3],2);
n_design=getinp('number of stimuli in each session','d',[25 37],37);
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
opts_knit.if_stats=1;
opts_knit.if_plot=1;
%
set_labels=cell(1,2);
set_labels{1}='avg'; %first set is average, later sets are individual subjects
iset=0;
if_offset=1;
if_scale=1;
opts_knit.allow_offset=if_offset;
opts_knit.allow_scale=if_scale;
opts_knit.if_log=0;
knit_string=sprintf('if offset: %1.0f, if scale: %1.0f',if_offset,if_scale);
iset=iset+1;
set_labels{iset+1}=knit_string; %first set is 'avg'        
data_ratings_knit=rs_knit_coordsets(data_ratings_all,setfield(aux_knit,'opts_knit',opts_knit));
set(gcf,'Name','consensus across subjects');
axes('Position',[0.25,0.04,0.01,0.01]); %for text
text(0,0,cat(2,'alignment across subjects with ',knit_string),'Interpreter','none','FontSize',8);
axis off;
%
data_ratings_knit.sets{1}.label=knit_string;
data_ratings_knit.sets{1}.label_long=cat(2,' affordances, ',knit_string);
data_ratings_combined=rs_concat_coordsets(data_ratings_combined,data_ratings_knit); %set 1 is avg, set 2 is consensus
%
%align the average and consensus spaces, with offset and scaling
%
[data_aligned,aux_aligned_out]=rs_knit_coordsets(data_ratings_combined,setfield(aux_knit,'opts_knit',setfields(opts_knit,{'allow_offset','allow_scale'},{1,1})));
set(gcf,'Name','aligning average and consensus');
axes('Position',[0.25,0.04,0.01,0.01]); %for text
text(0,0,'alignment of average and consensus','Interpreter','none','FontSize',8);
axis off;
data_ratings_aligned=aux_aligned_out.components; %set 1 is avg, set 2 is consensus, now Procrustes aligned with each other
%
% plot these, with labels for the stimuli for which we have images
%
aux_disp=struct();
opts_disp_basic=struct();
opts_disp_basic.data_label_interpreter='none';
opts_disp_basic.callout_amount=0.1;
%
opts_disp=opts_disp_basic;
opts_disp.connect_sets_color_mode='split';
opts_disp.connect_sets_method='all';
opts_disp.set_labels={'avg','consensus'};
opts_disp.set_colors={'k','r'};
%
image_nos_all=ratings.image_nos(:)';
image_nos_labelled=images_have334;
n_have=length(images_have334);
for k=1:n_have
    image_nos_ptrs(k)=find(ratings.image_nos==image_nos_labelled(k));
end
lab_string=sprintf('%3.0f points, %3.0f labelled',length(ratings.image_nos),length(image_nos_labelled));
opts_disp.data_label_method='list';
opts_disp.data_label_list=image_nos_ptrs;
%
for idim=dims_show
    opts_disp.dim_select=idim;
    aux_disp_aligned=rs_disp_coordsets(data_ratings_aligned,setfield(aux_disp,'opts_disp',opts_disp));
    set(gcf,'Name',cat(2,'aligning average and consensus, ',get(gcf,'Name')));
    axes('Position',[0.01,0.04,0.01,0.01]); %for text
    text(0,0,lab_string,'Interpreter','none','FontSize',8);
    axis off;
end
n_subsets=4;
%subdivide into don't have, have in set A, have in set B, have in both
image_ptrs=cell(1,n_subsets);
n_one=n_have-n_design;
n_both=n_design-n_one;
design_string=sprintf('%2.0f stimuli per session, containing %2.0f stimuli only in one session, %2.0f stimuli in both',n_design,n_one,n_both);
disp(design_string);
%
tags_sorted=[repmat(2,n_one,1);repmat(3,n_one,1);repmat(4,n_both,1)]'; %2 is set A only, 3 is set B only, 4 is both
%
%create a data structure with four copies of data_ratings_aligned, and labels
data_ratings_allsubsets=data_ratings_aligned;
for k=2:n_subsets
    data_ratings_allsubsets=rs_concat_coordsets(data_ratings_allsubsets,data_ratings_aligned);
end
set_labels_subsets=cell(1,2*n_subsets);
set_labels_base={'unused','A','B','both'};
for isub=1:n_subsets
    set_labels_subsets{2*isub-1}=cat(2,set_labels_base{isub},', avg');
    set_labels_subsets{2*isub}=cat(2,set_labels_base{isub},', consensus');
end %subset
%
n_random=1;
while (n_random>0)
    n_random=getinp('1 or more to make random subsets','d',[0 10000]);
    if n_random>0
        if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
        if (if_frozen~=0) 
            rng('default');
            if (if_frozen<0)
                rand(1,abs(if_frozen));
            end
        else
            rng('shuffle');
        end
        dstats=getinp('dimension for statistics','d',[1 dim_max_consensus],2);
        %calculate affordance coordinates for random subsets 
        rperm=zeros(n_random,n_have);
        aff_coords=zeros(n_design,dstats,2,n_random,2); %d1: each stim d2: the coords, d3: set A or set B, d4: which random set, d5: avg or consensus 
        image_iab=cell(n_random,2);
        for r=1:n_random
            rperm(r,:)=randperm(n_have);
            %
            tags=tags_sorted(rperm(r,:));
            for i_ab=1:2
                image_use=images_have334(find(tags==i_ab+1)); %tag=2 or 3, depending on set A or set B
                image_use=union(image_use,images_have334(find(tags==4))); %and the overlap set
                image_iab{i_ab}=zeros(1,length(image_use));
                for k=1:length(image_use)
                    image_iab{r,i_ab}(1,k)=find(ratings.image_nos==image_use(k));
                end
                for i_avcon=1:2 %avg or consensus
                    aff_coords(:,:,i_ab,r,i_avcon)=data_ratings_allsubsets.ds{i_avcon}{dstats}(image_iab{r,i_ab},:);
                end
            end
        end
        %calculate statistics
        %means, vector distance between groups
        aff_means_vecdiff=reshape(sqrt(sum(diff(mean(aff_coords,1),1,3).^2,2)),[n_random 2]); % d1: which random set, d2: avg or consensus
        %covariances, compared by the natural distance of for covariance matrices (sum of log-squareds of eigenvalues of inv(a)*b
        cov_dsq=zeros(n_random,2);
        for r=1:n_random
            for i_avcon=1:2
                cov_a=cov(aff_coords(:,:,1,r,i_avcon));
                cov_b=cov(aff_coords(:,:,2,r,i_avcon));
                cov_mtx=inv(cov_a)*cov_b;
                cov_dsq(r,i_avcon)=sum(log(eig(cov_mtx)).^2);
            end
        end
        %sift configurations based on quantiles of selected statistics
        sift_stats_all=[aff_means_vecdiff,cov_dsq];
        disp('stats are differences in means (avg, consensus) and differnces in covariances (avg, consensus)');
        sift_stats_sel=[1:size(sift_stats_all,2)];
        while (sift_stats_sel>0)
            sift_stats_sel=getinp('selection (0 to end)','d',[0 size(sift_stats_all,2)],sift_stats_sel);
            if all(sift_stats_sel>0)
                sift_stats=sift_stats_all(:,sift_stats_sel);
                quantile_sift=1.0;
                sifted=[1:n_random];
                while ~isempty(sifted) & quantile_sift>1/n_random
                    sifted_last=sifted;
                    quantile_sift=quantile_sift*0.5;
                    quantile_vals=quantile(sift_stats,quantile_sift,1);
                    sifted=find(all(sift_stats<=repmat(quantile_vals,[n_random,1]),2));
                    sifted_eachcount=sum(double(sift_stats<=repmat(quantile_vals,[n_random,1])),1);
                    disp(sprintf('quantile %7.4f            cutoffs: %s',quantile_sift,sprintf(' %8.5f',quantile_vals)));
                    disp(sprintf(' configurations remaining in each: %s',sprintf(' %8.0f',sifted_eachcount)));
                    disp(sprintf(' configurations remaining in all: %4.0f',length(sifted)));
                    disp(sprintf(' %6.0f',sifted(1:min(10,end))'));
                end
            end
        end
        %show affordance coordinates for selected partitions
        if_show=1;
        while (if_show>0)
            if_show=getinp('random subset to show (0 to end)','d',[0 n_random]);
            if if_show>0
                data_ratings_subsets=data_ratings_allsubsets;
                r=if_show;
                image_ptrs{1}=setdiff([1:ratings.n_imgs],union(image_iab{r,1},image_iab{r,2})); %images not in either set
                image_ptrs{2}=setdiff(image_iab{r,1},image_iab{r,2}); %images in A but not B
                image_ptrs{3}=setdiff(image_iab{r,2},image_iab{r,1}); %images in B but not A
                image_ptrs{4}=intersect(image_iab{r,1},image_iab{r,2});
                for isub=1:n_subsets
                    for i_avcon=1:2
                        exclude=setdiff([1:ratings.n_imgs],image_ptrs{isub});
                        i_avcon=2*(isub-1)+i_avcon;
                        for idim=1:dim_max_consensus
                            data_ratings_subsets.ds{i_avcon}{idim}(exclude,:)=NaN;
                        end %idim
                    end %avg or consensus
                end %subset
                opts_disp_subsets=opts_disp_basic;
                opts_disp_subsets.connect_sets_color_mode='first';
                opts_disp_subsets.connect_sets_method='list';
                opts_disp_subsets.connect_sets_list=reshape(1:2*n_subsets,2,n_subsets)';
                opts_disp_subsets.set_labels=set_labels_subsets;
                opts_disp_subsets.set_colors={'k','k','r','r',[0.3 0.3 1.0],[0.3 0.3 1.0],'m','m'};
                opts_disp_subsets.set_markers={'o','x'};
                opts_disp_subsets.callout_colors='set_colors';
                opts_disp_subsets.data_label_setsel_method='list';
                opts_disp_subsets.data_label_setsel_list=2*[1:n_subsets];
                %
                for idim=dims_show
                    opts_disp_subsets.dim_select=idim;
                    aux_disp_subset_out=rs_disp_coordsets(data_ratings_subsets,setfield(aux_disp,'opts_disp',opts_disp_subsets));
                    set(gcf,'Name',cat(2,'subsets, ',get(gcf,'Name')));
                    axes('Position',[0.01,0.04,0.01,0.01]); %for text
                    text(0,0,cat(2,design_string,sprintf(', choice %4.0f of %4.0f',if_show,n_random)),'Interpreter','none','FontSize',8);
                    axis off;
                end
            end
        end %if_show
    end
end %n_random>0
%need to sift based on quantiles of centroid shifts and variance differences

%hlid_orn_merge2: read multiple sets of raw orn data files with non-overlapping stimuli
% handle missing data, merge the orn reponses, and write a coordinate file
%
% Missing data within a set are filled in as in hlid_orn_merge,
% but then the sets (with non-overlapping stimuli) are appended, using corresponding ORNs
%
% the rois.glomeurli fields must match across sets
%
% Built on hlid_orn_merge, which only handles a single set of stimuli; missing data are 
% filled in based on keeping track of ORN's, and finding a multiplicate and additive constant for 
% each dataset
%s
%   See also:  HLID_SETUP, HLID_RASTM2COORDS_DEMO, HLID_RASTIM2COORDS_POOL, AFALWT, HLID_SVD_COORDS, HLID_PLOT_COORDS,
%   HLID_DA_STIMSELECT, HLID_ORN_MERGE.
%
hlid_setup;
if ~exist('opts_dasel')
    opts_dasel=struct;
end
if ~exist('hist_bins')  hist_bins=50; end
if ~exist('hist_quantiles') hist_quantiles=[.05 .25 .5 .75 .95]; end
nsets=getinp('number of sets of stimuli','d',[1 100]);
if_restore_size=getinp('1 to restore responses to match overall size of original data (0: legacy)','d',[0 1]);
if_submean=getinp('1 to subtract mean across responses before creating coordinates','d',[0 1],0);
%
filenames_short=cell(1,nsets);
pathname_short=cell(1,nsets);
nfiles=zeros(1,nsets);
nflies_use=zeros(1,nsets);
s=cell(1,nsets);
stimulus_names_set=cell(1,nsets);
stim_labels_set=cell(1,nsets);
nstims=zeros(1,nsets);
glomeruli_use=cell(1,nsets); %pointers to glomeruli used in each set
nglomeruli_use=zeros(1,nsets); %number of glomeruli used in each set
resps_raw=cell(1,nsets); %responses within each dataset prior to fill-in merge
resps=cell(1,nsets); %responses after fill-in merge within each set
files_use=cell(1,nsets); %files used within each set
dsid=cell(1,nsets);
%
for iset=1:nsets
    disp(sprintf('Enter set %1.0f',iset))
    [filenames_short{iset},pathname{iset}]=uigetfile('*fly*.mat',sprintf('Select raw ORN data files for set %1.0f',iset),'Multiselect','on');
    if ~iscell(filenames_short{iset})
        filenames_short{iset}=cellstr(filenames_short{iset});
    end
    nfiles(iset)=length(filenames_short{iset});
    files_use{iset}=[];
    nancols=cell(nfiles(iset),1);
    nanrows=cell(nfiles(iset),1);
    dsid{iset}=[];
    %
    %verify consistency of names and numbers of glomeruli and stimuli
    %
    for ifile=1:nfiles(iset)
        s{iset}{ifile}=load(cat(2,pathname{iset},filenames_short{iset}{ifile}));
        [s{iset}{ifile},optsused_dasel]=hlid_da_stimselect(s{iset}{ifile},opts_dasel);
        dsid_this=s{iset}{ifile}.meta.title;
        %'-'can be used within fields of file name
        dsid_this=strrep(dsid_this,'/','-');
        dsid_this=strrep(dsid_this,'\','-');
        dsid_this=strrep(dsid_this,'_','-');
        while contains(dsid_this,'--')
            dsid_this=strrep(dsid_this,'--','-');
        end
        %
        dsid{iset}=strvcat(dsid{iset},dsid_this);
        if (ifile==1)
            if (iset==1)
                glomeruli=s{iset}{ifile}.rois.glomeruli;
                nglomeruli=length(glomeruli);
            end
            stimulus_names_set{iset}=s{iset}{ifile}.response_amplitude_stim.stim';
            nstims(iset)=length(stimulus_names_set{iset});
        end
        glomeruli_check=s{iset}{ifile}.rois.glomeruli;
        stimulus_names_check=s{iset}{ifile}.response_amplitude_stim.stim';
        resps_read=s{iset}{ifile}.response_amplitude_stim.mean_peak;
        nancols{ifile}=find(all(isnan(resps_read),1));
        nanrows{ifile}=find(all(isnan(resps_read),2));
        %
        disp(sprintf('file %2.0f (%20s) read, %3.0f glomeruli (%3.0f all NaN), %3.0f stimuli (%3.0f all NaN)',ifile,filenames_short{iset}{ifile},...
            length(glomeruli_check),length(nancols{ifile}),...
            length(stimulus_names_check),length(nanrows{ifile})));
        ifok=1;
        if length(glomeruli_check)~=nglomeruli
            ifok=0;
            disp('number of glomeruli does not match')
        elseif any(strcmp(glomeruli,glomeruli_check)==0)
            ifok=0;
            disp('names of glomeruli do not match')
        end
        if length(stimulus_names_check)~=nstims(iset)
            ifok=0;
            disp('number of stimuli does not match')
        elseif any(strcmp(stimulus_names_set{iset},stimulus_names_check)==0)
            ifok=0;
            disp('names of stimuli do not match')
        end
        if (ifok==1)
            files_use{iset}=[files_use{iset},ifile];
        end
    end
    %shorten stimulus names
    stim_labels_set{iset}=stimulus_names_set{iset};
    for istim=1:nstims(iset)
        if contains(stim_labels_set{iset}{istim},'@')
            stim_labels_set{iset}{istim}=stim_labels_set{iset}{istim}(1:min(find(stim_labels_set{iset}{istim}=='@')-1));
        end
        stim_labels_set{iset}{istim}=deblank(stim_labels_set{iset}{istim});
    end
    %
    % select datasets
    %
    files_use{iset}=getinp('list of files to use','d',[1 nfiles(iset)],files_use{iset});
    nfiles_use(iset)=length(files_use{iset});
    resps_raw{iset}=zeros(nstims(iset),nglomeruli,nfiles_use(iset));
    glomeruli_missing=zeros(nglomeruli,nfiles_use(iset));
    for ifile_ptr=1:nfiles_use(iset)
        ifile=files_use{iset}(ifile_ptr);
        glomeruli_missing(nancols{ifile},ifile_ptr)=1;
        resps_raw{iset}(:,:,ifile_ptr)=s{iset}{ifile}.response_amplitude_stim.mean_peak;
    end
    for ipresent=0:nfiles_use(iset)
        disp(sprintf('number of glomeruli present in at least %2.0f datasets: %3.0f (missing in %2.0f or less)',...
            ipresent,sum(sum(glomeruli_missing,2)<=(nfiles_use(iset)-ipresent)),nfiles_use(iset)-ipresent));
    end
    min_present=getinp('minimum number of preps that a glomerulus must be present in','d',[0 nfiles_use(iset)],1);
    glomeruli_use{iset}=find(sum(glomeruli_missing,2)<=nfiles_use(iset)-min_present);
    nglomeruli_use(iset)=length(glomeruli_use{iset});
    resps_gu=resps_raw{iset}(:,glomeruli_use{iset},:);
    %
    %fill in missing data by affine interpolation
    % 
    resps_gur=reshape(resps_gu,[nstims(iset)*nglomeruli_use(iset),nfiles_use(iset)]);
    if ~exist('afalwt_opts') afalwt_opts=struct;end
    resps_tofill=isnan(resps_gur);
    if_canfill=1;
    if any(all(resps_tofill==1,1))
        disp('cannot fill in missing data, no stimuli present for some (glomerulus,file) pair')
        if_canfill=0;
    end
    if any(all(resps_tofill==1,2))
        disp('cannot fill in missing data, no (glomerulus,file) pair present for some stimulus')
        if_canfill=0;
    end
    if if_canfill
        [afalwt_fit,afalwt_b_change,afalwt_optsused]=afalwt(resps_gur,1-resps_tofill,afalwt_opts);
        resps_gur_fitted=(afalwt_fit.x_true*afalwt_fit.b_norm+repmat(afalwt_fit.a,size(resps_gur,1),1)); %interpolated data
        resps_gur_filled=resps_gur;
        resps_gur_filled(resps_tofill)=resps_gur_fitted(resps_tofill);
        resps_gu_filled=reshape(resps_gur_filled,[nstims(iset) nglomeruli_use(iset) nfiles_use(iset)]);
    else
        resps_gu_filled=resps_gu;
    end
    %
    disp(sprintf('%4.0f NaN values filled in',sum(resps_tofill(:))));
    %
    %show data, before glomerulus selection and after selection, normalized within each dataset
    %
    for ifig=1:3
        switch ifig
            case 1
                figname='all raw data';
                resps_plot=resps_raw{iset};
                glomeruli_plot=[1:nglomeruli];
            case 2
                figname='raw data from selected glomeruli';
                resps_plot=resps_gu;
                glomeruli_plot=glomeruli_use{iset};
            case 3
                figname='raw data from selected glomeruli with missing data filled in';
                resps_plot=resps_gu_filled;
                glomeruli_plot=glomeruli_use{iset};
        end
        figname=cat(2,sprintf('set %2.0f: ',iset),figname);
        figure;
        set(gcf,'Position',[50 100 1800 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',figname);
        [nr,nc]=nicesubp(nfiles_use(iset));
        for ifile_ptr=1:nfiles_use(iset)
            ifile=files_use{iset}(ifile_ptr);
            subplot(nr,nc,ifile_ptr);
            imagesc(resps_plot(:,:,ifile_ptr));
            minmax=[min(min(resps_plot(:,:,ifile_ptr),[],'omitnan')),max(max(resps_plot(:,:,ifile_ptr),[],'omitnan'))];
            title_string=sprintf('set %1.0f file %2.0f: %s  [%6.3f %6.3f]',iset,ifile,strrep(filenames_short{iset}{ifile},'.mat',''),minmax);
            title(title_string,'Interpreter','none');
            set(gca,'FontSize',7);
            set(gca,'XTick',[1:length(glomeruli_plot)]);
            set(gca,'XTickLabel',glomeruli(glomeruli_plot));
            set(gca,'YTick',[1:nstims(iset)]);
            set(gca,'YTickLabel',stim_labels_set{iset});
        end
        axes('Position',[0.01,0.02,0.01,0.01]); %for text
        text(0,0,figname,'Interpreter','none');
        axis off;
    end %ifig
    %
    if any(isnan(resps_gu_filled(:)))
        disp('Cannot proceed. Not all NaNs have been filled in.')
    end
    stim_labels_set{iset}=stimulus_names_set{iset};
    resps_set{iset}=reshape(afalwt_fit.x_true,[nstims(iset) nglomeruli_use(iset)]); %use regression slope as response measure
    if if_restore_size %restorr size if needed
        resps_set{iset}=resps_set{iset}*geomean(afalwt_fit.b_norm);
    end
    %
    %condition the data and stimulus names
    %
    if (if_submean)
        resps_set{iset}=resps_set{iset}-repmat(mean(resps_set{iset},1),nstims(iset),1);
    end
    %
    for istim=1:nstims(iset)
        if contains(stim_labels_set{iset}{istim},'@')
            stim_labels_set{iset}{istim}=stim_labels_set{iset}{istim}(1:min(find(stim_labels_set{iset}{istim}=='@')-1));
        end
        stim_labels_set{iset}{istim}=deblank(stim_labels_set{iset}{istim});
    end
end
%
%calculations common to all datasets
%
resp_range=[Inf -Inf];
for iset=1:nsets
    disp(sprintf('set %2.0f has %3.0f files used (of %3.0f), %3.0f stimuli, and %3.0f glomeruli used',iset,nfiles_use(iset),nfiles(iset),nstims(iset),nglomeruli_use(iset)))
    resp_range(1)=min(resp_range(1),min(resps_set{iset}(:)));
    resp_range(2)=max(resp_range(2),max(resps_set{iset}(:)));
end
dsid_show=sprintf('affine merging %2.0f files: %s,...,%s',sum(nfiles),deblank(dsid{1}(files_use{1}(1),:)),deblank(dsid(files_use{nsets}(nfiles_use(nsets)),:)));
%
%overall response histogram and quantiles
%
figure;
set(gcf,'Position',[50 100 800 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name','histogram');
for iset=1:nsets
    disp(sprintf('quantiles for set %1.0f',iset))
    subplot(2,nsets,iset);
    hist(resps_set{iset}(:),hist_bins);
    xlabel('response')
    ylabel('counts');
    set(gca,'XLim',resp_range);
    title(sprintf('set %1.0f',iset))
    quantiles=quantile(resps_set{iset}(:),hist_quantiles);
    for k=1:length(hist_quantiles)
        axes('Position',[0.01+(iset-1)/nsets,0.06+0.04*k,0.01,0.01]); %for text
        qt=sprintf(' quantile %5.3f: %7.3f',hist_quantiles(k),quantiles(k));
        disp(qt)
        text(0,0,qt,'Interpreter','none');
        axis off;
    end
end
axes('Position',[0.01,0.02,0.01,0.01]);
text(0,0,sprintf('restore size=%2.0f, subtract mean=%2.0f',if_restore_size,if_submean));
axis off;
%
%quantities calculated from all sets
%
maxdim_allowed=Inf;
stimulus_names=cell(0);
stim_labels=cell(0);
resp_range=[Inf -Inf];
for iset=1:nsets
    maxdim_allowed=min(maxdim_allowed,min(size(resps_set{iset}))-if_submean);
    stimulus_names=[stimulus_names;stimulus_names_set{iset}];
    stim_labels=[stim_labels;stim_labels_set{iset}];
    resp_range=[min(resp_range(1),min(resps_set{iset}(:))) max(resp_range(2),max(resps_set{iset}(:)))];
end
if length(unique(stim_labels))~=sum(nstims)
    warning('Not all stimulus names are unique');
end
if length(unique(stim_labels))~=sum(nstims)
    warning('Not all stim labels are unique');
end
maxdim=getinp('maximum number of dimensions for coordinate file','d',[1 maxdim_allowed],maxdim_allowed);
maxdim_use=maxdim_allowed;
%
%initialize the quantities to save
%
f_base=struct;
for iset=1:nsets
    f_base.metadata{iset}=s{iset}{files_use{iset}(1)}.meta; %original metadata from Hong Lab
    f_base.dsid{iset}=dsid{iset}(files_use{iset},:); %data set ID, with special chars turned into -
end
f_base.stimulus_names=strvcat(stimulus_names); %original stimulus names
f_base.stim_labels=strvcat(stim_labels); %shortened names for plotting
f_base.resps=resps_set; %original responses
if if_restore_size==0
    f_base.coord_opts.resp_type='response_amplitude_stim, after affine filling in'; %original field for responses from Hong Lab
else
    f_base.coord_opts.resp_type='response_amplitude_stim, after affine filling in and restoring to original size';
end
%
%combined set either has intersection of glomeruli in components, or glomeruli in union, with missing responses set to 0
%
ncombm=2; %number of combination modes: 1: intersection, 2: union, set missing to 0
resps_combs=cell(1,ncombm);
f_combm=cell(1,ncombm);
s_diag_all_combm=cell(1,ncombm);
u_full_combm=cell(1,ncombm);
v_full_combm=cell(1,ncombm);
s_full_combm=cell(1,ncombm);
coords_all_combm=cell(1,ncombm);
%
for icombm=1:2
    stims_sofar=0;
    switch icombm
        case 1             
            comb_label='intersection';
            glomeruli_combined=[1:nglomeruli];
            for iset=1:nsets
                glomeruli_combined=intersect(glomeruli_combined,glomeruli_use{iset});
            end
            resps_combined{icombm}=zeros(sum(nstims),length(glomeruli_combined));
            for iset=1:nsets
                resps_combined{icombm}(stims_sofar+[1:nstims(iset)],:)=resps_set{iset}(:,find(ismember(glomeruli_use{iset},glomeruli_combined)));
                stims_sofar=stims_sofar+nstims(iset);
            end
        case 2
            comb_label='union, missing set to 0';
            glomeruli_combined=[];
            for iset=1:nsets
                glomeruli_combined=union(glomeruli_combined,glomeruli_use{iset});
            end
            resps_combined{icombm}=zeros(sum(nstims),length(glomeruli_combined));
            for iset=1:nsets
                resps_combined{icombm}(stims_sofar+[1:nstims(iset)],ismember(glomeruli_combined,glomeruli_use{iset}))=resps_set{iset};
                stims_sofar=stims_sofar+nstims(iset);
            end
    end
    disp(sprintf(' combination method %20s: responses from %2.0f glomeruli',comb_label,length(glomeruli_combined)))
    %
    %create coords by SVD and add metadata and plot
    %
    [f,s_diag_all,u_full,v_full,s_full,coords_all]=hlid_coords_svd(f_base,resps_combined{icombm},maxdim,maxdim_use,if_submean);
    %save combined variables
    f_combm{icombm}=f;
    s_diag_all_combm{icombm}=s_diag_all;
    u_full_combm{icombm}=u_full;
    v_full_combm{icombm}=v_full;
    s_full_combm{icombm}=s_full;
    coords_all__combm{icombm}=coords_all;
    %
    if getinp(sprintf('1 if ok to write a coordinate file for %s',comb_label),'d',[0 1])
        data_fullname_write_def=strrep(hlid_opts.coord_data_fullname_write_def,'dsid','merged');
        data_fullname_write_def=strrep(data_fullname_write_def,'kc_soma_nls','orn_merged');
        data_fullname_write=getinp('coordinate file name to write','s',[],data_fullname_write_def);
        save(data_fullname_write,'-struct','f');
        disp(sprintf('wrote %s',data_fullname_write));
    end
    resps=resps_combined{icombm};
    roi_names=glomeruli(glomeruli_combined);
    hlid_coords_plot;
    %
    %special plots for combining datasets
    %
    figname=sprintf('combined sets via %s',comb_label);
    figure;
    set(gcf,'Position',[50 100 1800 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',figname);
    subplot(1,2,1)
    imagesc(resps_combined{icombm},resp_range);
    hold on;
    for iset=1:nsets-1
        stimct=sum(nstims(1:iset));
        plot([0 length(glomeruli_combined)]+[0.5 0.5],stimct+[0.5 0.5],'k','LineWidth',2);
    end
    title_string=figname;
    title(title_string,'Interpreter','none');
    set(gca,'FontSize',7);
    set(gca,'XTick',[1:length(glomeruli_combined)]);
    set(gca,'XTickLabel',glomeruli(glomeruli_combined));
    set(gca,'YTick',[1:sum(nstims)]);
    set(gca,'YTickLabel',stim_labels);
    colorbar;
    % axes('Position',[0.01,0.02,0.01,0.01]); %for text
    % text(0,0,figname,'Interpreter','none');
    % axis off;
    %
end %icombm

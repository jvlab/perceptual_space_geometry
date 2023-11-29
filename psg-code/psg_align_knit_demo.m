%psg_align_knit_demo: demonstration of alignment and knitting together of multiple datasets
% that have partially overlapping stimuli
%
% all datasets must have dimension lists beginning at 1 and without gaps
% aligned datasets and metadata (ds_align,sas_align) will have a NaN where there is no match
%
%  See also: PSG_ALIGN_COORDSETS, PSG_COORD_PIPE_PROC, PSG_GET_COORDSETS, PSG_READ_COORDDATA,
%    PROCRUSTES_CONSENSUS, PROCRUSTES_CONSENSUS_PTL_TEST, PSG_FINDRAYS.
%

%main structures and workflow:
%ds{nsets},            sas{nsets}: original datasets and metadata
%ds_align{nsets},      sas_align{nsets}: datasets with NaN's inserted to align the stimuli
%d_knitted,            sa_pooled: consensus rotation of ds_align, all stimuli, and metadata
%ds_components{nsets}, sas_align{nses}: components of d_knitted, corrsponding to original datasets, but with NaNs -- these are Procrustes transforms of ds_align
%ds_nonan{nsets}       sas_nonan{nsets}: components stripped of NaNs
%
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_align') opts_align=struct(); end %for psg_align_coordsets
if ~exist('opts_nonan') opts_nonan=struct(); end %for psg_remnan_coordsets
if ~exist('opts_pcon') opts_pcon=struct(); end % for procrustes_consensus
if ~exist('pcon_dim_max') pcon_dim_max=3; end %dimensions for alignment
%
if ~exist('color_list') color_list='rmbcg'; end
%
disp('This will attempt to knit together two or more coordinate datasets.');
%
opts_read=filldefault(opts_read,'input_type',0); %either experimental data or model
opts_align=filldefault(opts_align,'if_log',1);
opts_nonan=filldefault(opts_nonan,'if_log',1);
%
opts_pcon=filldefault(opts_pcon,'allow_reflection',1);
opts_pcon=filldefault(opts_pcon,'allow_offset',1);
opts_pcon=filldefault(opts_pcon,'allow_scale',0);
opts_pcon=filldefault(opts_pcon,'max_niters',1000); %nonstandard max
%
nsets=getinp('number of datasets','d',[1 100]);
[sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,[],nsets); %get the datasets
[sets_align,ds_align,sas_align,ovlp_array,sa_pooled,opts_align_used]=psg_align_coordsets(sets,ds,sas,opts_align); %align the stimulus names
nstims_all=sets_align{1}.nstims;
disp(sprintf('total stimuli: %3.0f',nstims_all));
%
for iset=1:nsets
    if (iset==1)
        dim_list_all=sets{iset}.dim_list;
    else
        dim_list_all=intersect(dim_list_all,sets{iset}.dim_list);
    end
    if length(dim_list_all)~=length(1:max(dim_list_all))
        disp(sprintf('some dimensions are missing in set %1.0f',iset))
    end
end
pcon_dim_max=getinp('maximum dimension for consensus alignment','d',[1 max(dim_list_all)],pcon_dim_max);
opts_pcon.allow_scale=getinp('1 to allow scaling for consensus','d',[0 1],opts_pcon.allow_scale);
%
consensus=cell(pcon_dim_max,1);
z=cell(pcon_dim_max,1);
znew=cell(pcon_dim_max,1);
ts=cell(pcon_dim_max,1);
details=cell(pcon_dim_max,1);
opts_pcon_used=cell(pcon_dim_max,1);
%
opts_pcon.overlaps=ovlp_array;
%
d_knitted=cell(pcon_dim_max,1);
ds_components=cell(1,nsets); %partial datasets, aligned via Procrustes
%
disp('overlap matrix')
disp(ovlp_array'*ovlp_array);
%
for ip=1:pcon_dim_max
    z{ip}=zeros(nstims_all,ip,nsets);
    for iset=1:nsets
        z{ip}(:,:,iset)=ds_align{iset}{ip};
    end
    [consensus{ip},znew{ip},ts{ip},details{ip},opts_pcon_used{ip}]=procrustes_consensus(z{ip},opts_pcon);
    disp(sprintf(' did Procrustes consensus for dim %1.0f, iterations: %4.0f, final total rms dev: %8.5f',...
        ip,length(details{ip}.rms_change),sqrt(sum(details{ip}.rms_dev(:,end).^2))));
    d_knitted{ip}=consensus{ip};
    for iset=1:nsets
        ds_components{iset}{1,ip}=znew{ip}(:,:,iset);
    end
end
%
%find the ray descriptors
%
[rays_knitted,opts_rays_knitted_used]=psg_findrays(sa_pooled.btc_specoords,opts_rays_used{1}); %ray parameters based on first dataset; finding rays only depends on metadata
%
[sets_nonan,ds_nonan,sas_nonan,opts_nonan_used]=psg_remnan_coordsets(sets_align,ds_components,sas_align,ovlp_array,opts_nonan); %remove the NaNs
%find the rays for sets with nan's removed (since the order has been changed) and use these to plot
rays_nonan=cell(nsets,1);
for iset=1:nsets
    [rays_nonan{iset},opts_rays_nonan_used{iset}]=psg_findrays(sas_nonan{iset}.btc_specoords,opts_rays_used{iset});
end
disp('created ray descriptors for knitted and nonan datasets');
%
%plot knitted data and individual sets
dim_con=1;
while max(dim_con)>0
    dim_con=getinp('consensus dimension to plot (0 to end)','d',[0 pcon_dim_max]);
    if max(dim_con)>0
        dims_to_plot=getinp('dimensions to plot','d',[1 dim_con]);
        tstring=sprintf('consensus dim %1.0f [%s]',dim_con,sprintf('%1.0f ',dims_to_plot));
        opts_plot=struct;
        figure;
        set(gcf,'Position',[100 100 1200 800]);
        set(gcf,'Name',cat(2,'knitted ',tstring));
        set(gcf,'NumberTitle','off');
        opts_plot_used=psg_plotcoords(d_knitted{dim_con},dims_to_plot,sa_pooled,rays_knitted,opts_plot);
        axis equal;
        axis vis3d;
        xlims=get(gca,'XLim');
        ylims=get(gca,'YLim');
        zlims=get(gca,'ZLim');
        axes('Position',[0.01,0.05,0.01,0.01]); %for text
        text(0,0,cat(2,'knitted ',tstring),'Interpreter','none','FontSize',10);
        axis off;
        %
        opts_rays_nonan_used=cell(nsets,1);
        opts_plot_nonan_used=cell(nsets,1);
        for iset=1:nsets
            tstringc=sprintf(' component set %1.0f: %s, %s',iset,sets{iset}.label);
            figure;
            set(gcf,'Position',[100 100 1200 800]);
            set(gcf,'Name',cat(2,tstringc,' ',tstring));
            set(gcf,'NumberTitle','off');
            opts_plot_nonan_used{iset}=psg_plotcoords(ds_nonan{iset}{dim_con},dims_to_plot,sas_nonan{iset},rays_nonan{iset},opts_plot);
            axis equal
            axis vis3d
            set(gca,'XLim',xlims);
            set(gca,'YLim',ylims);
            set(gca,'ZLim',zlims);
            axes('Position',[0.01,0.05,0.01,0.01]); %for text
            text(0,0,cat(2,tstringc,' ',tstring),'Interpreter','none','FontSize',10);
            axis off;
        end
        %plot knitted with compoennts
        figure;
        set(gcf,'Position',[100 100 1200 800]);
        set(gcf,'Name',cat(2,'composite ',tstring));
        set(gcf,'NumberTitle','off');
        %
        opts_plot_components=cell(1,nsets);
        opts_plot_components=cell(1,nsets);
        %need to customize this, and then plot, onsame axes, each component
        %using color_order
        opts_plot_knitted=struct;
        opts_plot_knitted.marker_noray='';
        opts_plot_knitted.color_origin='k';
        opts_plot_knitted.color_nearest_nbr='k';
        opts_plot_knitted.noray_connect=0;
        opts_plot_knitted_used=psg_plotcoords(d_knitted{dim_con},dims_to_plot,sa_pooled,setfield(rays_knitted,'nrays',0),opts_plot_knitted);
        for iset=1:nsets
            opts_plot_components{iset}=opts_plot_knitted;
            opts_plot_components{iset}.axis_handle=opts_plot_knitted_used.axis_handle;
            pcolor=color_list(1+mod(iset-1,length(color_list)));
            opts_plot_components{iset}.color_origin=pcolor;
            opts_plot_components{iset}.color_nearest_nbr=pcolor;
            opts_plot_components_used=psg_plotcoords(ds_nonan{iset}{dim_con},dims_to_plot,sas_nonan{iset},...
                setfield(rays_nonan{iset},'nrays',0),opts_plot_components{iset});
        end
        axis equal
        axis vis3d
        set(gca,'XLim',xlims);
        set(gca,'YLim',ylims);
        set(gca,'ZLim',zlims);
        axes('Position',[0.01,0.05,0.01,0.01]); %for text
        text(0,0,cat(2,'composite ',tstring),'Interpreter','none','FontSize',10);
        axis off;
    end
end
if getinp('1 to write a file with knitted coordinate data and metadata','d',[0 1])
    opts_write=struct;
    opts_write.data_fullname_def='[paradigm]pooled_coords_ID.mat';
    opts_write_used=psg_write_coorddata([],d_knitted,setfield([],'stim_labels',strvcat(sa_pooled.typenames)),opts_write);
    %
    metadata_fullname_def=opts_write_used.data_fullname;
    metadata_fullname_def=metadata_fullname_def(1:-1+min(find(cat(2,metadata_fullname_def,'_')=='_')));
    if isfield(sa_pooled,'nsubsamp')
        metadata_fullname_def=cat(2,metadata_fullname_def,sprintf('%1.0f',sa_pooled.nsubsamp));
    end
    metadata_fullname=getinp('metadata file name','s',[],metadata_fullname_def);
    s=sa_pooled;
    save(metadata_fullname,'s');
end
%pipeline field needs to be added to data file
%example of pipeline field:
% 
% pipeline
% pipeline =
%   1×1 cell array
%     {1×1 struct}
% pipeline{1}
% ans = 
%   struct with fields:
% 
%     type: 'qform2coord'
%     sets: [1×1 struct]
%     opts: [1×1 struct]
% pipeline{1}.sets
% ans = 
%   struct with fields:
% 
%           type: 'qform'
%       dim_list: [1 2 3 4 5 6 7 8 9 10]
%         nstims: 25
%     label_long: './psg_data/bc6pt9.mat c:aug q:../stim/btc_allraysfixedb_avg_100surrs_madj.mat m:qform'
%          label: 'psg_data/bc6pt9 c:aug q:avg_madj m:qform'
% pipeline{1}.opts
% ans = 
%   struct with fields:
% 
%      opts_read_used: [1×1 struct]
%     opts_qpred_used: [1×1 struct]
% pipeline{1}.opts.opts_read_used
% ans = 
%   struct with fields:
% 
%                 input_type: 2
%                     if_log: 1
%               if_justsetup: 1
%          data_fullname_def: './psg_data/bgca3pt_coords_BL_sess01_04.mat'
%         setup_fullname_def: './psg_data/bgca3pt9.mat'
%                permutes_ok: 1
%      faces_mpi_atten_indiv: 1
%        faces_mpi_atten_age: 1
%     faces_mpi_atten_gender: 1
%        faces_mpi_atten_emo: 1
%        faces_mpi_atten_set: 0.2000
%                   permutes: [1×1 struct]
%            permute_raynums: []
%              data_fullname: []
%             setup_fullname: './psg_data/bc6pt9.mat'
%                 type_class: 'btc'
%            input_type_desc: 'qform model'
% pipeline{1}.opts.opts_qpred_used
% ans = 
%   struct with fields:
% 
%     qform_datafile_def: '../stim/btc_allraysfixedb_avg_100surrs_madj.mat'
%        qform_modeltype: 12
%                 if_log: 0
%        negeig_makezero: 1
%                 negtol: 1.0000e-07
%        if_pca_centroid: 1
%                 var_ex: [62.0821 110.7713 113.8419 113.8419 113.8419 113.8419 113.8419 113.8419 113.8419 113.8419]
%                var_tot: 113.8419
%             predcoords: [25×10 double]
% 

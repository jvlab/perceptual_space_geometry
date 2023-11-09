%psg_consensus_demo: demonstrate visualization of psg coordinates after
%Procrustes consensus calculation; coordinaes can be experimental or predicted
%
% Main results in results_consensus.
%
% Also demonstrates plotting with offset, plotting with pca rotation,
% plotting with connecting corresponding points, and plotting rings.
%
% 18Jun23: add question for if_nearest_neighbor
% 31Jul23: add option to define rays by single points (so faces data will be plotted correctly)
% 26Sep23: remove dim_list_all
% 04Oct23: add more flexible offset specifications and line width specifications
% 31Oct23: sanity check for consensus inputs
%
%  See also: PSG_GET_COORDSETS, PSG_FINDRAYS, PSG_QFORMPRED, PSG_PLOTCOORDS, PSG_VISUALIZE_DEMO, PROCRUSTES,
%    PSG_COLORS_LEGACY, PROCRUSTES_CONSENSUS, PSG_PROCRUSTES_DEMO.
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
if ~exist('plotformats')
    plotformats=[2 2;3 2;3 3;4 3;4 4;5 3]; %dim model, number of dims at a time
end
%
if_spray=getinp('1 to define rays by single points, suppress ray angle calculations, and suppress ray angle plots','d',[0 1],0);
if (if_spray)
    opts_rays.ray_minpts=1;
    if_plotrays=0;
end
%
[sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred);
nsets=length(sets);
nstims=size(ds{1}{1},1);
%
% below here, code differs c/w psg_procrustes_demo.
%
if ~exist('opts_pcon') opts_pcon=struct(); end %for procrustes_consensus
%
% choose datasets and dimensions
%
disp('datasets available:');
for iset=1:nsets
    disp(sprintf(' set %2.0f: dim range [%3.0f %3.0f] label: %s',iset,min(sets{iset}.dim_list),max(sets{iset}.dim_list),sets{iset}.label));
end
ncons=0;
cons_table=zeros(0,3); %dataset, dimension choice, and pointer in dim_list to dimension choice
file_strings='';
if_ok=0;
while if_ok==0
    for icons=1:ncons
        disp(sprintf('dataset %1.0f for consensus calc is %2.0f dim model of %s',...
            icons,cons_table(icons,2),sets{icons}.label));
    end
    choice_min=-1;
    if (ncons==0)
        choice_min=1;
    end
    choice=getinp(sprintf('dataset number to use as dataset %1.0f for consensus calculation, 0 to end, -1 to restart',ncons+1),'d',[choice_min nsets]);
    if (choice>0)
        ncons=ncons+1;
        cons_table(ncons,1)=choice;
        cons_table(ncons,2)=getinp(sprintf('dimension to use from dataset %1.0f (%1.0f)',ncons,cons_table(ncons,1)),'d',...
            [1 max(sets{cons_table(ncons,1)}.dim_list)]);
        if ismember(cons_table(ncons,2),sets{cons_table(ncons,1)}.dim_list)
            cons_table(ncons,3)=find(cons_table(ncons,2)==sets{cons_table(ncons,1)}.dim_list);
            if (ncons>1)
                file_strings=cat(2,file_strings,'+');
            end
            file_strings=cat(2,file_strings,sets{ncons}.label);
        else
            disp('this model dimension is not available.')
            ncons=ncons-1;
        end
    elseif choice==-1
        ncons=0;
        cons_table=zeros(0,3);
        file_strings='';
    else
        if_ok=1;
    end
end
%
cons_maxdim=max(cons_table(:,2));
z=zeros(nstims,cons_maxdim,ncons);
%retrieve datasets for consensus and pad coords if needed
for icons=1:ncons
    z(:,[1:cons_table(icons,2)],icons)=ds{cons_table(icons,1)}{cons_table(icons,3)};
end
%do consensus calculation without and with scaling allowed
results_consensus=struct;
initialize_set=min(find(max(cons_table(:,2)==cons_table(:,2)))); %use dataset with greatest dimensions for initialization
opts.pcon_initialize_set=initialize_set;
for allow_scale=0:1
    [consensus,znew,ts,details,opts_pcon_used]=procrustes_consensus(z,setfield(opts_pcon,'allow_scale',allow_scale));
    disp(sprintf('consensus calculation done, allow_scale=%1.0f, %4.0f iterations, final rms change %12.8f',...
        allow_scale,length(details.rms_change),details.rms_change(end)));
    results_consensus.consensus(:,:,1+allow_scale)=consensus;
    results_consensus.znew(:,:,:,1+allow_scale)=znew;
    results_consensus.ts{1+allow_scale}=ts;
    results_consensus.details{1+allow_scale}=details;
    results_consensus.opts_pcon_used{1+allow_scale}=opts_pcon_used;
end
%create lib{} with consensus, original datasets, and original datasets aligned to consensus
lib=cell(1,2+3*ncons); %consensus with and without scaling, then each dataset to be adjusted, then adjusted without and with scaling allowed
for allow_scale=0:1
    if (allow_scale==0)
        scale_label='consensus, no scaling';
    else
        scale_label='consensus, with scaling';
    end
    lib{1+allow_scale}.name=scale_label;
    lib{1+allow_scale}.d=results_consensus.consensus(:,:,1+allow_scale);
    lib{1+allow_scale}.sa=sas{initialize_set};
    lib{1+allow_scale}.rays=rayss{initialize_set};
    lib{1+allow_scale}.opts_vis=setfield(opts_vis,'file_string',cat(2,scale_label,'<-{',file_strings,'}'));
end
for icons=1:ncons
    for allow_scale=0:2 %2 is original library
         switch allow_scale
            case 0 %scale not allowed in adjustment
                coords=results_consensus.znew(:,:,icons,1+allow_scale);
                scale_label='adj, no scaling';
            case 1 %scale allowed in adjustment
                coords=results_consensus.znew(:,:,icons,1+allow_scale);
                scale_label=sprintf('adj, sc=%7.4f',results_consensus.details{1+allow_scale}.ts_cum{end}{icons}.scaling);
            case 2 %unadjusted
                coords=z(:,:,icons);
                scale_label='unadj';
        end
        lib{icons*3+allow_scale}.name=...
            sprintf('dataset %1.0f dim %2.0f %s',cons_table(icons,1),cons_table(icons,2),scale_label);
        lib{icons*3+allow_scale}.sa=sas{icons}; %changed from iset, 08Nov23
        lib{icons*3+allow_scale}.rays=rayss{icons}; %changed from iset, 08Nov23
        lib{icons*3+allow_scale}.d=coords;
        lib{icons*3+allow_scale}.opts_vis=setfield(opts_vis,'file_string',cat(2,sprintf(dfmt,sets{cons_table(icons,1)}.label,cons_table(icons,2)),scale_label));
    end %allow_scale
end %icons
%
%if requested, plot one or more datasets together
%
if_ok_plot=0;
while if_ok_plot==0
    disp('selections for plotting');
    for il=1:length(lib)
        disp(sprintf('%2.0f->%35s (%s)',il,lib{il}.name,lib{il}.opts_vis.file_string));
    end
    lib_list=getinp('choice(s), 0 to end','d',[0 length(lib)]);
    if min(lib_list)==0
        if_ok_plot=1;
    else
        nlib=length(lib_list);
        dm=cell(0);
        sam=cell(0);
        raysm=cell(0);
        opts_vism=cell(0);
        dims_lib=[];
        for il=1:nlib
            dm{il}{1}=lib{lib_list(il)}.d;
            sam{il}=lib{lib_list(il)}.sa;
            raysm{il}=lib{lib_list(il)}.rays;
            opts_vism{il}=lib{lib_list(il)}.opts_vis;
            dims_lib(il)=size(dm{il}{1},2);
        end
        nconds=size(dm{1}{1},1);
        %default offset pointer points to the origin
        offset_ptr=find(raysm{1}.whichray==0);
        if length(offset_ptr)~=1
            offset_ptr=0;
        end
        opts_plotm=struct;
        if isfield(opts_plot,'colors')
            opts_plotm.colors=opts_plot.colors;
        end
        opts_multm=struct;
        opts_multm.if_fit_range=double(nsets>1);
        if_replot=1;
        line_widths=[1:nlib];
        while (if_replot==1)
            opts_multm.if_fit_range=getinp('1 to fit range, 0 for standard range','d',[0 1],opts_multm.if_fit_range);
            line_widths=getinp(sprintf('line widths (%3.0f values)',nlib),'d',[1 10],line_widths);
            line_widths=line_widths(1+mod([0:nlib-1],nlib));
            if_pcrot=getinp('1 to apply pc rotations','d',[0 1]);
            if (if_pcrot)
                offset_norot=getinp('1 to avoid applying pcrot to offset, -1 to apply only to data','d',[-1 1],-1);
                if_pcrot_whichuse=getinp('which dataset to use for pcs (0 to rotate each separately)','d',[0 nlib]);
            else
                offset_norot=0;
                if_pcrot_whichuse=0;
            end
            if_offset=getinp('dimension for offsets (0 for none, -1 to specify each vector)','d',[-1 max(dims_lib)]);
            if (if_offset>0)
                offset_size=getinp('offset size','f',[0.1 10],1);
            elseif (if_offset==0)
                offset_size=0;
            else
                offset_vecs=zeros(nlib,max(dims_lib));
                for il=1:nlib
                   disp(sprintf('for plot %2.0f->%35s (%s)',lib_list(il),lib{lib_list(il)}.name,lib{lib_list(il)}.opts_vis.file_string));
                   offset_vecs(il,:)=getinp(sprintf('vector of length %1.0f to be subtracted',max(dims_lib)),'f',[-Inf Inf],zeros(1,max(dims_lib)));
                end
            end
            offset_ptr=getinp('condition to plot at 0 (-1 for centroid, 0 for none)','d',[-1 nconds]);
            if_connect=getinp('1 for star connection, 2 for circuit, 3 for custom, 0 for none, negative to only show connections','d',[-3 3]);
            opts_multm.connect_only=double(if_connect<0);
            switch if_connect
                case 0
                    opts_multm.connect_specs=[];
                case {-1,1}
                    opts_multm.connect_specs='star';
                case {-2,2}
                    opts_multm.connect_specs='circuit';
                case {-3,3}
                    connect_spec_pairs=getinp('a sequence of pairs, e.g., [1 2 2 3] for [1 2;2 3]','d',[1 nlib],[1 2]);
                    opts_multm.connect_specs=reshape(connect_spec_pairs,2,length(connect_spec_pairs)/2)';
            end
            if (if_connect~=0)
                opts_multm.connect_line_width=getinp('width for connection line','d',[1 4],1);
                opts_multm.connect_line_type=getinp('style for connection line (non-neg directions only, neg dirs will be --)','s',[],'-');
            end
            opts_plotm.if_rings=getinp('1 to plot rings','d',[0 1],0);
            opts_plotm.if_nearest_neighbor=getinp('1 to connect nearest neighbors, 0 not, -1 if unassigned points','d',[-1 1],-1);
            for il=1:nlib
                opts_vism_use{il}=opts_vism{il};
                if if_offset>=0
                    opts_vism_use{il}.offset=offset_size*double([1:max(dims_lib)]==if_offset)*(il-mean([1:nlib]));
                else
                    opts_vism_use{il}.offset=offset_vecs(il,:);                    
                end
                opts_vism_use{il}.offset_ptr=offset_ptr;
                opts_vism_use{il}.if_pcrot=if_pcrot;
                opts_vism_use{il}.offset_norot=offset_norot;
            end
            opts_multm.if_pcrot_whichuse=if_pcrot_whichuse;
            opts_multm.line_widths=line_widths;
            [opts_vism_used,opts_plotm_used,opts_multm_used]=psg_visualize(plotformats,dm,sam,raysm,opts_vism_use,opts_plotm,opts_multm);
            if_replot=getinp('1 to replot','d',[0 1]);
        end %if_replot
    end
end %if_ok_plot

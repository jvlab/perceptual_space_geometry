% psg_align_disp_enh_btc_demo_rs: demonstrate alignment and display of datasets across human subjects and DreamSim (per Muya Zhou),
% with offsets and other enhanced plotting options
%
% Based on rs_align_disp_enh_coordsets_demo.
%
%%section to force btc defaults, even if rs_aux_defaults.mat has been created or modified
if ~exist('aux_force_filename') aux_force_filename='rs_aux_defaults_btc.mat'; end
auxs_force=struct;
opts_needed={'opts_read','opts_rays','opts_check','opts_align','opts_qpred','opts_knit','opts_disp','opts_disp_enh'};
for k=1:length(opts_needed)
    auxs_force.(opts_needed{k})=rs_aux_force(opts_needed{k},[],aux_force_filename);
end
%
filename_paradigms{1}={...
        '../psg/psg_data/bgca3pt_coords_DS_sess01_10.mat',... 
        '../psg/psg_data/bgca3pt_coords_BL_sess01_10.mat',... 
        '../psg/psg_data/bgca3pt_coords_MC_sess01_10.mat',... 
        '../psg/psg_data/bgca3pt_coords_NF_sess01_10.mat',...
        '../psg/psg_data/bgca3pt_coords_SAW_sess01_10.mat',...
        '../psg/psg_data/bgca3pt_coords_SN_sess01_10.mat',...
        '../psg/psg_data/bgca3pt_coords_ZK_sess01_10.mat'};
filename_paradigms{2}={...
        '../psg/psg_data/bc6pt_coords_DS_sess01_10.mat',... 
        '../psg/psg_data/bc6pt_coords_BL_sess01_10.mat',... 
        '../psg/psg_data/bc6pt_coords_CME_sess01_10.mat',... 
        '../psg/psg_data/bc6pt_coords_MC_sess01_10.mat',...
        '../psg/psg_data/bc6pt_coords_SAW_sess01_10.mat',...
        '../psg/psg_data/bc6pt_coords_ZK_sess01_10.mat'};
filename_paradigms{3}={...
        '../psg/psg_data/bcpm3pt_coords_DS_sess01_10.mat',... 
        '../psg/psg_data/bcpm3pt_coords_BL_sess01_10.mat',... 
        '../psg/psg_data/bcpm3pt_coords_MC_sess01_10.mat',...
        '../psg/psg_data/bcpm3pt_coords_SAW_sess01_10.mat',...
        '../psg/psg_data/bcpm3pt_coords_ZK_sess01_10.mat'};
filename_paradigms{4}={...
        '../psg/psg_data/bdce3pt_coords_DS_sess01_10.mat',... 
        '../psg/psg_data/bdce3pt_coords_BL_sess01_10.mat',... 
        '../psg/psg_data/bdce3pt_coords_MC_sess01_10.mat',...
        '../psg/psg_data/bdce3pt_coords_SAW_sess01_10.mat',...
        '../psg/psg_data/bdce3pt_coords_SN_sess01_10.mat',...
        '../psg/psg_data/bdce3pt_coords_ZK_sess01_10.mat'};
filename_paradigms{5}={...
        '../psg/psg_data/tvpm3pt_coords_DS_sess01_10.mat',... 
        '../psg/psg_data/tvpm3pt_coords_BL_sess01_10.mat',... 
        '../psg/psg_data/tvpm3pt_coords_MC_sess01_10.mat',...
        '../psg/psg_data/tvpm3pt_coords_NF_sess01_10.mat',...
        '../psg/psg_data/tvpm3pt_coords_SAW_sess01_10.mat',...
        '../psg/psg_data/tvpm3pt_coords_SN_sess01_10.mat',...
        '../psg/psg_data/tvpm3pt_coords_ZK_sess01_10.mat'};
filename_paradigms{6}={...
        '../psg/psg_data/detv3pt_coords_DS_sess01_10.mat',... 
        '../psg/psg_data/detv3pt_coords_MC_sess01_10.mat',... 
        '../psg/psg_data/detv3pt_coords_NF_sess01_10.mat',... 
        '../psg/psg_data/detv3pt_coords_ZK_sess01_10.mat'};
filename_paradigms{7}={...
        '../psg/psg_data/gtva3pt_coords_DS_sess01_10.mat',... 
        '../psg/psg_data/gtva3pt_coords_MC_sess01_10.mat',... 
        '../psg/psg_data/gtva3pt_coords_NF_sess01_10.mat',... 
        '../psg/psg_data/gtva3pt_coords_ZK_sess01_10.mat'};
filename_paradigms{8}={...
        '../psg/psg_data/btuv3pt_coords_DS_sess01_10.mat',... 
        '../psg/psg_data/btuv3pt_coords_MC_sess01_10.mat',... 
        '../psg/psg_data/btuv3pt_coords_NF_sess01_10.mat',... 
        '../psg/psg_data/btuv3pt_coords_ZK_sess01_10.mat'};
filename_paradigms{9}={...
        '../psg/psg_data/dgea3pt_coords_DS_sess01_10.mat',... 
        '../psg/psg_data/dgea3pt_coords_BL_sess01_10.mat',... 
        '../psg/psg_data/dgea3pt_coords_MC_sess01_10.mat',...
        '../psg/psg_data/dgea3pt_coords_NF_sess01_10.mat',...
        '../psg/psg_data/dgea3pt_coords_SN_sess01_10.mat',...
        '../psg/psg_data/dgea3pt_coords_ZK_sess01_10.mat'};
%
nparas=length(filename_paradigms);
ncgps=1; %number of coordinate groups
label_maxlength=6; %max length of a stimulus label
aux_outs=cell(nparas,1);
%
dim_select=3;
margin_amounts=ones(nparas,dim_select);
margin_amounts(1,:)=repmat(2,1,dim_select);
margin_amounts(3,:)=repmat(2,1,dim_select);
if ~exist('para_list') para_list=[1:nparas]; end
%
plot_marg=0.5; %global plot margin
for ipara_ptr=1:length(para_list)
    ipara=para_list(ipara_ptr);
    filenames=filename_paradigms{ipara};
    nfiles=length(filenames);
    aux_in=auxs_force;
    aux_in.opts_read=setfields(aux_in.opts_read,{'input_type','if_auto','if_log'},{1,1,0});
    aux_in.nsets=nfiles;
    disp(sprintf(' group %1.0f: %2.0f files',ipara,nfiles))
    [data_read,aux_read]=rs_get_coordsets(filenames,aux_in);
    paradigm_name=data_read.sets{1}.paradigm_name;
    if_ok=1;
    for ifile=1:nfiles
        rays=aux_read.rayss{ifile};
        if isempty(rays)
            disp(sprintf(' file %1.0f: %70s: ray structure not created',ifile,filenames{ifile}))
            if_ok=0;
        else
            disp(sprintf(' file %1.0f: %70s: %3.0f rays, %3.0f rings, %3.0f pairs',ifile,filenames{ifile},rays.nrays,rays.nrings,rays.npairs))
        end
    end
    if (if_ok)
        %align data, rotate to consensus, and rotate consensus into pca coords
        aux_align_def=auxs_force.opts_align;
        aux_align_def.opts_align.if_log=0;
        [data_align,aux_align]=rs_align_coordsets(data_read,aux_align_def);
        aux_knit_def=auxs_force.opts_knit;
        aux_knit_def.opts_knit.if_log=0;
        aux_knit_def.opts_knit.if_pca=1; %rotate to PCA
        [data_consensus,aux_knit]=rs_knit_coordsets(data_align,aux_knit_def);
        data_disp=aux_knit.components;
        rays_use=aux_knit.rayss{1};
        %
        hfig=figure;
        opts_disp=auxs_force.opts_disp;
        opts_disp.fig_handle=hfig;
        opts_disp.fig_position=[50 80 1200 800];
        %
        set_label_string=[];
        opts_disp.fig_name=sprintf('group %1.0f: %s',ipara,paradigm_name);
        for ifile=1:nfiles
            opts_disp.set_labels{ifile}=data_read.sets{ifile}.subj_id;
            set_label_string=cat(2,set_label_string,data_read.sets{ifile}.subj_id,' ');
        end
        opts_disp.data_label_font_size=7;
        opts_disp.axis_label_prefix='pc';
        opts_disp.dim_select=dim_select;
        opts_disp.coord_group_method='keeplow';
        %
        opts_disp.set_offsets='margin_amount';
        opts_disp.set_offsets_margin_amount=margin_amounts(ipara,:);
        opts_disp.set_offsets_coordchoices=2;
        %
        opts_disp_enh=auxs_force.opts_disp_enh;
        opts_disp_enh.if_points=0;
        opts_disp_enh.if_rays=1;
        opts_disp_enh.if_rings=0;
        opts_disp_enh.if_nbrs=0;
        %
        aux_outs{ipara}=rs_disp_enh_coordsets(data_disp,setfields(struct,{'opts_disp','opts_disp_enh'},{opts_disp,opts_disp_enh}),rays_use);
        %
        %make tickmars at 2 JND, no labels
        xlims=plot_marg*[-1 1]+[floor(min(get(gca,'XLim'))),ceil(max(get(gca,'XLim')))];
        set(gca,'XLim',xlims);
        set(gca,'XTick',[xlims(1):2:xlims(2)]);
        set(gca,'XTickLabel',[]);
        %
        ylims=plot_marg*[-1 1]+[floor(min(get(gca,'Ylim'))),ceil(max(get(gca,'Ylim')))];
        set(gca,'Ylim',ylims);
        set(gca,'YTick',[ylims(1):2:ylims(2)]);
        set(gca,'YTickLabel',[]);
        %
        zlims=plot_marg*[-1 1]+[floor(min(get(gca,'ZLim'))),ceil(max(get(gca,'ZLim')))];
        set(gca,'ZLim',xlims);
        set(gca,'ZTick',[xlims(1):2:xlims(2)]);
        set(gca,'ZTickLabel',[]);
        %
        axes('Position',[0.01,0.02,0.01,0.01]); %for text
        text(0,0,cat(2,paradigm_name,' ',set_label_string,'(ticks: 2 JND),',data_read.sets{1}.subj_id,' is labelled'),'Interpreter','none','FontSize',8);
        axis off;
        %
    end %if_ok
end %ipara

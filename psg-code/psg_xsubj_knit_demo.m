%psg_xsubj_knit_demo: consistency across subjects, within individual sets of coordinates
% typically with four axes explored in each
%
% results saved in r
%
% 25Nov25: add option to augment the dimension of the knitted set, compare to its constituents
% 29Nov25: add if_pca, better labeling for plots
%
%  See also:  RS_READ_COORDSETS, RS_KNIT_COORDSETS, RS_KNIT_COORDSETS_DEMO, PSG_PCAOFFSET.
%
coord_string=getinp('coordinate string, e.g., bgca3pt','s',[],'bgca3pt');
d_max=getinp('max embedding dimension to analyze','d',[2 10],7);
if_pca=getinp('1 to do initial pca within each dimension','d',[0 1]);
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
nshuffs=getinp('number of shuffles','d',[0 1000],100);
if (if_frozen~=0) 
    rng('default');
    if (if_frozen<0)
        rand(1,abs(if_frozen));
    end
else
    rng('shuffle');
end
opts_read=struct;
filenames=[];
nsets=length(filenames);
opts_read.if_gui=1;
opts_read.if_log=0;
opts_read.if_warn=1;
opts_read.if_auto=1;
opts_read.if_symaug=0;
opts_read.input_type=1;
opts_read.ui_filter=cat(2,'./psg_data/',coord_string,'_coords_*_sess01_10.mat');
%
opts_align=struct;
opts_align.min=1;
opts_align.if_log=0;
%
opts_knit=struct;
opts_knit.if_log=1;
opts_knit.allow_reflection=1;
opts_knit.allow_offset=1;
opts_knit.allow_scale=0;
opts_knit.if_normscale=1; %nondefault but irrelevant if allow_scale=0
opts_knit.if_pca=1; %nondefault
opts_knit.if_stats=1;
opts_knit.nshuffs=nshuffs;
opts_knit.dim_max_in=d_max;
%
aux_read=struct;
aux_read.opts_read=opts_read;
aux_read.nsets=nsets;
%
aux_align=struct;
aux_align.opts_align=opts_align;
%
aux_knit=struct;
aux_knit.opts_knit=opts_knit;
%
%read the data
[data_read,aux_read_out]=rs_get_coordsets(filenames,aux_read);
nsets=length(data_read.sets);
if if_pca
    for iset=1:nsets
        nds=length(data_read.ds{iset});
        for idim=1:nds
            coords=data_read.ds{iset}{idim};
            coords=coords-repmat(mean(coords,1),size(coords,1),1); %subtract centroid
            data_read.ds{iset}{idim}=psg_pcaoffset(coords); %replace by pca
        end
        disp(sprintf('pca and centering applied for dataset %2.0f',iset));
    end
end
%align
[data_align,aux_align_out]=rs_align_coordsets(data_read,aux_align);
%
%shorten the label for plotting
subj_ids=cell(0);
subj_id_string='';
for k=1:length(data_align.sets)
    subj_id=data_align.sets{k}.subj_id;
    subj_ids{k}=subj_id;
    data_align.sets{k}.label=cat(2,data_align.sets{k}.paradigm_name,':',subj_id);
    if ~contains(subj_id_string,subj_id)
        if ~isempty(subj_id_string) 
            subj_id_string=cat(2,subj_id_string,'+');
        end
        subj_id_string=cat(2,subj_id_string,subj_id);
    end
end
%knit
[data_knit,aux_knit_out]=rs_knit_coordsets(data_align,aux_knit);
set(gcf,'Name',sprintf('%s (%s): no scaling',coord_string,subj_id_string));
set(gcf,'Position',[50 100 1500 700]);
axes('Position',[0.25,0.04,0.01,0.01]); %for text
text(0,0,sprintf('if_pca: %1.0f',if_pca),'Interpreter','none','FontSize',8);
axis off;

%also knit with allowing a scaling between datasets
aux_knit2=aux_knit;
aux_knit2.opts_knit.allow_scale=1;
aux_knit2.opts_knit.if_normscale=1;
[data_knit2,aux_knit2_out]=rs_knit_coordsets(data_align,aux_knit2);
set(gcf,'Name',sprintf('%s (%s): scaling+norm',coord_string,subj_id_string));
set(gcf,'Position',[50 100 1500 700]);
axes('Position',[0.25,0.04,0.01,0.01]); %for text
text(0,0,sprintf('if_pca: %1.0f',if_pca),'Interpreter','none','FontSize',8);
axis off;
%
r=struct;
r.data_read=data_read;
r.nsets=nsets;
r.if_pca=if_pca;
r.coord_string=coord_string;
r.subj_ids=subj_ids;
r.nshuffs=nshuffs;
r.knit_stats_setup=aux_knit_out.knit_stats_setup;
r.knit_stats_noscale=aux_knit_out.knit_stats;
r.knit_stats_scale_renornm=aux_knit2_out.knit_stats;

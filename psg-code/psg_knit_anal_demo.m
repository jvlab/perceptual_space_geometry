%psg_knit_anal_demo: participation_ratio analysis of knitted-together experiments
% typically with four axes explored in each
%
% Has built-in selection of sets to combine and colors
% Allows for adding Gaussian noise
%
% results saved in r
%
%  See also:  RS_READ_COORDSETS, RS_KNIT_COORDSETS.
%
if ~exist('subset_defs')
    subset_defs{1}.contains=cell(0);
    subset_defs{1}.name='all';
    subset_defs{1}.color=[0 0 0];
    %
    subset_defs{2}.contains={'g'};
    subset_defs{2}.name='has gamma';
    subset_defs{2}.color=[0.7 0.7 0.7];
    %
    subset_defs{3}.contains={'b','c'};
    subset_defs{3}.name='has beta-card';
    subset_defs{3}.color=[0 0 1];
    %
    subset_defs{4}.contains={'d','e'};
    subset_defs{4}.name='has beta-diag';
    subset_defs{4}.color=[0 1 0];
    %
    subset_defs{5}.contains={'t','u','v','w'};
    subset_defs{5}.name='has theta';
    subset_defs{5}.color=[0.5 0.5 0];
    %
    subset_defs{6}.contains={'a'};
    subset_defs{6}.name='has alpha';
    subset_defs{6}.color=[1 0 0];
    %
end
if ~exist('btc_colorspec')
    btc_colorspec=struct;
    btc_colorspec.bcpm=[0.0 0.2 0.8];
    btc_colorspec.bdce=[0.0 0.7 0.8];
    btc_colorspec.bgca=[1.0 0.0 0.2];
    btc_colorspec.btuv=[0.5 0.5 0.2];   
    btc_colorspec.detv=[0.2 0.8 0.0];
    btc_colorspec.dgea=[0.7 0.3 0.0];
    btc_colorspec.gtva=[0.2 0.2 0.0];
    btc_colorspec.tvpm=[0.8 0.8 0.0];
end
%
subj_id=upper(getinp('subject ID','s',[],'MC'));
d_max=getinp('max embedding dimension to analyze','d',[2 10],7);
sym_list_avail=psg_btcsyms();
for k=1:length(sym_list_avail)
    disp(sprintf(' %1.0f-> symmetry %s',k,sym_list_avail{k}));
end
sym_ptr_list=getinp('symmetries to apply','d',[1 length(sym_list_avail)],[1:length(sym_list_avail)]);
nsymaugs=length(sym_ptr_list);
sym_list=cell(1,nsymaugs);
for k=1:nsymaugs
    sym_list{k}=sym_list_avail{sym_ptr_list(k)};
end
noise_jit=getinp('amount of noise to add, in rms jnd''s','f',[0 5],0);
if noise_jit>0
    if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
    if (if_frozen~=0) 
        rng('default');
        if (if_frozen<0)
            rand(1,abs(if_frozen));
        end
    else
        rng('shuffle');
    end
end
opts_read=struct;
opts_read.if_gui=1;
opts_read.if_log=0;
opts_read.if_warn=1;
opts_read.if_auto=1;
opts_read.if_symaug=0;
opts_read.input_type=1;
%
opts_align=struct;
opts_align.if_log=0; %non-default
%
opts_knit=struct;
opts_knit.if_log=0; %non-default 
opts_knit.allow_reflection=1;
opts_knit.allow_offset=1;
opts_knit.allow_scale=0;
opts_knit.if_normscale=1; %nondefault but irrelevant if allow_scale=0
opts_knit.if_pca=1; %nondefault
opts_knit.if_stats=0;
opts_knit.dim_max_in=d_max;
%
if strcmp(upper(subj_id),'QFM')
    opts_read.ui_filter=cat(2,'./psg_data/*3pt_coords_',subj_id,'_sess01_01.mat');
else
    opts_read.ui_filter=cat(2,'./psg_data/*3pt_coords_',subj_id,'_sess01_10.mat');
end
%
aux_read=struct;
aux_read.opts_read=opts_read;
aux_align=struct;
aux_align.opts_align=opts_align;
aux_knit=struct;
aux_knit.opts_knit=opts_knit;
%
%do a dummy read to get file names
%
aux_read.nsets=0;
[data_in,aux_read_out]=rs_get_coordsets([],aux_read);
nfiles=length(aux_read_out.opts_read);
fullnames=cell(1,nfiles);
for ifile=1:nfiles
    fullnames{ifile}=aux_read_out.opts_read{ifile}.data_fullname;
end
%
fullnames=sort(fullnames); %sort file names in alphab order
shortnames=fullnames;
btc_names=cell(1,nfiles);
for ifile=1:nfiles
    maxsep=max(union(find(fullnames{ifile}=='\'),find(fullnames{ifile}=='/')));
    if ~isempty(maxsep)
        shortnames{ifile}=fullnames{ifile}(maxsep+1:end);
    end
    disp(sprintf(' file %2.0f: %s',ifile,shortnames{ifile}));
    btc_names{ifile}=shortnames{ifile}(1:4);
end
%
%specify subsets
%
subset_defs_used=cell(0);
for isubset=1:length(subset_defs)
    sd=subset_defs{isubset};
    if length(sd.contains)>0
        sub_contains=[];
        for k=1:length(sd.contains)
            sub_contains=[sub_contains,find(contains(btc_names,sd.contains{k}))];
        end
        sub_contains=unique(sub_contains);
    else
        sub_contains=[1:nfiles];
    end
    if ~isempty(sub_contains)
        subset_defs_used{end+1}=sd;
        subset_defs_used{end}.fileptrs=sub_contains;
        disp(sprintf('subset %2.0f',length(subset_defs_used)));
        disp(subset_defs_used{end})
    end
end
getinp('anything to proceed','d',[-Inf Inf],0);
nsubsets=length(subset_defs_used);
nprocs=nfiles+nsubsets;
%
aux_read.nsets=nfiles;
aux_read.opts_read.if_warn=0; %turn off warnings when augmenting by symmetry
aux_read.opts_read.if_symaug=1;
aux_read.opts_read.if_log=0;
%
data_in_aug=cell(nsymaugs,1);
aux_read_aug=cell(nsymaugs,1);
%
%
data_proc_align=cell(nsymaugs,nprocs);
data_proc_knit=cell(nsymaugs,nprocs);
aux_align_out=cell(nsymaugs,nprocs);
%aux_knit_out=cell(nsymaugs,nprocs); %this is very large, don't save each one separately
%
log_list=[];
%
%read, align and knit in two stages:
% first, within the files created from a primary file by augmenting by symmetry,
%    making data_sym_knit{isel}, one for each primary file that is combined
% then, knitting across the data_sym_knit files that descended from the primary files
%
%  If sym_apply='none', the first stage is trivial
%  If the subset is a singleton, the second stage is trivial
%
for isymaug=1:nsymaugs
    %
    sym_apply=sym_list{isymaug};
    disp('************');
    msg=sprintf('processing symmetry %1.0f: %s',isymaug,sym_apply);
    disp(msg);
    log_list=strvcat(log_list,msg);
    %
    %read with symmetry augmentation
    aux_read.opts_read.sym_apply=sym_apply;
    [data_in_aug{isymaug},aux_read_aug{isymaug}]=rs_get_coordsets(fullnames,aux_read);
    nfiles_aug=length(aux_read_aug{isymaug}.opts_read);
    data_filenames=cell(1,nfiles_aug);
    for ifile=1:nfiles_aug
        data_filenames{ifile}=aux_read_aug{isymaug}.opts_read{ifile}.data_fullname;
    end
    %add noise if requested
    if noise_jit>0
        for ifile=1:nfiles_aug
            data_jit=data_in_aug{isymaug}.ds{ifile};
            for k=1:length(data_jit)
                if ~isempty(data_jit{k})
                    nds=size(data_jit{k},2);
                    data_jit{k}=data_jit{k}+(noise_jit/sqrt(nds))*randn(size(data_jit{k},1),nds); %rms jitter across all dimensions
                end
            end
            data_in_aug{isymaug}.ds{ifile}=data_jit;
        end
    end
    sym_sets=cell(1,nfiles);
    ptrs_nosym=zeros(1,nfiles);
    disp(sprintf(' %3.0f original files, %3.0f after symmetry augmentation',nfiles,nfiles_aug));
    for isym_set=1:nfiles
        sym_sets{isym_set}=strmatch(fullnames{isym_set},data_filenames);
        disp(sprintf(' symmetry set for %2.0f (%30s): %2.0f augmented files, %3.0f to %3.0f',...
            isym_set,shortnames{isym_set},length(sym_sets{isym_set}),min(sym_sets{isym_set}),max(sym_sets{isym_set})));
    end
    %find the files needed for each computation
    for iproc=1:nprocs
        isubset=max(0,iproc-nfiles); %0 for a single file, >=1 for a subset
        if isubset==0
            select=iproc;
        else
            select=subset_defs_used{isubset}.fileptrs;
        end
        data_sym_in=struct;
        data_sym_knit=cell(length(select),1);
        %create data_sym_knit, knitting the symmetrized files for primary file select(isel)
        if_ok=1;
        for isel=1:length(select)
            select_aug=sym_sets{select(isel)}(:); %the augmented files
            data_sym.ds=data_in_aug{isymaug}.ds(select_aug);
            data_sym.sas=data_in_aug{isymaug}.sas(select_aug);
            data_sym.sets=data_in_aug{isymaug}.sets(select_aug);
            %align and knit
            [data_sym_align,aux_sym_align_out]=rs_align_coordsets(data_sym,aux_align);
            [data_sym_knit{isel},aux_sym_knit_out]=rs_knit_coordsets(data_sym_align,aux_knit);
            msg=sprintf(' iproc %2.0f, symmetry %10s, primary file %2.0f (%20s) has %2.0f files augmented by symmetry %8s; aligning and knitting ->',...
                iproc,sym_apply, select(isel),shortnames{select(isel)},length(select_aug),sym_apply);
            if ~isempty(data_sym_knit)
                log_msg=cat(2,msg,': success');
            else
                log_msg=cat(2,msg,': failure');
                if_ok=0;
            end
            disp(log_msg);
            log_list=strvcat(log_list,log_msg);
        end
        if if_ok==1
            msg=sprintf(' iproc %2.0f, symmetry %10s, combining',iproc, sym_apply);
            %format conversion
            data_knit=struct;
            for isel=1:length(select)
                data_knit.ds{isel}=data_sym_knit{isel}.ds{1};
                data_knit.sas{isel}=data_sym_knit{isel}.sas{1};
                data_knit.sets{isel}=data_sym_knit{isel}.sets{1};
            end
            %knit together across different primary files
            [data_proc_align{isymaug,iproc},aux_align_out{isymaug,iproc}]=rs_align_coordsets(data_knit,aux_align);
            if ~isempty(data_proc_align{isymaug,iproc})
                [data_proc_knit{isymaug,iproc},aux_knit_out]=rs_knit_coordsets(data_proc_align{isymaug,iproc},aux_knit);
                if ~isempty(data_proc_knit{isymaug,iproc})
                    log_msg=cat(2,msg,': success');
                else
                    log_msg=cat(2,msg,': knit failure');
                end
            else
                log_msg=cat(2,msg,': align failure');
            end
            disp(log_msg);
            log_list=strvcat(log_list,log_msg);
        end
    end
end
disp('*************');
%
%for each file, compute the participation ratio as a function of dimension
%
r=struct;
r.subj_id=subj_id;
r.d_max=d_max;
r.nfiles=nfiles;
r.nsubsets=nsubsets;
r.nprocs=nprocs;
r.nsymaugs=nsymaugs;
r.sym_list=sym_list;
r.noise_jit=noise_jit;
r.subset_defs_used=subset_defs_used;
r.log_list=log_list;
r.datasets=cell(nprocs,1);
r.datasets(1:nfiles)=shortnames;
r.btc_colorspec=btc_colorspec;
for isubset=1:nsubsets
    r.datasets{nfiles+isubset}=r.subset_defs_used{isubset}.name;
end
r.data_proc_knit=data_proc_knit;
%
%from here down, only quantities in r are used
%
%calculations
r.participation_ratio=zeros(r.nprocs,r.d_max,r.nsymaugs);
r.participation_ratio_dims={'d1: procs (dataset or subsets), d2: embedding dimension, d3: symmetry augmentation'};
for isymaug=1:r.nsymaugs
    sym_apply=r.sym_list{isymaug};
    for iproc=1:r.nprocs
        disp(sprintf('analyzing symmetry augmentation %10s of %s',sym_apply,r.datasets{iproc}));
        coord_sets=r.data_proc_knit{isymaug,iproc}.ds{1};
        for idim=1:r.d_max
            coords=coord_sets{idim};
            coords=coords-repmat(mean(coords,1),size(coords,1),1); %center
            [u,s,v]=svd(coords);
            pr=(sum(diag(s)).^2)/sum(diag(s).^2);
            r.participation_ratio(iproc,idim,isymaug)=pr;
        end
    end
end
%plots
figure;
set(gcf,'Position',[50 100 1300 800]);
set(gcf,'NumberTitle','off');
tstring=sprintf('subj: %3s, noise jitter  %7.4f (rms)',r.subj_id,r.noise_jit);
set(gcf,'Name',tstring);
[nrows,ncols]=nicesubp(r.nsymaugs,0.7);
for isymaug=1:r.nsymaugs
    sym_apply=r.sym_list{isymaug};
    subplot(nrows,ncols,isymaug);
    for iproc=1:r.nprocs
        hp=plot([1:r.d_max],r.participation_ratio(iproc,:,isymaug)','LineWidth',2);
        hold on;
        if iproc>r.nfiles
            isubset=iproc-r.nfiles;
            set(hp,'LineStyle','--');
            set(hp,'Color',r.subset_defs_used{isubset}.color);
        else
            btc_name=r.datasets{iproc}(1:4);
            if isfield(r.btc_colorspec,btc_name)
                set(hp,'Color',r.btc_colorspec.(btc_name));
            end
        end
    end
    set(gca,'XLim',[1 r.d_max]);
    set(gca,'XTick',[1:r.d_max]);
    set(gca,'YLim',[1 r.d_max]);
    set(gca,'YTick',[1:r.d_max]);
    xlabel('embedding dimension');
    ylabel('participation ratio');
    legend_text=strrep(strrep(r.datasets,'_coords',''),'.mat','');
    legend(legend_text,'Interpreter','none','Location','NorthWest','FontSize',7);
    title(sprintf('subj: %3s, symm aug: %8s',r.subj_id,sym_apply));
end
axes('Position',[0.01,0.02,0.01,0.01]); %for text
text(0,0,tstring,'Interpreter','none');
axis off;

%psg_knit_anal_demo: participation_ratio analysis of knitted-together experiments
% typically with four axes explored in each
%
%  See also:  RS_READ_COORDSETS, RS_KNIT_COORDSETS.
%
subj_id=upper(getinp('subject ID','s',[],'MC'));
if_symaug=getinp('1 to apply a symmetry','d',[0 1],1);
d_max=getinp('max dimension to analyze','d',[2 7],7);
%
opts_read.if_symaug=if_symaug;
if if_symaug
    sym_list=psg_btcsyms();
    disp(sym_list);
    if_match=0;
    while (if_match~=1)
        sym_apply=getinp('a symmetry to apply prior to knitting','s',[],'full');
        if_match=length(strmatch(sym_apply,sym_list,'exact'));
    end
else
    sym_apply='none';
end
%
opts_read=struct;
opts_read.if_gui=1;
opts_read.if_log=1;
opts_read.if_warn=1;
opts_read.if_auto=1;
opts_read.if_sym_aug=if_symaug;
opts_read.sym_apply=sym_apply;
%
opts_align=struct;
opts_align.if_log=1;
%
opts_knit=struct;
opts_knit.if_log=1;
opts_knit.allow_reflection=1;
opts_knit.allow_offset=1;
opts_knit.allow_scale=0;
opts_knit.if_normscale=1; %nondefault but irrelevant if allow_scale=0
opts_knit.if_pca=1; %nondefault
opts_knit.if_stats=0;
opts_knit.dim_max_in=d_max;
opts_read.ui_filter=cat(2,'./psg_data/*3pt_coords_',subj_id,'_sess01_10.mat');
%
aux_read=struct;
aux_read.opts_read=opts_read;
aux_align=struct;
aux_align.opts_align=opts_align;
aux_knit=struct;
aux_knit.opts_knit=opts_knit;
%
[data_in_aug,aux_read]=rs_get_coordsets([],aux_read);
%
%find range of each symmetry set
%
nfiles_aug=length(aux_read.opts_read);
data_filenames=cell(1,nfiles_aug);
for ifile=1:nfiles_aug
    data_fullname=aux_read.opts_read{ifile}.data_fullname;
    maxsep=max(union(find(data_fullname=='\'),find(data_fullname=='/')));
    if ~isempty(maxsep)
        data_fullname=data_fullname(maxsep+1:end);
    end
    data_filenames{ifile}=data_fullname;
end
%data_filename_list=unique(data_filenames,'stable');
data_filename_list=cellstr(unique(strvcat(data_filenames),'rows','stable')); %stable doesn't work with cell arrays
nfiles=length(data_filename_list);
sym_sets=cell(1,nfiles);
ptrs_nosym=zeros(1,nfiles);
for isym_set=1:nfiles
    sym_sets{isym_set}=strmatch(data_filename_list{isym_set},data_filenames);
    disp(sprintf(' symmetry set for file %2.0f (%30s): %2.0f elements, augmented files %3.0f to %3.0f',...
        isym_set,data_filename_list{isym_set},length(sym_sets{isym_set}),min(sym_sets{isym_set}),max(sym_sets{isym_set})));
    ptrs_nosym(isym_set)=min(sym_sets{isym_set});
end
%specify subsets
if_ok=0;
while (if_ok==0)
    nsubsets=getinp('number of subsets','d',[0 Inf]);
    subsets=cell(1,nsubsets);
    for isubset=1:nsubsets
        subsets{isubset}=getinp(sprintf('symmetry sets for subset %2.0f',isubset),'d',[1 nfiles]);
    end
    if_ok=getinp('1 if ok, and proceed to analysis','d',[0 1]);
end
%create a knit-together dataset from all original files, for each subset
data_in_nosym=cell(1,nsubsets);
data_align_nosym=cell(1,nsubsets);
data_knit_nosym=cell(1,nsubsets);
for isubset=1:nsubsets
    data_in_nosym{isubset}=struct;
    data_in_nosym{isubset}.ds=data_in_aug.ds(ptrs_nosym(subsets{isubset}));
    data_in_nosym{isubset}.sas=data_in_aug.sas(ptrs_nosym(subsets{isubset}));
    data_in_nosym{isubset}.sets=data_in_aug.sets(ptrs_nosym(subsets{isubset}));    
    [data_align_nosym{isubset},aux_align_nosym]=rs_align_coordsets(data_in_nosym{isubset},aux_align);
    [data_knit_nosym{isubset},aux_knit_nosym]=rs_knit_coordsets(data_align_nosym{isubset},aux_knit);
end
%
%for each file, compute the participation ratio as a function of dimension
%
r=struct;
r.subj_id=subj_id;
r.dataset=data_filename_list;
for isubset=1:nsubsets
    r.dataset{nfiles+isubset}=sprintf('subset %2.0f',isubset);
end
r.subsets=subsets;
r.participation_ratio=zeros(nfiles+1,d_max);
r.participation_ratio_dims={'d1: dataset, d2: dimension of space'};
for isym_set=1:nfiles+nsubsets
    disp(sprintf('analyzing %s',r.dataset{isym_set}));
    if isym_set<=nfiles
        coord_sets=data_in_aug.ds{ptrs_nosym(isym_set)};
    else
        coord_sets=data_knit_nosym{isym_set-nfiles}.ds{1};
    end
    for idim=1:d_max
        coords=coord_sets{idim};
        coords=coords-repmat(mean(coords,1),size(coords,1),1); %center
        [u,s,v]=svd(coords);
        pr=(sum(diag(s)).^2)/sum(diag(s).^2);
        r.participation_ratio(isym_set,idim)=pr;
    end
end

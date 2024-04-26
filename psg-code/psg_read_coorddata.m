function [d,sa,opts_used,pipeline]=psg_read_coorddata(data_fullname,setup_fullname,opts)
% [d,sa,opts_used,pipeline]=psg_read_coorddata(data_fullname,setup_fullname) reads
% coordinate data (typically from a multidimensional-scaling experiment),
% from mat-files transferred from Python fitting procedure
%
% data_fullname: full path and file name of data file, contains fields such as
%  dim1, dim2, ...; requested if not supplied or empty
% setup_fullname: full path and file name of a setup file, typically written by
%  psg_spokes_setup, with field s, and subfields s.nstims,s.specs,s.spec_labels, etc;  requested if not supplied
% opts.if_log: 1 to log, defaults to 0, can be omitted
% opts.if_justsetup: 1 to only read setup, defaults to 0
% opts.data_fullname_def: default data file
% opts.setup_fullname_def: default setup file
% opts.permutes: suggested ray permutations, e.g., permutes.bgca=[2 1 3 4]
% opts.permutes_ok: defaults to 1, to accept suggested ray permutations
% opts.domain_list: a cell array of names for domain experiments, defaults to
%    {'texture','intermediate_texture','intermediate_object','image','word'}
% opts.setup_suffix: suffix added to data fullname to create setup_fullname, defaults to '9' for compatibility with btc datasets
%
% d: d{k} has size [s.nstims k], and is the fits for the k-dimensional model
%   coordinates are re-ordered to match the order found in the setup, i.e., ordered by s.typenames.  
%   The original data file field stim_labels is checked to verify that all
%   stimuli in s.typenames are present.
% sa: selected fields of s, and also "stimulus coordinates", of two kinds. 
%   These are used for identifying rays and rings, determining nearest neighbors, etc.
%    btc_augcoords, augmented coords, of size [s.nstims nbtc] 
%    btc_specoords, specified coords, of size [s.nstims nbtc]
%  for btc, stimulus coords are taken from spec
%  for faces_mpi, stimulus coords are determined by psg_typenames2colors (age,gender,emo,set), but not individual
%  for other paradigms (materials, domains) stimulus coords are an identity matrix
% opts_used: options used
%    opts_used.dim_list: list of dimensions available
%    opts_used.data_fullname: data file full name used
%    opts_used.setup_fullname: setup file full name used
%    opts_used.type_class: 'btc', 'faces_mpi'
%  pipeline: pipeline field of data_fullname, if present
%
% 17Dec22: add logic for permute_raynums and allow a list; adapt setup file name to data file name
% 04Apr23: add opts.permutes_ok
% 01Jul23: modifications for compatibility with faces_mpi; add type_class; add face_prefix_list
% 03Jul23: add optional attenuations for each faces_mpi coord type and one-hot coords for emo and indiv
% 09Jul23: changes for irgb
% 27Sep23: add pipeline
% 08Nov23: changes for materials
% 28Nov23: if btc_augcoords or btc_specoords is present in s of metadata, it is copied to sa
% 22Feb24: localization params now from psg_localopts
% 23Feb24: documentation fixes
% 23Feb24: add domains experiment
% 26Apr24: make search for type_class start at first char, add type_class_aux
%
% See also: PSG_DEFOPTS, BTC_DEFINE, PSG_FINDRAYS, PSG_SPOKES_SETUP, BTC_AUGCOORDS, BTC_LETCODE2VEC,
%    PSG_VISUALIZE_DEMO, PSG_PLOTCOORDS, PSG_QFORMPRED_DEMO, PSG_TYPENAMES2COLORS, PSG_LOCALOPTS.
%
opts_local=psg_localopts;
xfr_fields={'nstims','nchecks','nsubsamp','specs','spec_labels','opts_psg','typenames','btc_dict','if_frozen_psg',...
    'spec_params','paradigm_name','paradigm_type',...
    'metadata','bsetModelLL','rawLLs','debiasedRelativeLL','biasEstimate'}; %these are for domains experiment
%
% face_prefix_list={'fc_','gy_'}; %prefixes on stimulus labels to be removed to match typenames (fc=full color, gy=gray)
%
dim_text='dim'; %leadin for fields of d
nrgb=3;
if (nargin<3)
    opts=struct;
end
if (nargin<1)
    data_fullname=[];
end
if (nargin<2)
    setup_fullname=[];
end
opts=filldefault(opts,'if_log',0);
opts=filldefault(opts,'if_justsetup',0);
opts=filldefault(opts,'data_fullname_def',opts_local.coord_data_fullname_def);
opts=filldefault(opts,'setup_fullname_def',opts_local.setup_fullname_def);
opts=filldefault(opts,'permutes_ok',1);
opts=filldefault(opts,'setup_suffix',opts_local.setup_setup_suffix);
%to allow for attenuation of any coordinate type in faces_mpi
opts=filldefault(opts,'faces_mpi_atten_indiv',1); %factor to attenuate "indiv" by in computing faces_mpi coords
opts=filldefault(opts,'faces_mpi_atten_age',1); %factor to attenuate "age" by in computing faces_mpi coords
opts=filldefault(opts,'faces_mpi_atten_gender',1); %factor to attenuate "gender" by in computing faces_mpi coords
opts=filldefault(opts,'faces_mpi_atten_emo',1); %factor to attenuate "emo" by in computing faces_mpi coords
opts=filldefault(opts,'faces_mpi_atten_set',0.2); %factor to attenuate "set" by in computing faces_mpi coords
%
opts=filldefault(opts,'domain_list',opts_local.domain_list_def); %domain list
if ~isfield(opts,'type_class_aux') 
    type_class_aux=[];
else
    type_class_aux=opts.type_class_aux;
end
%defaults to change plotting order
permutes=struct();
permutes.bgca=[2 1 3 4]; % permute ray numbers to the order gbca
% permutes.bdce=[1 3 2 4]; % permute ray numbers to the order bcde
opts=filldefault(opts,'permutes',permutes);
opts_used=opts;
pipeline=struct;
type_class=opts_local.type_class_def; %assumed type class
need_setup_file=1; %assume setup file is needed
coord_string='_coords';
if ~opts.if_justsetup
    if isempty(data_fullname)
        data_fullname=getinp('full path and file name of data file','s',[],opts.data_fullname_def);
    end
    underscore_sep=min(strfind(data_fullname,coord_string));
    if ~isempty(underscore_sep)
        %
        domain_match=0;
        for id=1:length(opts.domain_list)
            if contains(data_fullname,cat(2,opts.domain_list{id},'_coords'))
                domain_match=id;
            end
        end
        %
        data_shortname=data_fullname;
        seps=max(union(strfind(data_fullname,'/'),strfind(data_fullname,'\')));
        if ~isempty(seps)
            data_shortname=data_shortname(seps+1:end);
        end
        if ismember(1,strfind(data_shortname,'faces_mpi'))
            opts.setup_fullname_def=cat(2,data_fullname(1:underscore_sep-1),'.mat');
            type_class='faces_mpi';
        elseif ismember(1,strfind(data_shortname,'irgb'))
            opts.setup_fullname_def=cat(2,data_fullname(1:underscore_sep-1),'.mat');
            type_class='irgb';
        elseif ismember(1,strfind(data_shortname,'mater'))
            opts.setup_fullname_def=cat(2,data_fullname(1:underscore_sep-1),'.mat');
            type_class='mater';
        elseif ismember(1,strfind(data_shortname,type_class_aux))
            opts.setup_fullname_def=data_fullname;
            type_class=type_class_aux;
            need_setup_file=0; %no setup file needed
        elseif domain_match>0
            opts.setup_fullname_def='';
            type_class='domain';
            need_setup_file=0; %no setup file needed
            %find subject ID
            subjid_string=data_fullname(underscore_sep+length(coord_string)+1:end); %should be something like EVF.mat
            subjid_string=strrep(subjid_string,'.mat',''); %remove extension if present
            if isfield(opts_local.domain_sigma,subjid_string)
                domain_sigma=opts_local.domain_sigma.(subjid_string);
            else
                domain_sigma=1;
                warning(sprintf('subject ID %2.0s from file name %s not recognized',data_fullname,subjid_string));
            end
        else
            opts.setup_fullname_def=cat(2,data_fullname(1:underscore_sep-1),opts.setup_suffix,'.mat');
        end
    end
else
    data_fullname=[];
end
if need_setup_file
    if isempty(setup_fullname)
        setup_fullname=getinp('full path and file name of psg setup file','s',[],opts.setup_fullname_def);
    end
else
    setup_fullname='[unused]';
end
if ~opts.if_justsetup
    d_read=load(data_fullname);
    d_fields=fieldnames(d_read);
    if isfield(d_read,'pipeline')
        pipeline=d_read.pipeline;
    end
    if (opts.if_log)
        disp(sprintf('%3.0f different stimulus types found in  data file %s',length(d_read.stim_labels),data_fullname));
    end
end
% permute the ray order?
opts.permute_raynums=[];
if isstruct(opts.permutes)
    perm_list=fieldnames(opts.permutes);
    for iperm=1:length(perm_list)
        if strfind(setup_fullname,perm_list{iperm})
            disp(sprintf('suggested ray permutation for %s:',perm_list{iperm}))
            disp(opts.permutes.(perm_list{iperm}));
            if ~opts.permutes_ok
                if getinp('1 if ok','d',[0 1],1)
                    opts.permute_raynums=permutes.(perm_list{iperm});
                end
            else
                opts.permute_raynums=permutes.(perm_list{iperm});
            end
        end
    end
end
%
opts.data_fullname=data_fullname;
opts.setup_fullname=setup_fullname;
opts_used=opts;
opts_used.type_class=type_class;
%
%read the setup; retrieve specified and augmented coordinates
%
if need_setup_file %setup file needed for metadata
    s=getfield(load(setup_fullname),'s');
    if strcmp(type_class,'mater') %install nstims and adjoin setup file name to typename
        s.nstims=size(s.typenames,1)
        setup_basename=strrep(setup_fullname,'.mat','');
        setup_basename=setup_basename(max(strfind(cat(2,'/',setup_basename),'/')):end);
        setup_basename=setup_basename(max(strfind(cat(2,'\',setup_basename),'\')):end);
        for istim=1:s.nstims
             s.typenames{istim}=cat(2,setup_basename,'-',s.typenames{istim});   
        end
    end
else %create metadata without a file
    s=d_read; %capture metadata and likelihoods this way
    s.nstims=size(d_read.stim_labels,1);
    s.typenames=cell(s.nstims,1);
    for istim=1:s.nstims
        s.typenames{istim}=deblank(d_read.stim_labels(istim,:));
    end 
end
%
if (opts.if_log)
    disp(sprintf('%3.0f different stimulus types found in setup file %s',s.nstims,setup_fullname));
end
sa=struct;
for ifield=1:length(xfr_fields)
    if isfield(s,xfr_fields{ifield})
        sa.(xfr_fields{ifield})=s.(xfr_fields{ifield});
    end
end
%make sure that sa has btc_specoords, and, for type_class='btc',also btc_augcoords
%files written at time of data collection won't have these fields
%files written in a processing pipeline might have these fields
switch type_class
    case 'btc'
        nbtc=length(s.btc_dict.codel);
        if isfield(s,'btc_specoords') & isfield(s,'btc_augcoords')
            if size(s.btc_specoords,2)==nbtc & size(s.btc_augcoords,2)==nbtc
                sa.btc_specoords=s.btc_specoords;
                sa.btc_augcoords=s.btc_augcoords;
            else
                error('btc_specooords or btc_augcoords cannot be retrieved from metadata file.');
            end
        else
            sa.btc_augcoords=zeros(s.nstims,nbtc);
            sa.btc_specoords=zeros(s.nstims,nbtc);
            for istim=1:s.nstims
                sa.btc_augcoords(istim,:)=s.btc_methods{istim}.vec;
                sa.btc_specoords(istim,:)=btc_letcode2vec(s.specs{istim},s.btc_dict);
            end
        end
    case 'faces_mpi'
        if isfield(s,'btc_specoords')
            sa.btc_specoords=s.btc_specoords;
        else
            %use coordinates extracted from type names, as modified by opts.faces_mpi_atten_[indiv|age|gender|emo|set]
            [rgb,symb,vecs,ou]=psg_typenames2colors(s.typenames,[]);
            all_coords=ou.faces_mpi.attrib_table_num;
            %order in attrib_table_num is given by attrib_table_order, after omitting indiv: {'indiv'  'age'  'gender'  'emo'  'set'}       
            indiv_col=strmatch('indiv',ou.faces_mpi.attrib_table_order);
            age_col=strmatch('age',ou.faces_mpi.attrib_table_order);
            gender_col=strmatch('gender',ou.faces_mpi.attrib_table_order);
            emo_col=strmatch('emo',ou.faces_mpi.attrib_table_order);
            set_col=strmatch('emo',ou.faces_mpi.attrib_table_order);
            %coords based on individuals: one-hot actually present
            indiv_unique=unique(all_coords(:,indiv_col));
            indiv_coords=zeros(s.nstims,length(indiv_unique)); %one-hot for individual coords
            for iu=1:length(indiv_unique)
                indiv_coords(find(all_coords(:,indiv_col)==indiv_unique(iu)),iu)=opts.faces_mpi_atten_indiv;
            end
            %coords based on emo: one-hot, based on total available list
            nemo=ou.faces_mpi.attrib_info.emo.nlevels;
            emo_coords=zeros(s.nstims,nemo);
            for iu=1:nemo
                emo_coords(find(all_coords(:,emo_col)==iu),iu)=opts.faces_mpi_atten_emo;
            end
            %coords based on age
            age_coords=opts.faces_mpi_atten_age*all_coords(:,age_col);
            %coords based on gender
            gender_coords=opts.faces_mpi_atten_gender*all_coords(:,gender_col);
            %coords based on set
            set_coords=opts.faces_mpi_atten_set*all_coords(:,set_col);
            sa.btc_specoords=[age_coords gender_coords set_coords emo_coords indiv_coords];
        end
        %
        %truncate suffixes on stim_labels if necessary
        %
        for ifp=1:length(face_prefix_list)
            d_read.stim_labels=char(strrep(cellstr(d_read.stim_labels),face_prefix_list{ifp},''));
        end
    case 'irgb'
        if isfield(s,'btc_specoords')
            sa.btc_specoords=s.btc_specoords;
        else
            %
            %use mean and diagonal of covariance as coordinates but if either
            %means are all identical or covs are all identical, then set to
            %zero so that rays can be found
            %
            irgb_coords=zeros(s.nstims,nrgb*2);
            for istim=1:s.nstims
                irgb_coords(istim,:)=[s.specs{istim}.mean_val,diag(s.specs{istim}.cov)'];
            end
            if all(min(irgb_coords(:,1:nrgb))==max(irgb_coords(:,1:nrgb)))
                irgb_coords(:,1:nrgb)=0;
            elseif  all(min(irgb_coords(:,nrgb+[1:nrgb]))==max(irgb_coords(:,nrgb+[1:nrgb])))
                irgb_coords(:,nrgb+[1:nrgb])=0;
            end
            sa.btc_specoords=irgb_coords;
        end
    case 'mater'
        if isfield(s,'btc_specoords')
            sa.btc_specoords=s.btc_specoords;
        else
            sa.btc_specoords=eye(s.nstims); %no meaningful coordinates
        end
    case type_class_aux
        sa.btc_specoords=eye(s.nstims);
    case 'domain'
        sa.btc_specoords=eye(s.nstims); %no meaningful coords
        sa.paradigm_name=opts.domain_list{domain_match};
        sa.sigma_orig=domain_sigma;
        sa.sigma_info='coordinates have been normalized to sigma=1';
end
%parse the data
d=cell(0);
if ~opts.if_justsetup
    %d=d_read;
    d_fields=fieldnames(d_read);
    stim_labels=d_read.stim_labels;
    stim_sorted=zeros(s.nstims,1);
    for istim=1:s.nstims
        stim_match=strmatch(s.typenames{istim},stim_labels,'exact');
        if length(stim_match)~=1
            stim_match=0;
        end
        stim_sorted(istim)=stim_match;
    end
    stim_found=find(stim_sorted>0);
    if (opts.if_log) | any(stim_sorted==0)
        disp(sprintf('%3.0f of %3.0f labels found',sum(stim_sorted>0),s.nstims));
    end
    if any(stim_sorted==0)
       warning(sprintf('%3.0f stimulus labels not found',sum(stim_sorted==0)));
    end
    nmods=length(d_fields);
    dim_list=[];
    for ifield=1:length(d_fields)
        fn=d_fields{ifield};
        ind=strfind(fn,dim_text);
        dimno=str2num(fn(ind+length(dim_text):end));
        if dimno>0
            coords=NaN(s.nstims,dimno);
            coords_read=d_read.(fn);
            coords(stim_found,:)=coords_read(stim_sorted(stim_found),:); %unscramble
            if strcmp(type_class,'domain')
                coords=coords/domain_sigma;
            end
            d{dimno}=coords;
            dim_list=sort([dim_list,dimno]);
        end
    end
    opts_used.dim_list=dim_list;
    if (opts.if_log)
        disp(sprintf('models with %2.0f to %2.0f dimensions read.',min(dim_list),max(dim_list)));
    end
end
return
end

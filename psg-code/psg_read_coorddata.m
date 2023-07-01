function [d,sa,opts_used]=psg_read_coorddata(data_fullname,setup_fullname,opts)
% [d,sa,opts_used]=psg_read_coorddata(data_fullname,setup_fullname) reads
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
%
% d: d{k} has size [s.nstims k], and is the fits for the k-dimensional model
%   coordinates are re-ordered to match the order found in the setup
% sa: selected fields of s, and also
%    btc_augcoords, augmented coords, of size [s.nstims nbtc] 
%    btc_specoords, specificed coords, of size [s.nstims nbtc]
%
% for btc, stimulus coords are taken from spec
% for faces_mpi, stimulus coords are determined by psg_typenames2colors (age,gender,emo,set), but not individual
%
% opts_used: options used
%    opts_used.dim_list: list of dimensions available
%    opts_used.data_fullname: data file full name used
%    opts_used.setup_fullname: setup file full name used
%    opts_used.type_class: 'btc', 'faces_mpi'
%
%
% 17Dec22: add logic for permute_raynums and allow a list; adapt setup file name to data file name
% 04Apr23: add opts.permutes_ok
% 01Jul23: modifications for compatibility with faces_mpi; add type_class; 
%
% See also: PSG_DEFOPTS, BTC_DEFINE, PSG_FINDRAYS, PSG_SPOKES_SETUP, BTC_AUGCOORDS, BTC_LETCODE2VEC,
%    PSG_VISUALIZE_DEMO, PSG_PLOTCOORDS, PSG_QFORMPRED_DEMO.
%
xfr_fields={'nstims','nchecks','nsubsamp','specs','spec_labels','opts_psg','typenames','btc_dict','if_frozen_psg'};
dim_text='dim'; %leadin for fields of d
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
opts=filldefault(opts,'data_fullname_def','./psg_data/bgca3pt_coords_BL_sess01_04.mat');
opts=filldefault(opts,'setup_fullname_def','./psg_data/bgca3pt9.mat');
opts=filldefault(opts,'permutes_ok',1);
%
%defaults to change plotting order
permutes=struct();
permutes.bgca=[2 1 3 4]; % permute ray numbers to the order gbca
% permutes.bdce=[1 3 2 4]; % permute ray numbers to the order bcde
opts=filldefault(opts,'permutes',permutes);
opts_used=opts;
type_class='btc'; %assume btc
if ~opts.if_justsetup
    if isempty(data_fullname)
        data_fullname=getinp('full path and file name of data file','s',[],opts.data_fullname_def);
    end
    underscore_sep=min(strfind(data_fullname,'_coords'));
    if ~isempty(underscore_sep)
        %
        if ~isempty(strfind(data_fullname,'faces_mpi'))
            opts.setup_fullname_def=cat(2,data_fullname(1:underscore_sep-1),'.mat');
            type_class='faces_mpi';
        else
            opts.setup_fullname_def=cat(2,data_fullname(1:underscore_sep-1),'9.mat');
        end
    end
else
    data_fullname=[];
end
if isempty(setup_fullname)
    setup_fullname=getinp('full path and file name of psg setup file','s',[],opts.setup_fullname_def);
end
if ~opts.if_justsetup
    d_read=load(data_fullname);
    d_fields=fieldnames(d_read);
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
s=getfield(load(setup_fullname),'s');
if (opts.if_log)
    disp(sprintf('%3.0f different stimulus types found in setup file %s',s.nstims,setup_fullname));
end
sa=struct;
for ifield=1:length(xfr_fields)
    if isfield(s,xfr_fields{ifield})
        sa.(xfr_fields{ifield})=s.(xfr_fields{ifield});
    end
end
switch type_class
    case 'btc'
        nbtc=length(s.btc_dict.codel);
        sa.btc_augcoords=zeros(s.nstims,nbtc);
        sa.btc_specoords=zeros(s.nstims,nbtc);
        for istim=1:s.nstims
            sa.btc_augcoords(istim,:)=s.btc_methods{istim}.vec;
            sa.btc_specoords(istim,:)=btc_letcode2vec(s.specs{istim},s.btc_dict);
        end
    case 'faces_mpi'
        %use coordinates extracted from type names
        %order is given by attrib_table_order, after omitting indiv: {'indiv'  'age'  'gender'  'emo'  'set'}       
        [rgb,symb,vecs,ou]=psg_typenames2colors(s.typenames,[]);
        sa.btc_specoords=ou.faces_mpi.attrib_table_num(:,2:end); % omit 'indiv'
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

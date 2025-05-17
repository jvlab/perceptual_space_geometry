function [model,if_warn,opts_used]=btc_soidfg_model(opts,btc_dict)
% [model,if_warn,opts_used]=btc_soidfg_model(opts,btc_dict) defines a model for figure-ground psychophysical data
%
% case by case to set up tables for each symmetry type -- a bit clunky
%
% input:
%  opts: options structure from btc_soidfg_define, [] if not supplied; btc_soidfg_define called to fill with defaults
%   uses the following fields:
%      coords:  the coordinate set to include in the model (e.g.,'gbcdetuvwa' or 'bcde')
%      sym_type: one of opts.sym_type (see btc_soidfg_define)
%      model_type: one of opts.model_type (see btc_soidfg_define)
%      model_exclude: model parameters to exclude by name, defaults to cell(0), an entry of e.g., '_d_e' to exclude (d,e) cross-terms
%      model_axes_probed: cell array  (1,n) of axes probed, defaults to cell(0); an entry {'b','c','d'} means that 
%           a self-model term is kept only if it (or a symmetry-equivalent) uses one of these axes
%      model_planes_probed: cell array (1,n) of planes probed, defaults to cell(0); an entry {'bc','bd','de'} means that
%           a cross-model term is kept only if it (or a symmetry-equivalent) uses one of these planes
%  btc_dict: structure returned by btc_define; can be omitted
%
% output:
%  model: structure defining the model
%    model.nparams: number of parameters
%    model.qforms:  array of size (20,20,nparams) indicating the
%        contribution of each parameter to a quadratic threshold model.
%        Each slice of model.qforms is symmetric, and element (j,k)
%        indicates contribution to (btc_j*btc_k), where btc_j and btc_k are
%        btc coordinates, first 10 in figure, second 10 in ground
%    model.param_name:  cell array of size (nparams,1) with parameter names
%    model.param_type: string label for parameter type
%        'self_fmg' (self-term, figure minus ground)
%        'cross_fmg' (cross-term, figure minus ground)
%        'self_fig' (self-term, figure)
%        'cross_fig (cross-term, figure)
%        'self_gnd' (self-term, ground)
%        'cross_gnd' (cross-term, ground)
%        'self_fig_gnd' (self-term, fig and ground linked)
%        'cross_fig_gnd' (cross-term, fig and gnd linked)
%    model.param_name_equiv: cell array of size (nparams,1) containing parameter names that are equivalent up to symmetry
%      below, list_name={'single','pair_ordered','pair_sorted'}
%    model.equiv_tables.(list_name): entries are a cell array of all equivalent items to equiv_tables.(list_name).(entry)
%    model.unique_tables.(list_name).(uequiv): unique elements up to symmetry, value is number of items that unique_tables.(list_name).(uequiv) maps to
%    model.counts.(list_name): number of entries in model.equiv_tables.(list_name)
%  if_warn: 1 if symmetry classes are not contained in coords
%  opts_used: options after defaults filled in
%  
% 27Apr20:  Add options to exclude based on axes and planes probed
% 10Jun20:  Pass btc_dict to btc_hflip, btc_vflip, btc_rotcode, btc_hvi, btc_hvirev to speed up
% 12Jun20:  Change parameter types in linked fig-gnd model from 'self_fig+self_gnd' to 'self_fig_plus_self_gnd' and same for cross
% 24Jun20:  Bug fix: change logic for exclusion based on axis or plane so that for allquad, fig_gnd_b_b, etc, is handled as an on-axis case
%      affects f-g,f,g(linked) and allquad
%
%  See also:  BTC_SOIDFG_DEFINE, BTC_DEFINE, BTC_SOIDFG_MODEL_TEST, BTC_SYMCLASSES,
%  BTC_HFLIP, BTC_VFLIP, BTC_HVI, BTC_HVIREV, BTC_ROTCODE, BTC_SOIDFG_MODEL_PROJ, BTC_SOIDFG_NAME2LETS.
%
if (nargin==0)
    opts=[];
end
if (nargin<=1)
    btc_dict=btc_define([]);
end
opts=btc_soidfg_define(opts);
opts_used=opts;
codel=btc_dict.codel;
btc_n=length(codel); %10
%
model=struct();
model.qforms=zeros(2*btc_n,2*btc_n,0);
model.param_name=cell(0);
model.param_type=cell(0);
model.param_name_equiv=cell(0);
[model.symclasses,if_warn]=btc_symclasses(opts,btc_dict);
%
btc_pos=nan(length(opts.coords),1);
btc_nvalid=0;
lists.single=cell(0);
for ilet=1:length(opts.coords)
    let=opts.coords(ilet);
    letpos=find(codel==let);
    if length(letpos)~=1
        warning(sprintf('opts.coords, %s, contains an unrecognized coordinate; coords must be in %s',opts.coords,codel));
    else
        btc_pos(ilet)=letpos;
        btc_nvalid=btc_nvalid+1;
        lists.single{btc_nvalid,1}=let;
    end
end
model.btcpos=btc_pos;
%
%create cell arrays of all single letters, all pairs, all ordered pairs
%
lists.pair_sorted=cell(btc_nvalid*(btc_nvalid-1)/2,1);
lists.pair_ordered=cell(btc_nvalid*(btc_nvalid-1),1);
pair_sorted=0;
pair_ordered=0;
for ilet=1:btc_nvalid
    for jlet=1:btc_nvalid
        ijlet=cat(2,lists.single{ilet},lists.single{jlet});
        pair_ordered=pair_ordered+1;
        lists.pair_ordered{pair_ordered}=ijlet;
        if jlet>ilet
            pair_sorted=pair_sorted+1;
            lists.pair_sorted{pair_sorted}=sort(ijlet);
        end
    end
end
model.lists=lists;
%
%create tables that give the symmetry-equivalents of the lists for the three fields of model.lists: 'single','pair_sorted','pair_ordered'
%  equiv_tables.(list_name): entries are a cell array of all equivalent items to equiv_tables.(list_name).(entry)
%  unique_tables.(list_name).(uequiv): ueqiv are elements up to symmetry, entries are the number of items that map to it
%  * to find the equivalent value of any item, use equiv_tables.(list_name).(entry){1}
equiv_tables=struct();
unique_tables=struct();
find_unique=struct();
list_fields=fieldnames(lists);
for ifield=1:length(list_fields)
    list_name=list_fields{ifield};
    equiv_tables.(list_name)=struct;
    list=lists.(list_name);
    needsort=strcmp(list_name,'pair_sorted');
    for ilist=1:length(list)
        entry=list{ilist};
        equiv_tables.(list_name).(entry){1}=entry;
        switch opts.sym_type
            case 'none'
            case {'hflip','vflip','rot180'}
                switch opts.sym_type
                    case 'hflip'
                        trans=btc_hflip(entry,btc_dict);
                    case 'vflip'
                        trans=btc_vflip(entry,btc_dict);
                    case 'rot180'
                        trans=btc_rotcode(entry,2,btc_dict);
                end
                if (needsort)
                    trans=sort(trans);
                end
                equiv_tables.(list_name).(entry){2}=trans;
                equiv_tables.(list_name).(entry)=unique(equiv_tables.(list_name).(entry));
            case 'rot90'
                rots=cell(1,4);
                rots{1}=entry;
                for irot=1:3
                    trans=btc_rotcode(entry,irot,btc_dict);
                    if (needsort)
                        trans=sort(trans);
                    end
                    rots{irot+1}=trans;
                end
                equiv_tables.(list_name).(entry)=unique(rots);
            case {'hvflip'} 
                mods=cell(1,4);
                mods{1}=entry;
                mods{2}=btc_hflip(entry,btc_dict);
                mods{3}=btc_vflip(entry,btc_dict);
                mods{4}=btc_rotcode(entry,2,btc_dict);
                if (needsort)
                    for imod=1:length(mods)
                        mods{imod}=sort(mods{imod});
                    end
                end
                equiv_tables.(list_name).(entry)=unique(mods);
            case {'diag','full'} %do the full cases by recognizing that these are the diagonal-flip cases, optionally composed with hflip
                diagmods=cell(1,4);
                diagmods{1}=entry;
                diagmods{2}=btc_hvi(entry,btc_dict);
                diagmods{3}=btc_hvirev(entry,btc_dict);
                diagmods{4}=btc_rotcode(entry,2,btc_dict);
                if strcmp(opts.sym_type,'full')
                    mods=cell(1,8);
                    for imod=1:4
                        mods{imod}=diagmods{imod};
                        mods{4+imod}=btc_hflip(diagmods{imod},btc_dict);
                    end
                else
                    mods=diagmods;
                end
                if (needsort)
                    for imod=1:length(mods)
                        mods{imod}=sort(mods{imod});
                    end
                end
                equiv_tables.(list_name).(entry)=unique(mods);
        end %sym_type
    end %list 
    %create table of unique-up-to-symmetry items, and how often each is used
    unique_tables.(list_name)=struct;
    entry_names=fieldnames(equiv_tables.(list_name));
    for ientry=1:length(entry_names)
        entry=entry_names{ientry};
        equivs=equiv_tables.(list_name).(entry);
        uequiv=equivs{1};
        if isfield(unique_tables.(list_name),uequiv)
            unique_tables.(list_name).(uequiv)=unique_tables.(list_name).(uequiv)+1;
        else
            unique_tables.(list_name).(uequiv)=1;
        end
        find_unique.(list_name).(entry)=uequiv;
    end
    model.counts.(list_name)=length(entry_names);
end %list
model.equiv_tables=equiv_tables;
model.unique_tables=unique_tables;
%
% now use the tables to set up model
%
%
%terms involving a single coord, in figure or ground
%
unames=fieldnames(model.unique_tables.single);
for iz=1:length(unames)
    uname=unames{iz};
    equivs=model.equiv_tables.single.(uname);
    if ~isempty(strfind(opts.model_type,'f-g')) %figure minus ground self-terms
        model=btc_model_self(model,equivs,[1 -1],'fmg_','self_fmg',btc_n,codel);
    end
    if strcmp(opts.model_type,'f-g,f') | strcmp(opts.model_type,'f-g,f,g(indep)') | strcmp(opts.model_type,'allquad') %figure self-terms
        model=btc_model_self(model,equivs,[1 0],'fig_','self_fig',btc_n,codel);
    end
    if strcmp(opts.model_type,'f-g,g') | strcmp(opts.model_type,'f-g,f,g(indep)') | strcmp(opts.model_type,'allquad') %ground self-terms
        model=btc_model_self(model,equivs,[0 1],'gnd_','self_gnd',btc_n,codel);
    end
    if strcmp(opts.model_type,'f-g,f,g(linked)')
        model=btc_model_self(model,equivs,[1 0;0 1],'fig_gnd_','self_fig_plus_self_gnd',btc_n,codel); %figure-self and linked ground-self
    end
end %iz
%
%cross-terms using sorted pairs 
%
unames=fieldnames(model.unique_tables.pair_sorted);
for iz=1:length(unames)
    uname=unames{iz};
    equivs=model.equiv_tables.pair_sorted.(uname);
    if ~isempty(strfind(opts.model_type,'f-g')) %figure minus ground cross-terms
        model=btc_model_cross(model,equivs,[1 -1],[1 -1],'fmg_','cross_fmg',btc_n,codel);
    end
    if strcmp(opts.model_type,'f-g,f') | strcmp(opts.model_type,'f-g,f,g(indep)') | strcmp(opts.model_type,'allquad') %figure cross-terms
        model=btc_model_cross(model,equivs,[1 0],[1 0],'fig_','cross_fig',btc_n,codel);
    end
    if strcmp(opts.model_type,'f-g,g') | strcmp(opts.model_type,'f-g,f,g(indep)') | strcmp(opts.model_type,'allquad') %ground cross-terms
        model=btc_model_cross(model,equivs,[0 1],[0 1],'gnd_','cross_gnd',btc_n,codel);
    end
    if strcmp(opts.model_type,'f-g,f,g(linked)')
        model=btc_model_cross(model,equivs,[1 0;0 1],[1 0;0 1],'fig_gnd_','cross_fig_plus_cross_gnd',btc_n,codel); %figure-cross and linked ground-cross
    end
end
%
%cross-terms using ordered pairs: first element of pair is figure, second is ground
%
unames=fieldnames(model.unique_tables.pair_ordered);
for iz=1:length(unames)
    uname=unames{iz};
    equivs=model.equiv_tables.pair_ordered.(uname);
    if strcmp(opts.model_type,'allquad')
        model=btc_model_cross(model,equivs,[1 0],[0 1],'fig_gnd_','cross_fig_gnd',btc_n,codel);
    end
end
%
%exclude parameters based on 
%  explicit list of names
%  whether axes are used
%  whether planes are used
%
keep=struct();
keep.qforms=zeros(2*btc_n,2*btc_n,0);
keep.param_name=cell(0);
keep.param_type=cell(0);
keep.param_name_equiv=cell(0);
ikeep=0;
exclude_name=0;
exclude_axes=0;
exclude_planes=0;
for iparam=1:size(model.qforms,3)
    if_keep_name=1;
    if_keep_axes=1;
    if_keep_planes=1;
    if ~isempty(opts.model_exclude) %exclude on basis of param name?
        for iex=1:length(opts.model_exclude)
            if_keep_name=and(if_keep_name,isempty(strfind(model.param_name{iparam},opts.model_exclude{iex})));
        end
        exclude_name=exclude_name+double(if_keep_name==0);
    end
    %to keep a self-parameter (single btc coord), then one of its symmetry-equivalents must be in opts.model_axes_probed
    %if ~isempty(opts.model_axes_probed) & ~isempty(strfind(model.param_type{iparam},'self')) 
    if ~isempty(opts.model_axes_probed) & length(unique(btc_soidfg_name2lets(model.param_name{iparam})))==1  %mod 24Jun20
        if_keep_axes=0; %assume not ok
        iequiv=0;
        while (if_keep_axes==0) & iequiv<length(model.param_name_equiv{iparam})
            iequiv=iequiv+1;
            eqname=model.param_name_equiv{iparam}{iequiv};
            %btc_lets=eqname(find(eqname=='_')+1); %extract the btc coords used in the equivalent name
            btc_lets=unique(btc_soidfg_name2lets(eqname)); %extract the btc coords used in the equivalent name  %mod 24Jun20
            if all(ismember(btc_lets,char(opts.model_axes_probed))) %are all of the coords listed in model_axes_probed?
                if_keep_axes=1;
            end
        end
        if (if_keep_axes==0)
            exclude_axes=exclude_axes+1;
        end
    end
    %to keep a cross-parameter (two btc coords), then one of its symmetry-equivalents must use a plane that is in opts.model_planes_probed
    %if ~isempty(opts.model_planes_probed) & ~isempty(strfind(model.param_type{iparam},'cross')) 
    if ~isempty(opts.model_planes_probed) & length(unique(btc_soidfg_name2lets(model.param_name{iparam})))==2 %mod 24Jun20
        if_keep_planes=0; %assume not ok
        iequiv=0;
        while (if_keep_planes==0) & iequiv<length(model.param_name_equiv{iparam})
            iequiv=iequiv+1;
            eqname=model.param_name_equiv{iparam}{iequiv};
            %btc_lets=eqname(find(eqname=='_')+1); %extract the btc coords used in the equivalent name
            btc_lets=btc_soidfg_name2lets(eqname); %extract the btc coords used in the equivalent name  %mod 24Jun20
            for iplane=1:length(opts.model_planes_probed)
                if (all(ismember(btc_lets,opts.model_planes_probed{iplane})))
                    if_keep_planes=1;
                end
            end
        end
        if (if_keep_planes==0)
            exclude_planes=exclude_planes+1;
        end

    end
    if_keep=if_keep_name & if_keep_axes & if_keep_planes;
    if (if_keep)
        ikeep=ikeep+1;
        keep.qforms(:,:,ikeep)=model.qforms(:,:,iparam);
        keep.param_name{ikeep}=model.param_name{iparam};
        keep.param_type{ikeep}=model.param_type{iparam};
        keep.param_name_equiv{ikeep}=model.param_name_equiv{iparam};
    end
end
if (opts.verbose)
    disp(sprintf(' model has %4.0f params, %4.0f excluded by param name, %4.0f excluded by axis, %4.0f excluded by plane, %4.0f kept',...
        iparam,exclude_name,exclude_axes,exclude_planes,ikeep));
end

model.qforms=keep.qforms;
model.param_name=keep.param_name;
model.param_type=keep.param_type;
model.param_name_equiv=keep.param_name_equiv;
model.nparams=size(model.qforms,3);
return

function model_new=btc_model_self(model,equivs,signs,name_prefix,param_type,btc_n,codel);
%create entries in the model structure for a self term
% model: model structure to be updated
% equivs: equivalent coords under allowed symmetries; each entry is a btc letter
% signs: one or more rows, each row is a two-vector of 0,+1,-1 indicating contribution of figure and ground
%    two rows allow for linking figure and ground
% name_prefix: prefix for parameter name
% param_type: entry for model.param_type
% btc_n: length(codel)
% codel: standard string of btc coords, from btc_define
%
q=zeros(btc_n*2,btc_n*2);
names_equiv=cell(0);
for iequiv=1:length(equivs)
    uname_equiv=equivs{iequiv};
    names_equiv{1,end+1}=cat(2,name_prefix,uname_equiv);
    hvec=zeros(1,btc_n);
    hvec(codel==uname_equiv)=1;
    for ilink=1:size(signs,1)
        qvec_x=[hvec*signs(ilink,1) hvec*signs(ilink,2)];
        q=q+qvec_x'*qvec_x; %outer product
    end %ilink
end
model.qforms(:,:,end+1)=q;
model.param_name{1,end+1}=cat(2,name_prefix,equivs{1});
model.param_type{1,end+1}=param_type;
model.param_name_equiv{1,end+1}=names_equiv;
%
model_new=model;
%
return

function model_new=btc_model_cross(model,equivs,signs_x,signs_y,name_prefix,param_type,btc_n,codel);
%create entries in the model structure for a cross-term
% model: model structure to be updated
% equivs: equivalent coords under allowed symmetries; each entry is a btc letter pair. 
% signs_x, signs_y: one or more rows, each row is a two-vector of 0,+1,-1 indicating contribution of figure and ground
%    two rows allow for linking figure and ground
% name_prefix: prefix for parameter name
% param_type: entry for model.param_type
% btc_n: length(codel)
% codel: standard string of btc coords, from btc_define
%
q=zeros(btc_n*2,btc_n*2);
names_equiv=cell(0);; 
for iequiv=1:length(equivs)
    uname_equiv=equivs{iequiv};
    names_equiv{1,end+1}=cat(2,name_prefix,uname_equiv(1),'_',uname_equiv(2));
    hvec_x=zeros(1,btc_n);
    hvec_x(codel==uname_equiv(1))=1;
    hvec_y=zeros(1,btc_n);
    hvec_y(codel==uname_equiv(2))=1;
    for ilink=1:size(signs_x,1)
        qvec_x=[hvec_x*signs_x(ilink,1) hvec_x*signs_x(ilink,2)];
        qvec_y=[hvec_y*signs_y(ilink,1) hvec_y*signs_y(ilink,2)];
        q=q+qvec_x'*qvec_y+qvec_y'*qvec_x; %outer product
    end %ilink
end
model.qforms(:,:,end+1)=q;
model.param_name{1,end+1}=cat(2,name_prefix,equivs{1}(1),'_',equivs{1}(2));
model.param_type{1,end+1}=param_type;
model.param_name_equiv{1,end+1}=names_equiv;
%
model_new=model;
%
return

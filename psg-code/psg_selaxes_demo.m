%psg_selaxes_demo: select and keep specific axes from a btc dataset
% 
%  This applies only to btc datasets, and produces a modified coords and setup file
%  Intended to process one or more files from a single setup file, e.g., bgca3pt
%
%   Writes new data files and setup files, but note that the stimuli in the data file :
%     are reordered to conform to ordering in setup (since this reordering is done when data files are read)
%
% Selection is based on specs field of metadata file
%   The stimuli that are "on axis" are the stimuli that have at most one nonzero value in specs
%   (other values may be not specified or within tol of zero)
%
%  See also: PSG_GET_COORDSETS, BTC_DEFINE.
%
if ~exist('opts_read') opts_read=struct();end %for psg_get_coordsets
if ~exist('opts_write') opts_write=struct(); end %for psg_write_coorddata
%
if ~exist('tol') tol=10^-5; end
%
meta_field_select={'specs','spec_labels','typenames','btc_augcoords','btc_methods'}; %fields that need to be selected based on retained axes
%
opts_read.input_type=1;
opts_read.if_log=1;
opts_read.if_warn=0;
opts_read.if_data_only=1;
opts_read.nfiles_max=Inf;
%
if_debug=getinp('1 for debugging mode','d',[0 1]);
if (if_debug)
    debug_string='-test';
else
    debug_string=[];
end
%
tag_coords='_coords_';
%
dict=btc_define;
%
if ~exist('axes_filters')
    axes_filters={'bgca3','dgea3','bdce3','bc6','bcpm3','bcpm24','bc55q','bcpp55q','bcpm55q','bcmp55q','bcmm55q'}; %file name has a pt at the end
end
for k=1:length(axes_filters)
    disp(sprintf(' %1.0f->axis filter %s',k,axes_filters{k}));
end
axes_filter=axes_filters{getinp('choice','d',[1 length(axes_filters)])};
ui_filter=cat(2,axes_filter,'pt*_*coords_*.mat');
%
axes_avail=[];
for k=1:length(axes_filter) %keep in original order
    if ismember(axes_filter(k),dict.codel)
        axes_avail=cat(2,axes_avail,axes_filter(k));
    end
end
axes_keep_sorted=getinp(sprintf('axes to keep (subset of %s)',axes_avail),'s',[]);
axes_keep=[];
for k=1:length(axes_filter)
    if ismember(axes_filter(k),axes_keep_sorted)
        axes_keep=cat(2,axes_keep,axes_filter(k));
    end
end
if length(axes_keep)==length(axes_avail)
    disp('Warning: the newly created files will have the same data and metadata as the original files.')
end
%    
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(setfield(opts_read,'ui_filter',ui_filter));
nsets=length(sets);
%
if_write=getinp('1 to write the files','d',[0 1]);
setup_written=[];
for iset=1:nsets
    disp(' ');
    nstims=size(ds{iset}{1},1);
    ndims=length(ds{iset});
    data_fullname=opts_read_used{iset}.data_fullname;
    data_fullname=strrep(data_fullname,'/','\');
    data_namestart=max([0,max(find(data_fullname=='/')),max(find(data_fullname=='\'))]);
    data_name=data_fullname((1+data_namestart):end);
    % typical data_name: bdce3pt_coords_ZK-gp_sess01_20.mat
    % typical data_new:  bc3pt_coords_ZK-gp_sess01_20-bdce.mat
    data_new=cat(2,strrep(data_name,axes_avail,axes_keep));
    data_new=strrep(data_new,'.mat',cat(2,'-',axes_filter,debug_string,'.mat')); %put axes available at end
    data_fullnew=cat(2,data_fullname(1:data_namestart),data_new);
    disp(sprintf(' original  data full name: %s',data_fullname));
    disp(sprintf(' original  data brief name: %s',data_name));
    disp(sprintf('processed  data brief name: %s',data_new));
    disp(sprintf('processed  data full name: %s',data_fullnew));
    %
    setup_fullname=opts_read_used{iset}.setup_fullname;
    setup_fullname=strrep(setup_fullname,'/','\');
    setup_namestart=max([0,(find(data_fullname=='/')),max(find(data_fullname=='\'))]);
    setup_name=setup_fullname((1+data_namestart):end);
    setup_new=cat(2,strrep(setup_name,axes_avail,axes_keep));
    setup_new=strrep(setup_new,'.mat',cat(2,debug_string,'.mat'));
    setup_fullnew=cat(2,setup_fullname(1:setup_namestart),setup_new);
    disp(sprintf(' original setup full name: %s',setup_fullname));
    disp(sprintf(' original setup brief name: %s',setup_name));
    disp(sprintf('processed setup brief name: %s',setup_new));
    disp(sprintf('processed setup full name: %s',setup_fullnew));
    %
    %the following are now aligned in the same order:
    % ds{iset}
    % sas{iset}.typenames
    % sas{iset}.spec_labels
    % sas{iset}.specs
    %
    %determine which data entries to keep
    nstims=size(sas{iset}.specs,1);
    keep_flags=zeros(nstims,1);
    coords_used=cell(nstims,1);
    for istim=1:nstims
        fns=fieldnames(sas{iset}.specs{istim});
        for ifn=1:length(fns)
            specified=fns{ifn};
            val=sas{iset}.specs{istim}.(specified);
            if abs(val)>tol
                coords_used{istim}=cat(2,coords_used{istim},specified);
            end
        end
        keep_flags(istim)=all(ismember(coords_used{istim},axes_keep));
        if (if_debug)
            disp(sprintf(' entry %2.0f: typename %15s, spec_label %15s, coords_used: %5s, keep_flag %2.0f',...
                istim,...
                sas{iset}.typenames{istim},sas{iset}.spec_labels{istim},coords_used{istim},keep_flags(istim)));
        end
    end
    %check consistency with prevoius
    if_ok=1;
    if (iset==1)
        keep_flags_init=keep_flags;
    else
        if (length(keep_flags)~=length(keep_flags_init))
            warning('different numbers of stimuli across datasets');
            if_ok=0;
        elseif any(keep_flags~=keep_flags_init)
            warning('inconsistent sets of stimuli to keep across datasets');
            if_ok=0;
        end
    end
    if if_ok==1
        inds_keep=find(keep_flags);
        %
        % data
        %
        ds_new=ds{iset};
        for idim=1:length(ds_new)
            if ~isempty(ds_new{idim})
                ds_new{idim}=ds_new{idim}(inds_keep,:);
            end
        end
        sout=struct;
        sout.stim_labels=strvcat(sas{iset}.typenames(inds_keep,:)); %since we have reordered the coordinates
        %metadata
        s_orig=getfield(load(setup_fullname),'s'); %reload setup file
        s_new=s_orig;
        for k=1:length(meta_field_select)
            fn=meta_field_select{k};
            s_new.(fn)=s_orig.(fn)(inds_keep,:);
        end
        s_new.nstims=length(inds_keep);
        sas_new.s=s_new;
        %
        if (if_write)
            opts_write_used=psg_write_coorddata(data_fullnew,ds_new,sout,opts_write);
            disp('data file written');
            if strcmp(setup_written,setup_fullnew)
                disp('setup file already written');
            else
                setup_written=setup_fullnew;
                save(setup_fullnew,'-struct','sas_new');
                disp('setup file written')
            end
        end
    else
        disp('keep flags for this dataset')
        disp(keep_flags(:)')
        disp('keep flags for initial dataset')
        disp(keep_flags_init(:)')
    end
end
%

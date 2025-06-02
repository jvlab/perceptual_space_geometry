function [sets_align,ds_align,sas_align,ovlp_array,sa_pooled,opts_used]=psg_align_coordsets(sets,ds,sas,opts)
% [ds_align,sas_align,ovlp_array,sa_pooled,opts_used]=psg_align_coordsets(ds,sas,opts)
% aligns datasets and metadata that have partially overlapping stimuli
% Here, "alignment" refers to matching up the stimuli at the level of the metadata, not a coordinate rotation
%
% sas{k}.typenames is used to establish stimulus identity
%  other fields of sas{k} corresponding to individual stimuli are reordered 
%  fields of sas{k} NOT corresponding to individual stimuli are unchanged
% coordinates of ds are reordered (sorted alphabetically) and values without overlaps are replaced by NaN's
%
% sets: cell array (one cell for each dataset) dataset descriptors (typically from psg_get_coordsets)
% ds: cell array (one cell for each dataset) of coordinate data (typically from psg_get_coordsets)
%  All datasets must have dimension lists beginning at 1 and without gaps
% sas: cell array of metadata (typically from psg_get_coordsets)
% opts: options
%   opts.if_log: 1 to log progress
%   opts.if_btc_specoords_remake: 1 to treat sas{iset}.btc_specoords as a pooled variable, 0 to not, [] to determine from
%       whether all of the sas.btc_specoords are square matrices of 0's and 1's of dimension equal to nstims_each(iset)
%       default is []; legacy behavior is 0; proper function for mater, domain and auxiliary classes requires [] or 1.
%       This will set remake btc_specoords to be the identity matrix of the pooled stimulus set.
%   opts.min: minimum number of datasets that must contain a stimulus, in order for the stimulus to be included
%       default is 1 (legacy behavior: all stimuli used), can also be 'any'; 
%       'all': stimuli must be present in all datasets to be kept
%
% sets_align: cell array (one cell for each dataset) dataset descriptors (typically from psg_get_coordsets) after alignment
% ds_align: cell array (one cell for each dataset) of coordinate data after alignment
% sas_align: cell array (one cell for each dataset) of metadata data after alignment
% ovlp_array: [nstims_kept nsets]: array of 1's for points in overlaps, 0 otherwise; stimulus order corresponds to sa_pooled
% sa_pooled: pooled metadata
%    ovlp_array, [sets|ds|sas]_align and sa_pooled only refer to kept stimuli by the opts.min criterion.
%    If opts.min='all', there are no NaN's in da_align, and sas_align, sa_pooled are all identical
%    If opts.min=1, typically there are NaN's in da_align corresponding to stimuli present in other datasets.
%    In either case sas_align{*} and sa_pooled match except possibly for b
%    (nstims_all refers to all stimuli present in any dataset, regardless of opts.min, nstims_kept are the number of stimuli kept)
%    
% opts_used: options used
%   opts_used.which_common [nstims_all,nsets] points to the original stimuli for each set, has zeros elsewhere
%   opts_used.ovlp_array_all=double(opts_used.which_common>0);
%   opts_used.which_common_kept point to the original stimuli for each set,
%     but only has rows corresponding to stimuli that meet the opts.min criterion
%      Note that ovlp_array=double(opts_used.which_common_kept>0)
%
% 05May24: added documentation about btc_specoords
% 24Nov24: removal of restriction on having same number of stimuli for mater, domain, and auxiliary classes
%     by detecting that btc_specoords is only 0s and 1s. Also add opts.if_btc_specoords_remake.
% 01Jun24: add opts.min
%
%  See also: PSG_ALIGN_KNIT_DEMO, PSG_GET_COORDSETS, PSG_READ_COORDDATA, PROCRUSTES_CONSENSUS_PTL_TEST, PSG_DEFOPTS,
%   PSG_REMNAN_COORDSETS.
%
opts_fields=psg_defopts();
%
if (nargin<=3) opts=struct; end
opts=filldefault(opts,'if_log',0);
opts=filldefault(opts,'if_btc_specoords_remake',[]);
opts=filldefault(opts,'min',1);
%
nsets=length(ds);
sets_align=sets;
ds_align=ds;
sas_align=sas;
if (opts.if_log)
    disp(sprintf('alignments attempted with %2.0f datasets',nsets));
end
typenames_all=cell(0);
nstims_each=zeros(nsets,1);
for iset=1:nsets
    nstims_each(iset)=sets{iset}.nstims;
    typenames_all=[typenames_all;sas{iset}.typenames];
end
typenames_all=unique(typenames_all);
nstims_all=length(typenames_all);
if (opts.if_log)
    for iset=1:nsets
        disp(sprintf(' set %3.0f: stimuli: %3.0f, type %10s, label %s',iset,sets{iset}.nstims,sets{iset}.type,sets{iset}.label));
    end
    disp(sprintf(' unique typenames: %3.0f',nstims_all));
end
which_common=zeros(nstims_all,nsets);
opts_used.warnings=[];
for iset=1:nsets
    for istim=1:nstims_each(iset)
        typename_to_match=sas{iset}.typenames{istim};
        find_common=strmatch(typename_to_match,typenames_all,'exact');
        if length(find_common)==1
            which_common(find_common,iset)=istim;
        elseif isempty(find_common)
            wmsg=sprintf('in set %3.0f: stim %s not found in common list',iset,typename_to_match);
            opts_used.warnings=strvcat(opts_used.warnings,wmsg);
            warning(wmsg);
        else
            wmsg=sprintf('in set %3.0f: stim %s has multiple matches',iset,typename_to_match);
            opts_used.warnings=strvcat(opts_used.warnings,wmsg);
            warning(wmsg);
        end
    end
end
%
switch opts.min
    case 'any'
        nmin=1;
    case 'all'
        nmin=nsets;
    otherwise
        nmin=opts.min;
end
ovlp_array_all=double(which_common>0);
noccur=sum(ovlp_array_all,2);
ptrs_kept=find(noccur>=nmin);
nstims_kept=length(ptrs_kept);
which_common_kept=which_common(ptrs_kept,:);
ovlp_array=double(which_common_kept>0);
if opts.if_log
    disp(sprintf(' typenames present in at least %2.0f datasets: %3.0f',nmin,nstims_kept));
end
%
%determine whether to treat if_btc_specoords as a pooled variable
%(do so if it is always a matrix of 0's and 1's, and of dimension equal to
%number of stimuli, or if specified)
if isempty(opts.if_btc_specoords_remake)
    if_btc_specoords_remake=1;
    for iset=1:nsets
        if isfield(sas{iset},'btc_specoords')
            z=sas{iset}.btc_specoords;
            %square, size nstims_each(iset), and only 0 and 1?
            if size(z,1)~=nstims_each(iset) | size(z,2)~=nstims_each(iset) | any(~ismember(z(:),[0 1]))
                if_btc_specoords_remake=0;
            end
         end
    end
    if (opts.if_log)
        disp(sprintf('if_btc_specoords_remake=%1.0f (determined from metadata)',if_btc_specoords_remake));
    end
else
    if_btc_specoords_remake=opts.if_btc_specoords_remake;
    if (opts.if_log)
        disp(sprintf('if_btc_specoords_remake=%1.0f (provided)',if_btc_specoords_remake));
    end
end
opts_used.if_btc_specoords_remake=if_btc_specoords_remake;
%
fields_align=opts_fields.fields_align;
fields_pool=opts_fields.fields_pool;
fields_remake=opts_fields.fields_remake;
if if_btc_specoords_remake %add to the remake list and remove from the align list
   fields_remake{end+1}='btc_specoords';
   idx=strmatch('btc_specoords',fields_align,'exact');
   fields_align=fields_align(setdiff(1:length(fields_align),idx));
end
opts_used=opts;
opts_used.fields_align=fields_align;
opts_used.fields_pool=fields_pool;
opts_used.fields_remake=fields_remake;
%
opts_used.which_common=which_common;
opts_used.which_common_kept=which_common_kept;
%adjust metadata
fields_all_vals=struct;
fields_all_vals.nstims=nstims_all;
fields_all_vals.typenames=typenames_all;
if if_btc_specoords_remake
    fields_all_vals.btc_specoords=eye(nstims_all);
end
%
typenames_kept=typenames_all(ptrs_kept);
fields_kept_vals=struct;
fields_kept_vals.nstims=nstims_kept;
fields_kept_vals.typenames=typenames_kept;
if if_btc_specoords_remake
    fields_kept_vals.btc_specoords=eye(nstims_kept);
end
opts_used.fields_all_vals=fields_all_vals;
opts_used.fields_kept_vals=fields_kept_vals;
opts_used.ovlp_array_all=ovlp_array_all;
%
sa_pooled=struct;
for iset=1:nsets
    %modify overall set descriptors
    sets_align{iset}.nstims=nstims_kept;
    %modify metadata
    fns=fieldnames(sas{iset});
    for ifn=1:length(fns)
        fn=fns{ifn};
        if ~isempty(strmatch(fn,fields_align,'exact')) %align the data from this field
            if iscell(sas{iset}.(fn)) %1d cell
                sas_align{iset}.(fn)=cell(nstims_kept,1);
                if ~isfield(sa_pooled,fn)
                    sa_pooled.(fn)=cell(nstims_kept,1);
                end
                for istim=1:nstims_kept
                    if which_common_kept(istim,iset)>0
                        to_fill=sas{iset}.(fn){which_common_kept(istim,iset)};
                        sas_align{iset}.(fn){istim}=to_fill;
                        sa_pooled.(fn){istim}=to_fill;
                    end
                end
            else %multi-dim array, first index is stimulus id
                if ~isfield(sa_pooled,fn)
                    sa_pooled.(fn)=[];
                end
                sas_align{iset}.(fn)=psg_align_coordsets_do(sas{iset}.(fn),which_common_kept(:,iset));
                sa_pooled.(fn)=psg_align_coordsets_do(sas{iset}.(fn),which_common_kept(:,iset),sa_pooled.(fn));
            end
        elseif ~isempty(strmatch(fn,fields_pool,'exact')) |  ~isempty(strmatch(fn,fields_remake,'exact'))%use data from all datasets combined, kept stimuli only
            sas_align{iset}.(fn)=fields_kept_vals.(fn);
            sa_pooled.(fn)=fields_kept_vals.(fn);
        else %other fields, just copy
            sas_align{iset}.(fn)=sas{iset}.(fn);
            if ~isfield(sa_pooled,fn)
                sa_pooled.(fn)=sas{iset}.(fn);
            end
        end
    end  
    % modify coordinates
    for idim=1:length(ds{iset})
        ds_align{iset}{idim}=psg_align_coordsets_do(ds{iset}{idim},which_common_kept(:,iset));
    end
end
return

function aligned=psg_align_coordsets_do(to_align,pointers,prev) %align a multidimensional variable based on first dimension
if (nargin<=2)
    prev=[];
end
nstims_each=size(to_align,1);
nstims_all=size(pointers,1);
dims_all=size(to_align);
dims_last=dims_all(2:end);
var_toalign=reshape(to_align,[nstims_each prod(dims_last)]); %make .(fn) 2-dimensional
if isempty(prev)
    var_aligned=nan(nstims_all,prod(dims_last));
else
    var_aligned=reshape(prev,[nstims_all prod(dims_last)]);
end
for istim=1:nstims_all
    if pointers(istim)>0 %works for >2-d variables
        var_aligned(istim,:)=var_toalign(pointers(istim),:);
    end
end
aligned=reshape(var_aligned,[nstims_all dims_last]);
return

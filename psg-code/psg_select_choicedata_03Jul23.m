function [d_select,sa_select,sel_used,desc_used,opts_used]=psg_select_choicedata(d,sa,sel,desc,opts)
% [d_select,sa_select,sel_used,desc_used,opts_used]=psg_select_choicedata(d,sa,sel,desc,opts)
% selects a subset of conditions from choice data 
%
% first 3 columns of d_select (the triad tokens) and all fields of sa are remapped accordingly
%
% d:  has size [ntriads 5], columns are [ref s1 s2 choices(d(ref,s1)<d(ref,s2)) total trials]
% sa: selected fields of stimulus specification structure s, and also
%    btc_augcoords, augmented coords, of size [s.nstims nbtc] 
%    btc_specoords, specificed coords, of size [s.nstims nbtc]
%   typically read by psg_read_choicedata
%
% sel: a string consisting of patterns in sa.typenames to be selected.  If more than one, concatenate and separate by |.
%   if empty, this is requested interactively
% desc: a description for the selection, best with no spaces, e.g, 'a_only',
%   if empty, this is requested interactively
% opts: options
%   opts.if_log: 1 to log
%
% d_select: array of 5 columns, containing the triads with only the selected stimli
% sa_select: sa, with all fields remapped
% sel_used: list_sel, or interactively-input list_sel
% desc_used: dec, or interactively-input desc
% opts_used:
%       opts_used.remapping : a list of length length(sa.typenames), indicating the mapping of each stimlus in d to d_select
%         Tokens that are not remapped have an entry of zero.
%       opts_used.original: a list of length equal to the number of typenames selected,
%         indicating (via an entry of 1 to length(sa.typenames) the source
%       nstims: number of stimuli (length(original)
%
% 03Jul23: modifications for compatibility with faces_mpi; defaults to not approve of an empty selection
% 
% See also: PSG_READ_CHOICEDATA, PSG_TENTLIKE_DEMO, PSG_UMI_TRIPLIKE_DEMO.
%
xfr_fields={'nchecks','nsubsamp','opts_psg','btc_dict'}; %fields to transfer
sel_fields={'specs','spec_labels','typenames','btc_augcoords','btc_specoords'}; %fields to permute
%
if (nargin<5)
    opts=struct;
end
if (nargin<3)
    sel=[];
end
if (nargin<4)
    desc=[];
end
opts=filldefault(opts,'if_log',0);
%
if isempty(sel)
    ifok=0;
    while (ifok==0)
        for k=1:length(sa.typenames)
            disp(sprintf(' typename %2.0f: %s',k,sa.typenames{k}));
        end
        sel=getinp('selection string or multiple strings, separated by |','s',[]);
        list_sel=psg_select_util(sel,sa);
        if isempty(list_sel)
            disp('selection yields empty list');
            ifok=getinp('1 if ok','d',[0 1],0);
        else
            disp('typenames selected:');
            disp(list_sel);
            ifok=getinp('1 if ok','d',[0 1]);
        end
    end
else
    list_sel=psg_select_util(sel,sa);
end
remapping=zeros(length(sa.typenames),1);
original=zeros(size(list_sel,1));
for ilist=1:length(list_sel)
    idx=strmatch(list_sel{ilist},sa.typenames,'exact');
    if ~isempty(idx)
        if length(idx)==1
            remapping(idx)=ilist;
            original(ilist)=idx;
        elseif length(idx)>1
            warning(sprintf('selection %s matches more than one typename',list_sel{ilist}));
        else
            warning(sprintf('selection %s does not match any typename',list_sel{ilist}));
        end
    end
end
if opts.if_log
    disp('remapping:');
    disp(remapping(:)');
    for k=1:length(sa.typenames)
        if remapping(k)==0
            disp(sprintf('%12s -> unselected',sa.typenames{k}));
        else
            disp(sprintf('%12s -> token %2.0f',sa.typenames{k},remapping(k)));
        end
    end
end
if isempty(desc)
    desc=getinp('descriptor','s',[]);
end
%
sa_select=struct;
sa_select.nstims=length(original);
%transfer some fields from sa without change
for ifield=1:length(xfr_fields)
    if isfield(sa,xfr_fields{ifield})
        sa_select.(xfr_fields{ifield})=sa.(xfr_fields{ifield});
    end
end
%
%select some fields from sa
%
for ifield=1:length(sel_fields)
    sel_name=sel_fields{ifield};
    if isfield(sa,sel_name)
        x=sa.(sel_name);
        sa_select.(sel_name)=x(original,:); %works for cells and numeric arrays
    end
end
%
%remap first three columns of d, containing the triad
%
triad_cols=[1:3];
d_select=d;
d_select(:,triad_cols)=remapping(d(:,triad_cols));
d_keep=find(~any(d_select(:,triad_cols)==0,2));
d_select=d_select(d_keep,:);
%
if opts.if_log
    disp(sprintf('selected %4.0f of %4.0f tokens, keeping %4.0f of %4.0f triads',length(original),sa.nstims,size(d_select,1),size(d,1)));
end
%
opts.remapping=remapping;
opts.original=original;
opts.nstims=length(original);
sel_used=sel;
desc_used=desc;
opts_used=opts;
return

function list_sel=psg_select_util(sel_parsed,sa)
%select typenames containing any of the substrings of sel_parsed, and keep the same order as in sa.typenames
list_sel_unsorted=cell(0);
while ~isempty(sel_parsed)
    [list_string,sel_parsed]=strtok(sel_parsed,'|');   
    for k=1:length(sa.typenames)
        if contains(sa.typenames{k},list_string)
            list_sel_unsorted{end+1,1}=sa.typenames{k};
        end
    end
end
list_sel_unsorted=unique(list_sel_unsorted);
%reorder list_sel_unsorted in the order in which they occurred in sa.typenames
list_sel=cell(0);
for k=1:length(sa.typenames)
    idx=strmatch(sa.typenames{k},list_sel_unsorted,'exact');
    if ~isempty(idx)
        list_sel{end+1}=sa.typenames{k};
    end
end
return

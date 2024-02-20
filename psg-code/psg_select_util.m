function [list_sel,inds]=psg_select_util(sel_parsed,sa)
% [list_sel,inds]=psg_select_util(sel_parsed,sa) is a utility to select typenames containing 
% any of the substrings of sel_parsed, and keep the same order as in sa.typenames
%
% sel_parsed: a selection string, consisting of one or more segments separated by |
% sa: a metadata structure,  typically read by psg_read_choicedata;
%   only field that matters is sa.typenames{:}
%
% list_sel: a cell array, with selected typenames in same order as sa.typenames
%   that contains any of the strings in sel_parsed
% inds: indexes of list_sel into sa.typenames
%
% 19Feb24: modularize and add inds
%
% See also:  PSG_SELECT_CHOICEDATA.
% 
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
inds=[];
for k=1:length(sa.typenames)
    idx=strmatch(sa.typenames{k},list_sel_unsorted,'exact');
    if ~isempty(idx)
        list_sel{end+1}=sa.typenames{k};
        inds=[inds,k];
    end
end
return

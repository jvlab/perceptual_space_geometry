function [sa_restrict,opts_used]=psg_meta_restrict(sa,include_flags,opts)
% [sa_restrict,opts_used]=psg-meta_restrict(sa,include_flags,opts) restricts a metadata
% structure sa (typically regurned by psg_read_coorddata)  to specific stimuli
% and optionally reorders
%
% sa: original dataset
% include_flags: a binary vector, length must be number of stimuli in sa
% opts:
%  opts.sort_order: sort_order(k) is pointer to source of output stimulus k (1:sum(include_flags if omitted)
%  opts.fields_cell: fields of sa to be treated as cells
%  opts.fields_array: fields of sa to be treaded as 2d arrays
%
%   See also: PSG_COORD_PIPE_PROC, PSG_READ_COORDDATA.
%
if (nargin<=2)
    opts=struct;
end
opts=filldefault(opts,'fields_cell',{'specs','spec_labels','typenames'});
opts=filldefault(opts,'fields_array',{'typenames','specs','spec_labels','btc_augcoords','btc_specoords'});
opts=filldefault(opts,'sort_order',[1:sum(include_flags)]);
opts_used=opts;
sa_restrict=sa;
sa_restrict.nstims=sum(include_flags);
for k=1:length(opts.fields_cell)
    fn=opts.fields_cell{k};
    if isfield(sa,fn)
        temp=sa.(fn)(include_flags>0);
        sa_restrict.(fn)=temp(opts.sort_order);
    end
end
for k=1:length(opts.fields_array)
    fn=opts.fields_array{k};
    if isfield(sa,fn)
        temp=sa.(fn)(include_flags>0,:);
        sa_restrict.(fn)=temp(opts.sort_order,:);
    end
end
return
end
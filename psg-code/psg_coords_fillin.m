function [s_filled,ds_filled]=psg_coords_fillin(s,ds,opts);
%[s_filled,ds_filled]=psg_coords_fillin(s,ds,opts) is a utility
%to fill in lower-dimensional coordinates that could be missing from a coordinate dataset
% 
% s: the sets metadata structure read by psg_get_coordsets
% ds: a cell array of coordinates, ds{k}=the k-dimensional model
% opts: opts.if_log=1 to log any changes (default), -1 to log always, 0 not to log
%
% s_filled: s with s.dimlist updated to reflect presence of all coords
% ds_filled: ds, with each prevoiusly empty ds{k} filled from lower coords of next available higher dimension
%
%  See also:  PSG_COORD_PIPE_PROC, PSG_GET_COORDSETS.
%
if (nargin<=2)
    opts=struct;
end
opts=filldefault(opts,'if_log',1);
if isfield(s,'label_long') 
    label=s.label_long;
elseif isfield(s,'label')
    label=s.label;
else
    label='unlabeled dataset';
end
dim_max=max(s.dim_list);
s_filled=s;
ds_filled=ds;
for id=1:dim_max
    if ismember(id,s.dim_list) & ~isempty(ds{id})
        if opts.if_log==-1
            disp(sprintf('%s: coords for model dim %2.0f exist',label,id));
        end
    elseif ~ismember(id,s.dim_list) & isempty(ds{id})
        dim_fill=min(s.dim_list(s.dim_list>id));
        s_filled.dim_list=sort([id,s_filled.dim_list]);
        ds_filled{id}=ds{dim_fill}(:,1:id);
        if opts.if_log~=0
            disp(sprintf('%s: coords for model dim %2.0f filled from dim %2.0f',label,id,dim_fill));
        end
    else
        disp(sprintf('%s: coord metadata for model dim %2.0f are inconsistent with data',label,id))
    end
end %id
return

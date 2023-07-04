function [legstring,ou]=psg_spec2legend(sa,point,opts)
%[legstring,ou]=psg_spec2legend(sa,point,opts) creates a short legend
%string from a specifiction in a metadata structure sa
%
% if sa.spec_labels is present, it is used 
% if not and sa.typenames is present, it is used
% otherwise, 'data' followed by point index
%
% sa: metadata structure returned by psg_read_coorddata or psg_read_choicedata
% point: the index of the point to label
% opts: options (not used at present)
%
% legstring: legend string
% ou: options used
%
% See also: PSG_PLOTCOORDS, PSG_READ_COORDDATA, PSG_PLOTANGLES.
%
if (nargin<=2)
    opts=struct;
end
ou=opts;
legstring=[];
if isfield(sa,'spec_labels')
    legstring=sa.spec_labels{point};
    %some shortening if faces_mpi
    if ~isempty(findstr(legstring,'id=')) & ...
        ~isempty(findstr(legstring,' age=')) & ...
        ~isempty(findstr(legstring,' gender=')) & ...
        ~isempty(findstr(legstring,' emo=')) & ...
        ~isempty(findstr(legstring,' set='))
        legstring=legstring(findstr(legstring,' age='):end); %remove id=xxx
        legstring=strrep(legstring,'gender','gen'); %shorten gender
        legstring=legstring(1:findstr(legstring,' set=')); %remove set=
    end
elseif isfield(sa,'typenames')
    legstring=sa.typenames{point};
else
    legstring=sprintf('data 1.0f',point);
end
return

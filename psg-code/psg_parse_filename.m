function [paradigm_type,paradigm_name,subj_id,file_type]=psg_parse_filename(filename)
% [paradigm_type,paradigm_name,subj_id,file_type]=psg_parse_filename(filename)
% parses a file name to extract paradigm type, paradigm name, and subject id
%
% filename: file name
%
% paradigm_type: 'animals' or 'btc'
% paradigm_name: for animals: 'word','object','intermdiate_object','intermediate_texture','texture'
%                for btc: bgca3pt, bc6pt, ...
% subj_id: 2 or  chars
% file_type: 'coords' or 'choices'
%
%  See also:  PSG_LIKE_MAKETABLE, PSG_PROCRUSTES_REGR_DEMO.
underscore='_';
%
filename=cat(2,filename,'.mat');
filename=strrep(filename,'.mat.mat','.mat');
underscores=find(filename=='_');
paradigm_type=[];
paradigm_name=[];
subj_id=[];
file_type=[];
if ~isempty(strfind(filename,'pt_'))
    paradigm_type='btc';
    paradigm_name=filename(1:underscores(1)-1);
    subj_id=filename(underscores(2)+1:underscores(3)-1);
else
    paradigm_type='animals';
    paradigm_name=filename(underscores(1)+1:underscores(end)-1);
    subj_id=filename(1:underscores(1)-1);
end
if ~isempty(strfind(filename,'coords'))
    file_type='coords';
end
if ~isempty(strfind(filename,'choices'))
    file_type='choices';
end
return

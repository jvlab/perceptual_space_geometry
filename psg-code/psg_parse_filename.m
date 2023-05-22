function p=psg_parse_filename(filename)
% p=psg_parse_filename(filename) parses a file name to extract 
% paradigm type, paradigm name, and subject id
%
% filename: file name
%
% p.paradigm_type: 'animals' or 'btc'
% p.paradigm_name: for animals: 'word','object','intermdiate_object','intermediate_texture','texture'
%                for btc: bgca3pt, bc6pt, ...
% p.subj_id: 2 or more uppercase chars
% p.file_type: 'coords' or 'choices'
% p.extra: expected extra characters
%
%  See also:  PSG_LIKE_MAKETABLE, PSG_PROCRUSTES_REGR_DEMO.
%
underscore='_';
%
%remove path if present
filename=strrep(filename,'/',filesep);
filename=strrep(filename,'\',filesep);
filename=cat(2,filesep,filename);
filename=filename(max(find(filename==filesep))+1:end);
%remove .mat if present
filename=cat(2,filename,'.mat');
filename=strrep(filename,'.mat.mat','.mat');
filename=strrep(filename,'.mat','');
%
underscores=find(filename==underscore);
paradigm_type=[];
paradigm_name=[];
subj_id=[];
file_type=[];
extra=[];
if ~isempty(strfind(filename,'pt_')) & length(underscores)>=3 %filename like bcpm3pt_coords_SAW_sess01_10
    paradigm_type='btc';
    paradigm_name=filename(1:underscores(1)-1);
    subj_id=filename(underscores(2)+1:underscores(3)-1);
    extra=filename(underscores(3)+1:end);
elseif length(underscores)>=1
    paradigm_type='animals'; %filename like ZK_intermediate_object_choices
    paradigm_name=filename(underscores(1)+1:underscores(end)-1);
    subj_id=filename(1:underscores(1)-1);
end
if ~isempty(strfind(filename,'coords'))
    file_type='coords';
end
if ~isempty(strfind(filename,'choices'))
    file_type='choices';
end
p=struct;
p.paradigm_type=paradigm_type;
p.paradigm_name=paradigm_name;
p.subj_id=subj_id;
p.file_type=file_type;
p.extra=extra;
return

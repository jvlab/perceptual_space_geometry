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
% 24Jun23: add compatibility with btcsel paradigm type
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
if ~isempty(strfind(filename,'pt_')) & length(underscores)>=3 %filename  bcpm3pt_coords_SAW_sess01_10 or 'bc6pt_choices_MC_sess01_10_sel_cp_rand'
    paradigm_type='btc';
    paradigm_name=filename(1:underscores(1)-1);
    idx=strfind(filename,'_sel');
    if ~isempty(idx)
        paradigm_type='btcsel';
        paradigm_name=cat(2,paradigm_name,filename(idx+4:end));
    end
    subj_id=filename(underscores(2)+1:underscores(3)-1);
    extra=filename(underscores(3)+1:end);
elseif ~isempty(strfind(filename,'faces_')) %%file name like faces_mpi_en2_fc_choices_MC_sess01_10_sel__
    paradigm_type='faces';
    %paradigm name may have variable number of underscores, but ends with choices or coords
    filename_mod=strrep(strrep(filename,'choices','*'),'coords','*'); %now like faces_mpi_en2_fc_*_MC_sess01_10_sel__
    underscores_mod=find(filename_mod==underscore);
    paradigm_name=filename(underscores_mod(1)+1:(find(filename_mod=='*')-2));
    filename_end=filename_mod(find(filename_mod=='*')+2:end); %like MC_sess01_10_sel__
    underscores_end=find(filename_end==underscore);
    subj_id=filename_end(1:underscores_end(1)-1);
    extra=filename_end(underscores_end(1)+1:end);
    sel_start=min(strfind(filename_end,'sel_'));
    if ~isempty(sel_start)
        paradigm_name=cat(2,paradigm_name,'-',filename_end(sel_start+4:end));
    end
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

function sets=psg_make_setstruct(type_string,dim_list,label_long,nstims,pipeline,opts)
% sets=psg_make_setstruct(type_string,dim_list,label_long,nstims,pipeline,opts) is a
% utility function to create a metadata structure
%
% type_string: typically 'data' or 'qform'
% dim_list: list of dimensions available
% label_long: long dataset label, typically file name and path
% nstims: number of stimuli
% pipeline: a processing pipeline structure, can be omitted
% opts: options
%
%  See also: PSG_GET_COORDSETS.
%
if (nargin<=5)
    pipeline=struct;
end
if (nargin<=6)
    opts=struct;
end
sets=struct;
sets.type=type_string;
sets.dim_list=dim_list;
sets.nstims=nstims;
sets.label_long=label_long;
sets.label=sets.label_long;
sets.label=strrep(sets.label,'./','');
sets.label=strrep(sets.label,'.mat','');
sets.label=strrep(sets.label,'coords_','');
sets.pipeline=pipeline;
return
end

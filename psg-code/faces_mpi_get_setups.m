% faces_mpi_get_setups: generate a list of setups for faces_mpi psg experiments
%
%   See also:  FACES_MPI_PSG_SETUP, FACES_MPI_INVENTORY, FACES_MPI_SETUP_CREATE.
%
%  faces_mpi_attrib_info should be defined, e.g., by faces_mpi_inventory.
%
% attributes are set up in the order of attrib_info.psg_order={'gender','age','emo','set'}
%
% modeled loosely after spokes_setup_create.
%
faces_mpi_setups=cell(0);
%
%a variable not listed is taken to be the fujll list in faces_mpi_attrib_info.psg_order
%
faces_mpi_setups{1}.lists.emo={'n'};
faces_mpi_setups{1}.count_each=2;
%
faces_mpi_setups{2}=faces_mpi_setups{1};
faces_mpi_setups{2}.count_each=3;
%
faces_mpi_setups{3}.lists.gender={'f'};
faces_mpi_setups{3}.lists.emo={'n'};
faces_mpi_setups{3}.lists.set={'a'};
faces_mpi_setups{3}.count_each=8;
%
faces_mpi_setups{4}=faces_mpi_setups{3};
faces_mpi_setups{4}.count_each=12;
%
faces_mpi_setups{5}=faces_mpi_setups{3};
faces_mpi_setups{5}.lists.gender={'m'};
faces_mpi_setups{6}=faces_mpi_setups{4};
faces_mpi_setups{6}.lists.gender={'m'};
%
faces_mpi_setups{7}.lists.set={'a'};
faces_mpi_setups{7}.count_each=1;
%
for ig=1:length(faces_mpi_attrib_info.gender.vals)
    for ia=1:length(faces_mpi_attrib_info.age.vals)
        faces_mpi_setups{end+1}.lists.set={'a'};
        faces_mpi_setups{end}.lists.age{1}=faces_mpi_attrib_info.age.vals(ia);
        faces_mpi_setups{end}.lists.gender{1}=faces_mpi_attrib_info.gender.vals(ig);
        faces_mpi_setups{end}.count_each=6;
    end
end
%
name_string_header='faces_mpi';
list_names=faces_mpi_attrib_info.psg_order;
list_lengths=ones(1,length(list_names));
for isetup=1:length(faces_mpi_setups)
    name_string=name_string_header;
    name_string_brief=[];
%   fill in unsupplied fields
    for il=1:length(list_names)
        list_name=list_names{il};
        if ~isfield(faces_mpi_setups{isetup}.lists,list_name)
            faces_mpi_setups{isetup}.lists.(list_name)=cellstr(faces_mpi_attrib_info.(list_name).vals)';
        end
        list_vals=faces_mpi_setups{isetup}.lists.(list_name);
        name_string=cat(2,name_string,'_',list_name,'-',sprintf('%s',list_vals{:}));
        faces_mpi_setups{isetup}.name=name_string;
        list_lengths(il)=length(list_vals);
        if list_lengths(il)~=faces_mpi_attrib_info.(list_name).nlevels
            name_string_brief=cat(2,name_string_brief,list_name(1),sprintf('%s',list_vals{:}));
        end
    end
    name_string=cat(2,name_string,'_',sprintf('x%1.0f',faces_mpi_setups{isetup}.count_each));
    name_string_brief=cat(2,name_string_brief,sprintf('%1.0f',faces_mpi_setups{isetup}.count_each));
    faces_mpi_setups{isetup}.name=name_string;
    faces_mpi_setups{isetup}.name_brief=name_string_brief;
    faces_mpi_setups{isetup}.nstims=prod(list_lengths)*faces_mpi_setups{isetup}.count_each;
end
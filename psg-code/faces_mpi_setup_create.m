function [faces_mpi_setups,faces_mpi_info,opts_used]=faces_mpi_setup_create(faces_mpi_array,faces_mpi_attrib_info,faces_mpi_table,opts)
% [faces_mpi_setups,faces_mpi_info,opts_used]=faces_mpi_setup_create(faces_mpi_array,faces_mpi_attrib_info,faces_mpi_table,opts)
% create some face setups for mpi dataset
%
% faces_mpi_array: binary array, size [182 3 2 6 2] indicating availablility of stimuli by id_num, age group, gender, expression, and set
% faces_mpi_attrib_info: structure describing attributes (dims 2-5 of faces_mpi_array)
% faces_mpi_table: table with above information
%  [above input arguments typically found in faces_mpi_inventory.mat]
% opts: options
%   how_rand: typically 'default', can be 'shuffle'.  Note that random number generator state is restored at end
%   if_log: 1 to log (default)
%   id_num_exclude: list of face IDs to exclude, defaults to [];
%
% faces_mpi_setups: setups for psg experiments
% faces_mpi_info: information about the setups
% opts_used: options used
%
%   See also:  FACES_MPI_PSG_SETUP, FACES_MPI_INVENTORY, FACES_MPI_GET_SETUPS.
%
if nargin<=3
    opts=struct;
end
opts=filldefault(opts,'how_rand','default');
opts=filldefault(opts,'if_log',1);
opts=filldefault(opts,'id_num_exclude',[]);
rand_state=rng;
rng(opts.how_rand);
%
dim_list=struct();
lev_list=struct();
attribs=faces_mpi_attrib_info.psg_order;
natts=length(attribs);
for iatt=1:natts
    dim_list.(attribs{iatt})=faces_mpi_attrib_info.(attribs{iatt}).dim;
    lev_list.(attribs{iatt})=faces_mpi_attrib_info.(attribs{iatt}).nlevels;
end
%
% create faces_mpi_setups
%
faces_mpi_get_setups;
nsetups=length(faces_mpi_setups);
n_needed=zeros(lev_list.gender,lev_list.age,nsetups);
%log and determine how many individual ids are needed for each setup
for isetup=1:nsetups
    if opts.if_log
        disp(sprintf(' setup %3.0f: %4.0f stims, name (%12s): %s',isetup,faces_mpi_setups{isetup}.nstims,...
            faces_mpi_setups{isetup}.name_brief,faces_mpi_setups{isetup}.name));
    end
    lists=faces_mpi_setups{isetup}.lists;
    for igender=1:lev_list.gender
        for iage=1:lev_list.age
            if any(contains(lists.age,faces_mpi_attrib_info.age.vals(iage))) & ...
                any(contains(lists.gender,faces_mpi_attrib_info.gender.vals(igender))) 
                n_needed(igender,iage,isetup)=faces_mpi_setups{isetup}.count_each;  
            end
        end
    end
end
faces_mpi_info.n_needed=n_needed;
%
dim_names=fieldnames(dim_list);
if (opts.if_log)
    for idim=1:natts
        dn=dim_names{idim};
        disp(sprintf(' %12s is dimension %2.0f in faces_mpi_array and has %3.0f levels',dn,dim_list.(dn),lev_list.(dn)));
    end
end
%
%find which faces have all images available
%
f_array=sum(sum(faces_mpi_array,dim_list.emo),dim_list.set);
f_array=permute(f_array,[1,dim_list.gender,dim_list.age]);
f_need_all=size(faces_mpi_array,dim_list.emo)*size(faces_mpi_array,dim_list.set);
faces_haveall=cell(lev_list.gender,lev_list.age);
faces_haveall_cum=[];
n_excluded=zeros(lev_list.gender,lev_list.age);
n_haveall=zeros(lev_list.gender,lev_list.age);
for igender=1:lev_list.gender
    for iage=1:lev_list.age
        faces_haveall_tent=find(f_array(:,igender,iage)==f_need_all);
        faces_haveall_exclude=intersect(faces_haveall_tent,opts.id_num_exclude);
        faces_haveall{igender,iage}=setdiff(faces_haveall_tent,faces_haveall_exclude);
        faces_haveall_cum=[faces_haveall_cum;faces_haveall{igender,iage}];
        n_excluded(igender,iage)=length(faces_haveall_exclude);
        n_haveall(igender,iage)=length(faces_haveall{igender,iage});
        if (opts.if_log)
            disp(sprintf(' gender   %s age   %s  number of subjects with all individuals available: %4.0f (%3.0f excluded)',...
                faces_mpi_attrib_info.gender.vals(igender),faces_mpi_attrib_info.age.vals(iage),n_haveall(igender,iage),n_excluded(igender,iage)));
        end
    end
end
if (opts.if_log)
    disp(sprintf(' gender %s age %s  number of subjects with all individuals available: %4.0f (%3.0f excluded)',...
        'all','all',sum(n_haveall(:)),sum(n_excluded(:))));
end
faces_mpi_info.n_excluded=n_excluded;
faces_mpi_info.n_haveall=n_haveall;
faces_mpi_info.faces_haveall=faces_haveall;
% 
% restore random number generator state
rng(rand_state);
opts_used=opts;
%
return

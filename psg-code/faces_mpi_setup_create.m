function [faces_mpi_setups,faces_mpi_info,opts_used]=faces_mpi_setup_create(faces_mpi_array,faces_mpi_attrib_info,faces_mpi_table,opts)
% [faces_mpi_setups,faces_mpi_info,opts_used]=faces_mpi_setup_create(faces_mpi_array,faces_mpi_attrib_info,faces_mpi_table,opts)
% create some face setups for mpi dataset, including the randomization
%
% faces_mpi_array: binary array, size [182 3 2 6 2] indicating availablility of stimuli by id_num, age group, gender, expression, and set
% faces_mpi_attrib_info: structure describing attributes (dims 2-5 of faces_mpi_array)
% faces_mpi_table: table with above information
%  [above input arguments typically found in faces_mpi_inventory.mat]
% opts: options
%   how_rand: typically 'default', can be 'shuffle'.  Note that random number generator state is restored at end
%   if_log: 1 to log (default), 2 for more extensive logging
%   id_num_exclude: list of face IDs to exclude, defaults to [];
%
% faces_mpi_setups: setups for psg experiments
% faces_mpi_info: information about the setups
% opts_used: options used
%
%   See also:  FACES_MPI_PSG_SETUP, FACES_MPI_INVENTORY, FACES_MPI_GET_SETUPS, ZPAD.
%
if nargin<=3
    opts=struct;
end
opts=filldefault(opts,'how_rand','default');
opts=filldefault(opts,'if_log',1);
opts=filldefault(opts,'id_num_exclude',[]);
opts=filldefault(opts,'id_digits',3);
opts=filldefault(opts,'img_type','jpg');
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
    lists=faces_mpi_setups{isetup}.lists;
    for igender=1:lev_list.gender
        for iage=1:lev_list.age
            if any(contains(lists.age,faces_mpi_attrib_info.age.vals(iage))) & ...
                any(contains(lists.gender,faces_mpi_attrib_info.gender.vals(igender))) 
                n_needed(igender,iage,isetup)=faces_mpi_setups{isetup}.count_each;  
            end
        end
    end
    if opts.if_log
        disp(sprintf(' setup %3.0f: %4.0f stims (%12s) %50s, need: y(f,m) %2.0f %2.0f m(f,m) %2.0f %2.0f o(f,m) %2.0f %2.0f',...
            isetup,faces_mpi_setups{isetup}.nstims,...
            faces_mpi_setups{isetup}.name_brief,faces_mpi_setups{isetup}.name,...,
            reshape(n_needed(:,:,isetup),[1 lev_list.gender*lev_list.age])));
    end
end
faces_mpi_info.n_needed=n_needed;
%
dim_names=fieldnames(dim_list);
if (opts.if_log>=2)
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
ids_avail=cell(lev_list.gender,lev_list.age);
ids_avail_cum=[];
n_excluded=zeros(lev_list.gender,lev_list.age);
n_avail=zeros(lev_list.gender,lev_list.age);
for igender=1:lev_list.gender
    for iage=1:lev_list.age
        faces_haveall_tent=find(f_array(:,igender,iage)==f_need_all);
        faces_haveall_exclude=intersect(faces_haveall_tent,opts.id_num_exclude);
        ids_avail{igender,iage}=setdiff(faces_haveall_tent,faces_haveall_exclude);
        ids_avail_cum=[ids_avail_cum;ids_avail{igender,iage}];
        n_excluded(igender,iage)=length(faces_haveall_exclude);
        n_avail(igender,iage)=length(ids_avail{igender,iage});
        if (opts.if_log>=2)
            disp(sprintf(' gender   %s age   %s  number of indivs with all images available: %4.0f (%3.0f excluded)',...
                faces_mpi_attrib_info.gender.vals(igender),faces_mpi_attrib_info.age.vals(iage),n_avail(igender,iage),n_excluded(igender,iage)));
        end
    end
end
if (opts.if_log>=2)
    disp(sprintf(' gender %s age %s  number of indivs with all images available: %4.0f (%3.0f excluded)',...
        'all','all',sum(n_avail(:)),sum(n_excluded(:))));
end
faces_mpi_info.n_excluded=n_excluded;
faces_mpi_info.n_avail=n_avail;
faces_mpi_info.ids_avail=ids_avail;
%
%choose a subset of setups
%
ifok=0;
while (ifok==0)
    setup_choices=getinp('choice(s) for setups to make','d',[1 nsetups]);
    n_needed_tot=sum(n_needed(:,:,setup_choices),3);
    if all(n_needed_tot(:)<=n_avail(:))
        ifok=1;
    else
        disp('not enough face ids available.');
    end
end
faces_mpi_info.setup_choices=setup_choices;
%
%select the ids, removing them from the available list, and create file names
%
ids_avail_after=ids_avail;
ids_used=[];
for isetup_ptr=1:length(setup_choices)
    isetup=setup_choices(isetup_ptr);
    fms=faces_mpi_setups{isetup};
    disp(sprintf(' creating files for setup %2.0f->%s',isetup,fms.name));
    ids_use=cell(lev_list.gender,lev_list.age);
    filenames=cell(fms.nstims,1);
    typenames=cell(fms.nstims,1);
    specs=cell(fms.nstims,1);
    spec_labels=cell(fms.nstims,1);
    istim=0;
    for igender=1:lev_list.gender
        for iage=1:lev_list.age
            if n_needed(igender,iage,isetup)>0
                select_vec=zeros(1,length(ids_avail_after{igender,iage}));
                select_vec([1:n_needed(igender,iage,isetup)])=1;
                select_vec=select_vec(randperm(length(ids_avail_after{igender,iage})));
                ids_use{igender,iage}=ids_avail_after{igender,iage}(find(select_vec==1));
                ids_avail_after{igender,iage}=setdiff(ids_avail_after{igender,iage},ids_use{igender,iage});
                ids_used=[ids_used;ids_use{igender,iage}(:)];
                for iemo=1:length(fms.lists.emo)
                    for iset=1:length(fms.lists.set)
                        for iuse=1:length(ids_use{igender,iage})
                            istim=istim+1;
                                fbase=cat(2,zpad(ids_use{igender,iage}(iuse),opts.id_digits),'_',...
                                    faces_mpi_attrib_info.age.vals(iage),'_',...
                                    faces_mpi_attrib_info.gender.vals(igender),'_',...
                                    faces_mpi_attrib_info.emo.vals(iemo),'_',...
                                    faces_mpi_attrib_info.set.vals(iset));
                                filenames{istim}=cat(2,fbase,'.',opts.img_type);
                                typenames{istim}=fbase; %
                                specs{istim}.id=ids_use{igender,iage}(iuse);
                                specs{istim}.gender=faces_mpi_attrib_info.gender.vals(igender);
                                specs{istim}.age=faces_mpi_attrib_info.age.vals(iage);
                                specs{istim}.emo=faces_mpi_attrib_info.emo.vals(iemo);
                                specs{istim}.set=faces_mpi_attrib_info.set.vals(iset);
                                spec_labels{istim}=cat(2,'id=',zpad(ids_use{igender,iage}(iuse),opts.id_digits),' ',...
                                    'age=',faces_mpi_attrib_info.age.vals(iage),' ',...
                                    'gender=',faces_mpi_attrib_info.gender.vals(igender),' ',...
                                    'emo=',faces_mpi_attrib_info.emo.vals(iemo),' ',...
                                    'set=',faces_mpi_attrib_info.set.vals(iset));
                        end
                    end %iset
                end %iemo
            end %needed?
        end %iage
    end %igender
    faces_mpi_setups{isetup}.ids_use=ids_use;
    faces_mpi_setups{isetup}.filenames=filenames;
    faces_mpi_setups{isetup}.typenames=typenames;
    faces_mpi_setups{isetup}.specs=specs;
    faces_mpi_setups{isetup}.spec_labels=spec_labels;
end
faces_mpi_info.ids_avail_after=ids_avail_after;
faces_mpi_info.ids_used=ids_used;
%
%generate the file names and the lists of used ids
%
% 
% restore random number generator state
rng(rand_state);
opts_used=opts;
%
return

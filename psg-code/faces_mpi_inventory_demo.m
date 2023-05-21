%faces_mpi_inventory_demo: demonstrate and use faces_mpi_inventory.
% an inventory of the contents of the MPI faces database
% 
%Face images from https://faces.mpdl.mpg.de/imeji/  Main citation:
% Ebner, N. C., Riediger, M., & Lindenberger, U., (2010).
% FACES—A database of facial expressions in young, middle-aged, and older men and women:
% Development and validation. Behavior Research Methods, 42(1), 351-362. https://doi.org/10.3758/BRM.42.1.351
%
% See also: FILLDEFAULT, FACES_MPI_INVENTORY.
%
if ~exist('faces_mpi_opts') faces_mpi_opts=[]; end
faces_mpi_opts=filldefault(faces_mpi_opts,'faces_mpi_path','C:\Users\jdvicto\Documents\jv\EY7977\faces\mpi_faces\FacesSetASetB\');
faces_mpi_opts=filldefault(faces_mpi_opts,'if_log',1);
faces_mpi_opts=filldefault(faces_mpi_opts,'if_log_details',0);
if ~exist('inventory_filename') inventory_filename='faces_mpi_inventory.mat'; end
%
if getinp('1 to create inventory mat-file from scratch','d',[0 1],0)
    [faces_mpi_table,faces_mpi_array,faces_mpi_attrib_info,faces_mpi_opts_used]=faces_mpi_inventory(faces_mpi_opts);
    if getinp('1 to save inventory','d',[0 1],1)
        inventory_filename=getinp('inventory file name','s',[],inventory_filename);
        save(inventory_filename,'faces_mpi_table','faces_mpi_array','faces_mpi_attrib_info','faces_mpi_opts_used');
    end
else
    inventory_filename=getinp('inventory file name','s',[],inventory_filename);
    load(inventory_filename);
end
disp(sprintf('total entries in database from table: %5.0f  from array: %5.0f',...
    size(faces_mpi_table,1),sum(faces_mpi_array(:))));
disp(fieldnames(faces_mpi_attrib_info));
faces_mpi_eachsubj=sum(sum(faces_mpi_array,4),5);
n_ages=faces_mpi_attrib_info.age.nlevels;
n_genders=faces_mpi_attrib_info.gender.nlevels;
this_ag_any=cell(n_ages,n_genders);
this_ag_all=cell(n_ages,n_genders);
n_all=prod([size(faces_mpi_array,4),size(faces_mpi_array,5)]);
for iage=1:n_ages
    for igender=1:n_genders
        this_ag_any{iage,igender}=find(faces_mpi_eachsubj(:,iage,igender)>0);
        this_ag_all{iage,igender}=find(faces_mpi_eachsubj(:,iage,igender)==n_all);
        disp(sprintf(' for age %s gender %s, %3.0f subjects, %3.0f subjects have a full set of images',...
            faces_mpi_attrib_info.age.vals(iage),...
            faces_mpi_attrib_info.gender.vals(igender),...
            length(this_ag_any{iage,igender}),...
            length(this_ag_all{iage,igender})));
    end
end

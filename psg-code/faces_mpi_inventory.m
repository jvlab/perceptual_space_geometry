function [table_faces_mpi,avail_array,attrib_info,opts_used]=faces_mpi_inventory(opts)
% [table_faces_mpi,avail_array,attrib_info,opts_used]=faces_mpi_inventory(opts)
% shows and returns an inventory of the contents of the MPI faces database and returns
% 
%Face images from https://faces.mpdl.mpg.de/imeji/  Main citation:
% Ebner, N. C., Riediger, M., & Lindenberger, U., (2010).
% FACES—A database of facial expressions in young, middle-aged, and older men and women:
% Development and validation. Behavior Research Methods, 42(1), 351-362. https://doi.org/10.3758/BRM.42.1.351
%
% opts:
%    opts.faces_mpi_path: path to the faces files
%    opts.if_log: 1 to log
%    opts.if_log_details: 1 to log details (counts for face database by multiple criteria)
%
% table_faces_mpi: table of face files and attributs
% avail_array: binary array indicating what is available: dim 1: id number, dim 2-5: attributes
% attrib_info: information about the attributes (number of levels, level names)
% opts_used: options used
%
% See also: FILLDEFAULT, FACES_MPI_INVENTORY_DEMO, FACES_MPI_PSG_SETUP.
%
if (nargin<1)
    opts=struct;
end
opts=filldefault(opts,'faces_mpi_path','C:\Users\jdvicto\Documents\jv\EY7977\faces\mpi_faces\FacesSetASetB\');
opts=filldefault(opts,'if_log',1);
opts=filldefault(opts,'if_log_details',0);
%
opts_used=opts;
%
table_faces_mpi=table();
avail_array=[];
attribs={'age','gender','emo','set'};
attrib_info=struct();
under='_';
nattribs=length(attribs);
attrib_levels=zeros(1,nattribs);
nunder_expected=length(attribs);
%
filenames=[];
separator_locs=zeros(0,nunder_expected+1);
file_count=0;
id_nums=zeros(0);
%
doscmd=cat(2,'dir ',opts.faces_mpi_path,'*.jpg /b');
[errs,dirtext]=dos(doscmd);
opts_used.errs=errs;
%
if errs~=0
    disp('errors encountered');
    return
else
    disp(' parsing directory');
    remain=dirtext;
    while ~isempty(remain)
        [filename,remain]=strtok(remain);
        if ~isempty(remain)
            file_count=file_count+1;
            filenames=strvcat(filenames,filename);
            under_ptrs=findstr(filename,under);
            if (length(under_ptrs)~=nunder_expected)
                disp(sprintf(' file %4.0f has unexpected format: %s',file_count,filename));
            else
                separator_locs(file_count,1:nunder_expected)=under_ptrs;
                separator_locs(file_count,nunder_expected+1)=min(findstr(filename,'.jpg'));
                id_nums(file_count,1)=str2num(filename(1:under_ptrs(1)-1));
            end               
        end
    end
    nfiles=file_count;
    disp(sprintf(' %4.0f files found with valid names.  Creating a table.',nfiles));
    table_faces_mpi=[table(filenames),table(id_nums)];
    table_faces_mpi.Properties.VariableNames={'filename','id_num'};
    %parse attributes into table
    attlist=cell(nattribs,1);
    for iatt=1:length(attribs)
        attlist{iatt}=[];
        for ifile=1:nfiles
            attlist{iatt}=strvcat(attlist{iatt},filenames(ifile,separator_locs(ifile,iatt)+1:separator_locs(ifile,iatt+1)-1));
        end
        table_faces_mpi=[table_faces_mpi,table(attlist{iatt})];
        table_faces_mpi.Properties.VariableNames{end}=attribs{iatt};
        if (opts.if_log)
            disp(sprintf('for attribute %s, unique values are:',attribs{iatt}));
        end
        attlist_u=unique(attlist{iatt});
        counts=zeros(1,length(attlist_u));
        for iu=1:length(attlist_u)
            counts(iu)=sum(attlist{iatt}==attlist_u(iu));
            if (opts.if_log)
                disp(sprintf(' value %2s: count %5.0f',attlist_u(iu),counts(iu)));
            end
        end
        attrib_levels(iatt)=length(attlist_u);
        attrib_info.(attribs{iatt}).vals=attlist_u;
        attrib_info.(attribs{iatt}).nlevels=attrib_levels(iatt);
        attrib_info.(attribs{iatt}).counts=counts;
        %then report, for each attribute, pair, and triplet, how many of each kind
    end
    %fill avail_array
    opts_used.avail_array_dims=cat(2,'id_num',attribs);
    avail_array=zeros([max(id_nums) attrib_levels]);
    for ifile=1:nfiles
        for iatt=1:length(attribs)
            att_coord(iatt)=find(attlist{iatt}(ifile)==attrib_info.(attribs{iatt}).vals);
        end
        avail_array(id_nums(ifile),att_coord(1),att_coord(2),att_coord(3),att_coord(4))=1;
    end
    %details
    if opts.if_log_details==1
        for ncombs=1:nattribs
            disp(sprintf('tallies for combinations of %1.0f attributes',ncombs));
            comb_list=nchoosek([1:nattribs],ncombs);
            for icomb=1:size(comb_list,1);
                attrib_sel=comb_list(icomb,:);
                sel_header_string=[];
                for sel_ptr=1:ncombs
                    sel_header_string=cat(2,sel_header_string,...
                        sprintf('  %7s: %3.0f levels',attribs{attrib_sel(sel_ptr)},attrib_levels(attrib_sel(sel_ptr))));
                end
                disp(sel_header_string)
                %mult-index based on attrib_levels
                nsel=prod(attrib_levels(attrib_sel));
                att_level_sel=zeros(1,ncombs);
                for isel=1:nsel
                    sel_string=[];
                    att_comb_ptr=isel-1;
                    val_sel=cell(0);
                    table_sel=table_faces_mpi; %start with a full table
                    for sel_ptr=1:ncombs
                        att_level_sel(1,sel_ptr)=1+mod(att_comb_ptr,attrib_levels(attrib_sel(sel_ptr)));
                        att_comb_ptr=floor(att_comb_ptr/attrib_levels(attrib_sel(sel_ptr)));
                        val_sel{sel_ptr}=attrib_info.(attribs{attrib_sel(sel_ptr)}).vals(att_level_sel(1,sel_ptr)); %value to select for this attribute
                        sel_string=cat(2,sel_string,...
                            sprintf('  %7s  level %1.0f: %s',attribs{attrib_sel(sel_ptr)},...
                            att_level_sel(sel_ptr),val_sel{sel_ptr}));
                        %take intersection and then add this to the string
                        table_sel=table_sel(strmatch(val_sel{sel_ptr},table_sel.(attribs{attrib_sel(sel_ptr)}),'exact'),:);
                    end
                    %disp(att_level_sel);
                    disp(sprintf('%s:  total entries: %7.0f',sel_string,size(table_sel,1)));
                end
            end
        end
    end %if_log_details
end
return

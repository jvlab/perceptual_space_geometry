function [r,p_used,f_used,s_used,fname_read,t_used,t_string]=btc_soid_summ_read(path_name,file_template,subj_id,type,if_auto)
% [r,p_used,f_used,s_used,fname_read,t_used,t_string]=...
%    btc_soid_summ_read(path_name,file_template,subj_id,type,if_auto)
% reads files of quadratic fits to psychophsical data, written by btc_soid_demo.
%
% path_name: path name, '.' if empty
% file_template: file template, btc_allraysfixedb_*_100surrs_madj if empty
% subj_id: subject ID, 'avg' if empty
% type: setup type, 1 to 12, defaults to max available if empty (assume symmetric, combine similar axes, augment coords by maxent)
% if_auto:  if = 1, then no console input requeseted if there is just one answer or template is provided (defaults to 0)
%
%   all arguments requested at console if empty
%
% r: results structure from chosen setup type
% p_used: path_name used
% f_used: file template used
% s_used: subject ID used
% fname_read: file name read
% t_used: type used
% t_string: descriptor of type used\
%
%   See also:  BTC_SOID_DEMO.
%
r=[];
if (nargin<1) path_name=[]; end
if (nargin<2) file_template=[]; end
if (nargin<3) subj_id=[]; end
if (nargin<4) type_no=[]; end
if (nargin<5) if_auto=0; end
if isempty(path_name)
    path_name='.';
    path_name=getinp('path name for btc quadratic-fit data','s',[],path_name);
end
p_used=path_name;
if isempty(file_template)
    file_template='btc_allraysfixedb_*_100surrs_madj';
end
if (if_auto==0) | isempty(file_template)
    file_template=getinp('file name template for btc quadratic-fit data','s',[],file_template);
end
f_used=file_template;
%
doscmd=cat(2,'dir ',path_name,filesep,file_template,'*.mat/B');
if (if_auto==0)
    disp('looking for files...');
end
[dos_errs,all_files]=dos(doscmd);
if (dos_errs>0)
    disp(all_files)
    disp('cannot find files.');
    return
end
all_files=deblank(all_files);
ffl=cell(0); %file list
subj_ids=cell(0);
nfiles=0;
while ~isempty(all_files)
    [fname,all_files]=strtok(all_files);
    blpos=findstr(fname,'_');
    if length(blpos)>=3 %expect that subject name will follow second underscore
        subj_id_this=fname((blpos(2)+1):(blpos(3)-1));
        if isempty(subj_id)
            nfiles=nfiles+1;
            ffl{nfiles}=fname;
            subj_ids{nfiles}=subj_id_this;
        else
            if strcmp(subj_id,subj_id_this)
                nfiles=nfiles+1;
                ffl{nfiles}=fname;
                subj_ids{nfiles}=subj_id_this;
            end
        end
    end
end
if (nfiles<1)
    disp('cannot find files that match template (for requested subject id, if given');
    return
end
if isempty(subj_id)
    subj_id_list=unique(subj_ids);
    if (if_auto==0) | length(subj_id_list)>1
        for is=1:length(subj_id_list)
            disp(sprintf('%1.0f->%s',is,subj_id_list{is}));
        end
        subj_id_ptr=getinp('subject choice','d',[1 length(subj_id_list)]);
    else
        subj_id_ptr=1;
    end
    subj_id=subj_id_list{subj_id_ptr};
    %reselect from file list
    nfiles_subj=0;
    files_subj=cell(0);
    for ifile=1:nfiles
        if strcmp(subj_id,subj_ids{ifile})
            nfiles_subj=nfiles_subj+1;
            files_subj{nfiles_subj}=ffl{ifile};
        end
    end
else
    files_subj=ffl;
    nfiles_subj=nfiles;
end
s_used=subj_id;
if (nfiles_subj<1)
    disp('cannot find files that have this subject id.');
    return
else
    if (if_auto==0) | nfiles_subj>1
        for ifile=1:nfiles_subj
            disp(sprintf('%1.0f->%s',ifile,files_subj{ifile}));
        end
        ifile_ptr=getinp('choice','d',[1 nfiles_subj]);
    else
        ifile_ptr=1;
    end
    fname_read=files_subj{ifile_ptr};
    results=getfield(load(cat(2,path_name,filesep,fname_read)),'r');
end
if isempty(type)
    for itype=1:length(results)
        disp(sprintf('type %2.0f ->%s',itype,results{itype}.setup.label));
    end
    if (if_auto==0) | length(results)>1
        type=getinp('choice','d',[1 length(results)],length(results));
    else
        type=1;
    end
end
t_used=type;
t_string=results{type}.setup.label;
r=results{type};
return
end


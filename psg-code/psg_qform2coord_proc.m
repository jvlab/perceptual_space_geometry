%psg_qform2coord_proc creates a coordinate dataset in similar format to experimental data
% from a quadratic form model
%
% The output datafile has the fields
%   dim[a], where [a] is a character string 1,...,7,8,9,10,
%   stim_labels, as in the coord data files
%   pipeline, indicating the processing pipeline to create this file
%    pipeline{n} indicates the most recent processing step; here, n=1, since this is the first
%    step that creates a coord dataset
%
%  See also: PSG_GET_COORDSETS, PSG_QFORMPRED, PSG_READ_COORD_DATA. PSG_COORD_PIPE_PROC, PSG_WRITE_COORDDATA.
%
if ~exist('fname_suggest_base') fname_suggest_base='*_coords_QFM_sess01_01.mat'; end %suggested template for output file name
%
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_qpred') opts_qpred=struct(); end %for psg_qformpred
opts_qpred=filldefault(opts_qpred,'qform_datafile_def','../stim/btc_allraysfixedb_avg_100surrs_madj.mat');
opts_qpred=filldefault(opts_qpred,'qform_modeltype',12); %if_symm=1 (symmetrize around origin), if_axes=1 (symmetrize bc, de, tuvw); ifaug=1 (augmented coords)
%
disp('This will create a coordinate dataset from a btc setup file and a quadratic form model.');
opts_read.input_type=2;
[sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred,1);
%
ndim_max=length(ds{1});
%
%assemble the output file
%
sout=struct;
sout.stim_labels=strvcat(sas{1}.typenames);
sout.pipeline{1}.type='qform2coord';
sout.pipeline{1}.sets=sets{1};
sout.pipeline{1}.opts.opts_read_used=opts_read_used{1};
sout.pipeline{1}.opts.opts_qpred_used=opts_qpred_used{1};
disp(sprintf('output structure with %2.0f stimuli and up to %2.0f dimensions created via %s',size(sout.stim_labels,1),ndim_max,sets{1}.label));
%
setup_name=opts_read_used{1}.setup_fullname;
setup_name=strrep(setup_name,'9.mat','');
setup_name=strrep(setup_name,'.mat','');
fname_suggest=strrep(fname_suggest_base,'*',setup_name);
fname=getinp(sprintf('file name, e.g., %s',fname_suggest_base),'s',[],fname_suggest);
%
opts_write=struct;
opts_write.data_fullname_def=fname;
opts_write.if_log=1;
opts_write_used=psg_write_coorddata([],ds{1},sout,opts_write);

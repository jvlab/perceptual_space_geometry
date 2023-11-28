function opts_used=psg_write_coorddata(data_fullname,ds,sout,opts)
% opts_used=psg_write_coorddata(data_fullname,ds,sout,opts) writes a coordinate structure
%
% data_fullname: data file name including path; if empty, requested interactively
% ds: cell array, ds{k} is the scaling solution for the kth dimension
%   (typically, ds{k} has size [nstims k])
% sout: structure with metadata
%   sout.stim_labels: strvcat of the stimulus labels, size(sout.stim_labels)=nstims
%   sout.pipeline: optional
% opts:
%   opts.data_fullname_def: default file name to write
%   opts.if_log: 1 to log
%
% opts_used: options used
%
% 28Nov23: multiple ways to determine number of stimuli
%
% See also:  PSG_READ_COORDDATA, PSG_QFORM2COORD_PROC.
%
if nargin<=3
    opts=struct;
end
opts=filldefault(opts,'data_fullname_def','./psg_data/bgca3pt_coords_QFM_sess01_01.mat');
opts=filldefault(opts,'if_log',1);
%
if isempty(data_fullname)
    data_fullname=getinp('full path and file name of data file','s',[],opts.data_fullname_def);
end
opts.data_fullname=data_fullname;
opts_used=opts;
%
if isfield(sout,'stim_labels')
    nstims=size(sout.stim_labels,1);
elseif isfield(sout,'nstims')
    nstims=sout.nstims;
elseif isfield(sout,'typenames')
    nstims=sout.typenames;
end
dim_string=[];
%
for idimptr=1:length(ds)
    idim=size(ds{idimptr},2);
    dim_string=cat(2,dim_string,sprintf('%1.0f ',idim));
    if size(ds{idimptr},1)~=nstims
        warning(sprintf('for dim %2.0f (pointer=%2.0f), number of stimuli found is %2.0f, expected: %2.0f',idim,idimptr,size(ds{idimptr},1),nstims));
    end   
    dname=cat(2,'dim',sprintf('%1.0f',idim));
    sout.(dname)=ds{idim};
end
dim_string=deblank(dim_string);
save(data_fullname,'-struct','sout');
if opts.if_log
    disp(sprintf('%s written with %2.0f stimuli and dimensions %s.',data_fullname,nstims,dim_string));
end
return


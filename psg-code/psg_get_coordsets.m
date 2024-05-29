function [sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred,nsets)
% [sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred,nsets)
% reads in one or more coordinate sets, either from a data file or via quadratic form prediction
% and checks for consistency of number of stimuli and stimulus typenames
%
% opts_read: options for psg_read_coorddata, can be empty or omitted
%   opts_read.if_log: 1 to log
%   opts_read.input_type: 0 for either, 1 forces expemental data, 2 forces quadratic form, can be a scalar, or an array that is cycled throuigh for each dataset
%   opts_read.data_fullnames: cell array of data file full names; if empty, will be requested
%   opts_read.setup_fullnames: cell array of setup file full names; if empty, will be requested
%   opts_read.if_auto: 1 not to ask if ok, and to use all defaults for qform model
% opts_rays: options for psg_findrays, can be empty or omitted
% opts_qpred: options for psg_qformpred, can be empty or omitted. as well
%            as default options for qform_datafile and qform_modeltype
% nsets: if present, number of datasets. otherwise requested from console
%   if negative, or entered as negative, then a dialog box is used to load files
%
% sets: cell array, sets{iset} is a structure of that describes the dataset
% ds: cell array, ds{iset}{nd} is a structure of coordinates (ntrials x nd), other fields have available dimensions and a label
% sas: cell array, sas{iset} is the setup structure returned by psg_read_coorddata
% rayss: cell array, rays{iset} is the ray structure returned by psg_findrays
% opts_read_used, opts_rays_used, opts_qpred_used: cell arrays of options used for each set
%
% 05Jan23: shortened sets{iset}.label; preserved original label as sets{iset}.label_long
% 27Jun23: override mode with single-point for rays, and suppressing ray angle calculation and plotting (for bcpm24pt and similar)
% 28Jun23: invoke psg_findray_setopts for ray defaults
% 25Sep23: add opts_read.input_type, option to force either experimental data (1) or model data (2)
% 27Sep23: add pipeline field to sets (always empty for qform model)
% 22Feb24: localization params now from psg_localopts
% 23Feb24: add data_fullnames, setup_fullnames, if_auto
% 28Apr24: add if_data_only (set to allow only experimental data to be read)
% 04May24: nsets can be negative, allowing for a dialog box to load multiple datasets
% 29May24: fixed bug if nsets=-1
%
%  See also: PSG_PROCRUSTES_DEMO, PSG_FINDRAYS, PSG_QFORMPRED, PSG_READ_COORDDATA, PSG_VISUALIZE_DEMO,
% PSG_CONSENSUS_DEMO, PSG_FINDRAY_SETOPTS, PSG_LOCALOPTS, PSG_COORD_PIPE_PROC, PSG_COORDS_FILLIN.
%
if nargin<1
    opts_read=struct();
end
if nargin<2
    opts_rays=struct();
end
if nargin<3
    opts_qpred=struct();
end
if nargin<4
    nsets=[];
end
opts_local=psg_localopts;
input_types={'experimental data','qform model'};
%
opts_read=filldefault(opts_read,'if_log',1);
opts_read=filldefault(opts_read,'input_type',0);
opts_read=filldefault(opts_read,'data_fullnames',cell(0));
opts_read=filldefault(opts_read,'setup_fullnames',cell(0));
opts_read=filldefault(opts_read,'if_auto',0);
opts_read=filldefault(opts_read,'if_data_only',0);
%
opts_qpred=filldefault(opts_qpred,'qform_datafile_def',opts_local.model_filename_def); %model file name
opts_qpred=filldefault(opts_qpred,'qform_modeltype',opts_local.modeltype_def);  %model type
%
if opts_read.if_data_only==1 %if only data can be read, override requests for any other input type
    opts_read.input_type=1;
end
if_ok=0;
while (if_ok==0)
    if isempty(nsets) | nsets==0
        nsets=getinp('number of datasets (negative: use a dialog box for multiple datasets)','d',[-100 100]);
    end
    nsets_pos=abs(nsets);
    if_dialog=double(nsets<0);
    if (if_dialog)
        if_dialog_ok=0;
        while (if_dialog_ok==0)
            [filenames_short,pathname]=uigetfile('*coords*.mat',sprintf('Select %1.0f coordinate files',nsets_pos),'Multiselect','on');
            if ~iscell(filenames_short) filenames_short={filenames_short}; end
            nfiles_sel=length(filenames_short);
            if_dialog_ok=double(nfiles_sel==nsets_pos);
            opts_read.input_type=1;
            if ~iscell(filenames_short)
                if filenames_short==0 %allow for an exit
                    if_dialog_ok=1;
                    if_dialog=0;
                end
            end
        end
    end
    %
    sets=cell(1,nsets_pos);
    ds=cell(1,nsets_pos);
    sas=cell(1,nsets_pos);
    rayss=cell(1,nsets_pos);
    opts_read_used=cell(1,nsets_pos);
    opts_rays_used=cell(1,nsets_pos);
    opts_qform_used=cell(1,nsets_pos);
    d_qform=cell(1,nsets_pos);
    d_mds=cell(1,nsets_pos);
    data_fullname=[];
    setup_fullname=[];
    %
    %read the datasets and/or the qform models
    %   
    for iset=1:nsets_pos
        disp(' ');
        disp(sprintf(' entering set %2.0f of %2.0f:',iset,nsets_pos));
        input_type_use=opts_read.input_type(mod(iset-1,length(opts_read.input_type))+1);
        if input_type_use==0
            input_type_use=getinp('1 for experimental data, 2 for qform model','d',[1 2]);
        else
            disp(sprintf('dataset %1.0f is %s',iset,input_types{input_type_use}));
        end
        if (if_dialog)
            data_fullname=cat(2,pathname,filenames_short{iset});
        else
            if nsets_pos<=length(opts_read.data_fullnames)
                data_fullname=opts_read.data_fullnames{iset};
            end
        end
        if nsets_pos<=length(opts_read.setup_fullnames)
            setup_fullname=opts_read.setup_fullnames{iset};
        end
        switch input_type_use
            case 1
                sets{iset}.type='data';
                [ds{iset},sas{iset},opts_read_used{iset},pipeline]=psg_read_coorddata(data_fullname,setup_fullname,opts_read);
                %determine whether one of the strings in ray_minpts_default
                %is present in setup file name, and if so, use this to
                %determine the default for opts_rays
                opts_rays_use=psg_findray_setopts(opts_read_used{iset}.setup_fullname,opts_rays);
                %
                [rayss{iset},opts_rays_used{iset}]=psg_findrays(sas{iset}.btc_specoords,setfield(opts_rays_use,'permute_raynums',opts_read_used{iset}.permute_raynums));
                opts_read.setup_fullname_def=opts_read_used{iset}.setup_fullname;
                sets{iset}.dim_list=opts_read_used{iset}.dim_list;
                sets{iset}.nstims=sas{iset}.nstims;
                sets{iset}.label_long=opts_read_used{iset}.data_fullname;
                sets{iset}.label=sets{iset}.label_long;
                sets{iset}.label=strrep(sets{iset}.label,'./','');
                sets{iset}.label=strrep(sets{iset}.label,'.mat','');
                sets{iset}.label=strrep(sets{iset}.label,'coords_','');
                sets{iset}.pipeline=pipeline;
                opts_qpred_used{iset}=struct();
            case 2
                sets{iset}.type='qform';
                [d,sas{iset},opts_read_used{iset}]=psg_read_coorddata(data_fullname,setup_fullname,setfield(opts_read,'if_justsetup',1));
                nbtc=size(sas{iset}.btc_augcoords,2);
                opts_read.setup_fullname_def=opts_read_used{iset}.setup_fullname;
                %determine whether one of the strings in ray_minpts_default
                %is present in setup file name, and if so, use this to
                %determine the default for opts_rays
                opts_rays_use=psg_findray_setopts(opts_read_used{iset}.setup_fullname,opts_rays);
                %
                [rayss{iset},opts_rays_used{iset}]=psg_findrays(sas{iset}.btc_specoords,setfield(opts_rays_use,'permute_raynums',opts_read_used{iset}.permute_raynums));
                if opts_read.if_auto==0
                    if_aug_spe=getinp('1 to use augmented coords, 2 to use spec coords','d',[1 2],1);
                else
                    if_aug_spe=1;
                end
                switch if_aug_spe
                    case 1
                        btc_coords=sas{iset}.btc_augcoords;
                        aug_spe_string='aug';
                    case 2
                        btc_coords=sas{iset}.btc_specoords;
                        btc_coords(isnan(btc_coords))=0;
                        aug_spe_string='spec';
                end
                %get quadratic form
                if opts_read.if_auto==0
                    qform_source_type=getinp(' quadratic form choice: 1->use threshold data file, 2->use identity','d',[1 2],1);
                else
                    qform_source_type=1;
                end
                switch qform_source_type
                    case 1
                        if opts_read.if_auto==0
                            qform_datafile=getinp('data file with path','s',[0 1],opts_qpred.qform_datafile_def);
                        else
                            qform_datafile=opts_qpred.qform_datafile_def;
                        end
                        opts_qpred.qform_datafile_def=qform_datafile;
                        btc_thresh_data=getfield(load(qform_datafile),'r');
                        q=btc_thresh_data{opts_qpred.qform_modeltype}.results.qfit;
                        q_label=btc_thresh_data{opts_qpred.qform_modeltype}.setup.label;
                        disp(sprintf(' model %2.0f (%s) loaded from %s',...
                            opts_qpred.qform_modeltype,q_label,qform_datafile));
                        qform_source=qform_datafile;
                    case 2
                        q=eye(nbtc);
                        qform_source='[identity]';
                end
                %compute the model
                [d_qform{iset},d_mds{iset},opts_qpred_used{iset}]=psg_qformpred(q,btc_coords,rayss{iset},opts_qpred);
                if opts_read.if_auto==0
                    if_qform=getinp('1 to use qform, 2 for mds model (should be equivalent)','d',[1 2],1);
                else
                    if_qform=1;
                end
                switch if_qform
                    case  1
                        ds{iset}=d_qform{iset}; %use quadratic form model
                        qform_mds_string='qform';
                    case 2
                        ds{iset}=d_mds{iset};
                        qform_mds_string='mds';
                end
                sets{iset}.dim_list=1:length(d_qform{iset});
                sets{iset}.nstims=sas{iset}.nstims;
                sets{iset}.label_long=cat(2,opts_read_used{iset}.setup_fullname,' c:',aug_spe_string,' q:',qform_source,' m:',qform_mds_string);
                sets{iset}.label=sets{iset}.label_long;
                sets{iset}.label=strrep(sets{iset}.label,'../','');
                sets{iset}.label=strrep(sets{iset}.label,'./','');
                sets{iset}.label=strrep(sets{iset}.label,'stim/','');
                sets{iset}.label=strrep(sets{iset}.label,'btc_allraysfixedb_','');
                sets{iset}.label=strrep(sets{iset}.label,'100surrs_','');
                sets{iset}.label=strrep(sets{iset}.label,'.mat','');
                sets{iset}.pipeline=struct;
        end  %data or model
        if (iset==1)
            nstims=sas{iset}.nstims;
            typenames=sas{iset}.typenames;
        end
        opts_read_used{iset}.input_type=input_type_use;
        opts_read_used{iset}.input_type_desc=input_types{input_type_use};
        %check consistency of number of stimuli and typenames
        if sets{iset}.nstims~=nstims
            disp(sprintf('warning: expecting %3.0f stimuli, dataset %1.0f has %3.0f stimuli',nstims,iset,sets{iset}.nstims));
        else
            typenames_mismatch=[];
            for istim=1:nstims
                if ~strcmp(sas{iset}.typenames{istim},typenames{istim})
                    typenames_mismatch=[typenames_mismatch,istim];
                end
            end
            if ~isempty(typenames_mismatch)
                disp('warning: the following stimulus type names do not match:')
                for istimptr=1:length(typenames_mismatch)
                    istim=typenames_mismatch(istimptr);
                    disp(sprintf('   for stim %3.0f, expecting %20s found %20s',istim,typenames{istim},sas{iset}.typenames{istim}));
                end
            end
        end %number of stimuli match?
    end % isets
    %summarize and check
    disp(' ');
    disp('datasets selected:');
    for iset=1:nsets_pos
        disp(sprintf(' set %2.0f: dim range [%3.0f %3.0f] label: %s',iset,min(sets{iset}.dim_list),max(sets{iset}.dim_list),sets{iset}.label));
    end
    if opts_read.if_auto==0
        if_ok=getinp('1 if ok','d',[0 1]);
    else
        if_ok=1;
    end
end % if_ok
return

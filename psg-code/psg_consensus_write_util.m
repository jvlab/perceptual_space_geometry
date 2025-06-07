%psg_consensus_write_util
%utility to write consensus files
%
%  See also: PSG_ALIGN_STATS_DEMO, PSG_ALIGN_VARA_DEMO, PSG_WRITE_COORDDATA.
%
% 07Jun25: options added: embedded setup file, remove details from pipeline
%
for allow_scale=0:1
    ia=allow_scale+1;
    disp(sprintf(' calculations with allow_scale=%1.0f',allow_scale));
    if_write_consensus=getinp('1 to write a file with consensus (knitted) coordinate data and metadata  (-1 to embed setup metadata)','d',[-1 1]);
    if if_write_consensus~=0
         opts_write=struct;
        opts_write.data_fullname_def='[paradigm]pooled_coords_ID.mat';
        %
        sout_consensus=struct;
        sout_consensus.stim_labels=strvcat(sa_pooled.typenames);
        %
        opts=struct;
        opts.pcon_dim_max=pcon_dim_max; %maximum consensus dimension created   
        opts.pcon_dim_max_comp=pcon_dim_max; %maximum component dimension used
        opts.details=details(:,ia); %details of Procrustes alignment
        opts.opts_read_used=opts_read_used; %file-reading options
        opts.opts_align_used=opts_align_used; %alignment options
        opts.opts_pcon_used=opts_pcon_used(:,ia); %options for consensus calculation for each dataset
        if if_write_consensus==-1
            sout_consensus.setup=sa_pooled;
        end
        sout_consensus.pipeline=psg_coord_pipe_util('consensus',opts,sets);
        if getinp('1 to remove details from pipeline to shorten output file','d',[0 1])
            sout_consensus.pipeline.opts=rmfield(sout_consensus.pipeline.opts,'details');
        end
        opts_write_used=psg_write_coorddata([],ds_knitted{ia},sout_consensus,opts_write);
        %
        metadata_fullname_def=opts_write_used.data_fullname;
        metadata_fullname_def=metadata_fullname_def(1:-1+min(strfind(cat(2,metadata_fullname_def,'_coords'),'_coords')));
        if isfield(sa_pooled,'nsubsamp')
            metadata_fullname_def=cat(2,metadata_fullname_def,sprintf('%1.0f',sa_pooled.nsubsamp));
        end
        if getinp('1 to write metadata (setup file) for consensus dataset','d',[0 1])
            s=sa_pooled;
            metadata_fullname=getinp('metadata file name','s',[],metadata_fullname_def);
            save(metadata_fullname,'s');
        end
    end
end

%psg_consensus_write_util
%utility to write consensus files
%
%  See also: PSG_ALIGN_STATS_DEMO, PSG_ALIGN_VARA_DEMO.
%
for allow_scale=0:1
    ia=allow_scale+1;
    disp(sprintf(' calculations with allow_scale=%1.0f',allow_scale));
    if getinp('1 to write a file with consensus (knitted) coordinate data and metadata','d',[0 1])
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
        sout_consensus.pipeline=psg_coord_pipe_util('consensus',opts,sets);
        opts_write_used=psg_write_coorddata([],ds_knitted{ia},sout_consensus,opts_write);
        %
        metadata_fullname_def=opts_write_used.data_fullname;
        metadata_fullname_def=metadata_fullname_def(1:-1+min(strfind(cat(2,metadata_fullname_def,'_coords'),'_coords')));
        if isfield(sa_pooled,'nsubsamp')
            metadata_fullname_def=cat(2,metadata_fullname_def,sprintf('%1.0f',sa_pooled.nsubsamp));
        end
        metadata_fullname=getinp('metadata file name','s',[],metadata_fullname_def);
        s=sa_pooled;
        save(metadata_fullname,'s');
    end
end

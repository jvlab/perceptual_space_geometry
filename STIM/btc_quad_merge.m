%btc_quad_merge: merge one or more quad map resources into a base file
%
% The number of components in files to be merged must be equal to the number in the base file, or to 1
% The number of stimuli in files to be merged cannot exceed the number in the base file
%
% r.which_merged,r.mix_scenarios_[merged,files] keeps track of what is merged in, r.mix_scneario_name is updated
%
%
%   See also:  PSG_SPOKES_SETUP, SPOKES_SETUP_CREATE, BTC_QUAD_BGCA_MAKE, BTc_QUAD_BGCA_MAKE_ONEC.
%
if ~exist('res_file_def') res_file_def='../stim/btc_quad_bgca_make_t20_sc7_it25000.mat'; end
if ~exist('res_file_merge_def') res_file_merge_def='../stim/btc_quad_bgca_make_onec_t20_it5000'; end
res_file_base=getinp('base resource file','s',[],res_file_def);
r=getfield(load(res_file_base),'r');
nstims=length(r.spec_final);
disp(sprintf('base resource file loaded with %2.0f stimuli and %2.0f components.',nstims,r.ncomponents));
disp(r)
r.which_merged=zeros(1,nstims);
r_base=r;
if_done=0;
nmerge=0;
res_file_merges=cell(0);
%
must_match={'map_size_final','alpharange_tol'};
%one-dim cells
merge_1dcell={'spec_final','btc_specoords','btc_augcoords','filebase','maps_init','maps_final','opts_metro_used'};
merge_d3nc={'stats_init_ideal','stats_init','stats_final'}; %dim 3 is the number of components
merge_2c={'stats_final_ideal','stats_final_avg','stats_final_max_component_diff','stats_final_max_dev_ideal','stats_final_avg_dev_ideal'};
merge_scal2vec={'bc_target','ncomponents'}; %can change from scalar to vector
%
ncomponents_orig=r.ncomponents;
%
while (if_done==0)
    nmerge=nmerge+1;
    if_ok=0;
    while (if_ok==0)
        disp('Note: the number of components in files to be merged must be equal to the number in the base file, or to 1'); 
        disp('Note: the number of stimuli in files to be merged cannot exceed the number in the base file.');
        res_file_merges{nmerge}=getinp(sprintf('resource file %1.0f to merge',nmerge),'s',[],res_file_merge_def);
        r_merge=getfield(load(res_file_merges{nmerge}),'r');
        nstims_merge=length(r_merge.spec_final);
        disp(sprintf('to-be-merged file %1.0f loaded with %2.0f stimuli and %2.0f components',nmerge,nstims_merge,r_merge.ncomponents));
        disp(r_merge);
        if_ok=1;
        if (r_merge.ncomponents~=ncomponents_orig) & (r_merge.ncomponents~=1)
            if_ok=0;
            disp('bad number of components');
        end
        if nstims_merge>nstims
            if_ok=0;
            disp('bad number of stimuli');
        end
        for ifn=1:length(must_match)
            fn=must_match{ifn};
            if isfield(r,fn) & isfield(r_merge,fn)
                if r.(fn)~=r_merge.(fn)
                    if_ok=0;
                    disp('mismatch of %s',fn)
                end
            end
        end
        name_string=[];
        if (if_ok) %ok to merge
            stim_list=getinp('stimulus numbers to merge','d',[1 nstims_merge]);
            for iptr=1:length(stim_list)
                istim=stim_list(iptr);
                name_string=cat(2,name_string,' ',sprintf('%1.0f',istim));
                %1-d cells
                for ifn=1:length(merge_1dcell)
                    fn=merge_1dcell{ifn};
                    if ~isfield(r,fn) | ~isfield(r_merge,fn)
                        disp(sprintf('field %s missing from base file or to-be-merged file, ignoring',fn))
                    else
                        r.(fn){istim}=r_merge.(fn){istim};
                    end
                end %1-d cells
                %arrays in which third dim is number of components
                for ifn=1:length(merge_d3nc)
                    fn=merge_d3nc{ifn};
                    if ~isfield(r,fn) | ~isfield(r_merge,fn)
                        disp(sprintf('field %s missing from base file or to-be-merged file, ignoring',fn))
                    else
                        x=r_merge.(fn)(istim,:,:);
                        if size(x,3)<ncomponents_orig
                            x(:,:,size(x,3)+1:ncomponents_orig)=NaN;
                        end
                        r.(fn)(istim,:,:)=x;
                    end
                end %d3nc
                %2-d arrays (stim x coords)
                for ifn=1:length(merge_2c)
                    fn=merge_2c{ifn};
                    if ~isfield(r,fn) | ~isfield(r_merge,fn)
                        disp(sprintf('field %s missing from base file or to-be-merged file, ignoring',fn))
                    else
                        r.(fn)(istim,:)=r_merge.(fn)(istim,:);
                    end
                end %2c
                %scalars that may become vectors
                for ifn=1:length(merge_scal2vec)
                    fn=merge_scal2vec{ifn};
                    if ~isfield(r,fn) | ~isfield(r_merge,fn)
                        disp(sprintf('field %s missing from base file or to-be-merged file, ignoring',fn))
                    else
                        if length(r.(fn))>1 %if it's already a vector, then just merge in
                            r.(fn)(istim)=r_merge.(fn);
                        else
                            if r_merge.(fn)~=r.(fn) %new value, create a vector; otherwise nothing to do
                                r.(fn)=repmat(r.(fn),1,nstims);
                                r.(fn)(istim)=r_merge.(fn);
                            end
                        end
                    end
                end %2c
            end %stim_list
            r.which_merged(stim_list)=nmerge;
            r.mix_scenarios_merged{nmerge}=r_merge.mix_scenario;
            r.mix_scenarios_files{nmerge}=res_file_merges{nmerge};
            r.mix_scenario_name=cat(2,r.mix_scenario_name,sprintf(' +merge[%1.0f]: stims %s',nmerge,name_string));
            disp('merge completed');
            disp(r);
        end %ok to merge
    end
    if_done=getinp('1 if done','d',[0 1]);
end
if getinp('1 to write a resource file','d',[0 1],1);
    fn_write=getinp('file name to write, e.g., btc_quad_bgca_merged*','s',[]);
    save(fn_write,'r');
end

% psg_geomodels_task_summ: special-purpose script to summarize analysis of
% nested geometric models, for task ms
%
%  See also:  PSG_GEOMODELS_RUN, PSG_GEOMODELS_SUMM.

%
% each of these files is created by psg_geomodels_run, from consensus
% spaces across the subjects common to the tasks and coordinates compared
%
if ~exist('file_list')
    file_list={...
    '../iso/task_ms_geomodels_thresh2sim_bc6.mat',...
    '../iso/task_ms_geomodels_thresh2sim_bcpm3.mat',...
    '../iso/task_ms_geomodels_thresh2sim_bgca3.mat',...
    '../iso/task_ms_geomodels_thresh2sim_dgea3.mat',...
    '../iso/task_ms_geomodels_sim2br_bgca3.mat',...
    '../iso/task_ms_geomodels_sim2br_dgea3.mat',...
    '../iso/task_ms_geomodels_sim2gm_bgca3.mat',...
    '../iso/task_ms_geomodels_sim2gm_dgea3.mat',...
    '../iso/task_ms_geomodels_sim2wm1000_bgca3.mat',...
    '../iso/task_ms_geomodels_sim2wm1000_dgea3.mat'};
end
if ~exist('dim_list')
    dim_list=[2:5];
end
%first model is the simpler (nested) model, second model is the more complex one
model_pairs={...
    'procrustes_noscale','procrustes_scale';...
    'procrustes_scale','affine_offset';...
    'affine_offset','projective';...
    'affine_offset','pwaffine'};
npairs=size(model_pairs,1);
%
nfiles=length(file_list);
for ifile=1:nfiles
    fn=file_list{ifile};
    underscores=find(fn=='_');
    tasks=fn(underscores(3)+1:underscores(4)-1);
    coords=fn(underscores(4)+1:end-4);
    disp(sprintf(' file %40s: tasks: %12s, coordinates: %6s',fn,tasks,coords));
    load(fn);
    disp('loaded');
    for ipair=1:npairs
        disp(sprintf('%20s vs %20s',model_pairs{ipair,2},model_pairs{ipair,1}));
        for idim=dim_list
            model_types=results{idim,idim}.model_types_def.model_types;
            model_indices=zeros(npairs,2);
            for im=1:2
                mi=strmatch(model_pairs{ipair,im},model_types,'exact');
                if length(mi)==1
                    model_indices(ipair,im)=mi;
                else
                    disp(sprintf(' %s not found in dim %3.0f.',model_pairs{ipair,im},idim));
                end
            end
            if all(model_indices(ipair,:)>0)
                 nshuff=results{idim,idim}.nshuff;
                d_shuff=reshape(results{idim,idim}.d_shuff(model_indices(ipair,2),:,model_indices(ipair,1),:),[nshuff,2]);
                d=results{idim,idim}.d(model_indices(ipair,2));
                d_nested=results{idim,idim}.d(model_indices(ipair,1));
                d_shuff_mean=mean(d_shuff);
                d_shuff_count=squeeze(results{idim,idim}.surrogate_count(model_indices(ipair,2),model_indices(ipair,1),:));
                disp(sprintf(' dim %2.0f: d: %7.5f  nested d: %7.5f; d_surr (den_surr, den_orig): mean %7.5f %7.5f count %3.0f %3.0f',...
                    idim,d,d_nested,d_shuff_mean,d_shuff_count));
            end
        end %idim
    end %ipair
end


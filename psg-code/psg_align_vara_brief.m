%brief script to summarize psg_algin_vara_demo
%
%Run after psg_align_vara_demo.
%
% See also:  PSG_ALIGN_VARA_DEMO.
%
if ~exist('quantiles_summ') quantiles_summ=[.001 0.01 .05 .1 .25 .5 .75 .95]; end
%
for iset=1:results.nsets
    disp(sprintf('set %2.0f: %s',iset,results.sets{iset}.label));
end
disp('nsets')
disp(results.nsets);
disp('nshuffs')
disp(results.nshuffs);
disp('group list');
disp(results.gp_list);
disp('if_normscale')
disp(results.if_normscale);
disp('tags');
disp(results.opts_multi.tags)
for if_scale=0:1
    disp(sprintf('summary, allow_scale=%1.0f',if_scale));
    rmsdev_grpwise=reshape(results.rmsdev_grpwise(:,1,1+if_scale),[results.dim_max,1])';
    rmsdev_grpwise_shuff=reshape(results.rmsdev_grpwise_shuff(:,1,1+if_scale,1,:),[results.dim_max,results.nshuffs])';
    pvals=sum(double(rmsdev_grpwise_shuff<=rmsdev_grpwise))/results.nshuffs;
    equality_count=sum(double(rmsdev_grpwise_shuff==rmsdev_grpwise));
    disp('data, p-values, quantiles, equality_count')
    disp('        dimension')
    disp([[NaN 1:results.dim_max];[NaN rmsdev_grpwise];[NaN pvals];[quantiles_summ(:),quantile(rmsdev_grpwise_shuff,quantiles_summ)];[NaN equality_count]]);
end
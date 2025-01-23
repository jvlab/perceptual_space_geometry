%brief script to summarize psg_algin_vara_demoopts_pcon
%
%Run after psg_align_vara_demo.
%
% See also:  PSG_ALIGN_VARA_DEMO.
%
disp(' ');
if ~exist('quantiles_summ') quantiles_summ=[.001 0.01 .05 .1 .25 .5 .75 .95]; end
%
disp(sprintf('nsets: %4.0f, nshuffs: %5.0f, if_normscale: %1.0f',results.nsets,results.nshuffs,results.if_normscale));
for iset=1:results.nsets
    disp(sprintf('set %2.0f: %s',iset,results.sets{iset}.label));
end
disp('group list');
disp(results.gp_list);
disp('tags');
disp(results.opts_multi.tags)
for if_scaling=0:1
    disp(sprintf('summary, allow_scale=%1.0f',if_scaling));
    rmsdev_grpwise=reshape(results.rmsdev_grpwise(:,1,1+if_scaling),[results.dim_max,1])';
    rmsdev_grpwise_shuff=reshape(results.rmsdev_grpwise_shuff(:,1,1+if_scaling,1,:),[results.dim_max,results.nshuffs])';
    pvals=sum(double(rmsdev_grpwise_shuff<=rmsdev_grpwise))/results.nshuffs;
    equality_count=sum(double(rmsdev_grpwise_shuff==rmsdev_grpwise));
    disp('data, p-values, quantiles, equality_count')
    disp('        dimension')
    disp([[NaN 1:results.dim_max];[NaN rmsdev_grpwise];[NaN pvals];[quantiles_summ(:),quantile(rmsdev_grpwise_shuff,quantiles_summ)];[NaN equality_count]]);
end
function map_new=mlis_run_restoremv(map,desired_mean,desired_var)
% map_new=mlis_run_restoremv(map,desired_mean,desired_var) is a utility function for mlis_run and family
% to rescale a map to give a desired mean and variance
%
%  See also:  MLIS_RUN, MLIS_RUN_SUB,  MLIS_WHITENING_FILTER.
%
map=map-mean(map(:));
map_new=desired_mean+map*sqrt(desired_var/mean(map(:).^2));
return


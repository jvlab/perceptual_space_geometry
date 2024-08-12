function [results,opts_used]=psg_majaxes(d_ref,sa_ref,d_adj,sa_adj,results_geo,opts)
% [results,opts_used]=psg_majaxes(d_ref,sa_ref,d_adj,sa_adj,results_geo,opts) analyzes an affine geometric model to determine major axes
%
%  d_ref, sa_ref: data and sa (labeling) structure for reference dataset, typically from psg_read_coorddata
%  d_adj, sa_adj: data and sa (labeling) structure for adjustable dataset, typically from psg_read_coorddata
%  results_geo: results field from psg_geomodels_run
%  opts: options (can be empty)
%
%  results: analysis results
%  opts_used: options used
%
%   See also: PSG_GEOMODELS_RUN, PSG_GEOMODELS_DEFINE.
%
if nargin<=5 opts=struct; end
model_types_def=psg_geomodels_define();
opts.model_types_def=model_types_def;
results=struct;
%
opts_used=opts;
return

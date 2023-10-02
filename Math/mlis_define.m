function mlis_opts=mlis_define(mlis_opts_def)
% mlis_opts=mlis_define(mlis_opts_def) sets up the options for routine to
% modify the local statistics in a map by jittering the phases
%
% mlis_opts_def: supplied options, may be empty
% mlis_opts: full options structure
%  ngens_max: maximum number of generations
%  phajit_max: maximum phase jitter, as fraction of cycle (1=full jitter)
%
% 03Mar20: added option (btc_down) to allow target statistics to be calculated from subsampled checks
% 03Mar20: added options (phmod_sfmin_nyq, phmod_sfmax_nyq) to restrict range of frequencies that are modified
% 05May20: added options (phmod_avoid_onaxis, phmod_avoid_nyq) to avoid jittering phases of on-axis (pure horiz or pure vert) frequencies
% 05May20: added options (if_mirror,if_mirror_btc) related to mirror-flipping to reduce spectral leakage
% 05May20: added distance_success: distance (in btc units) to declare success
% 14Sep20: added options for regions: portions within image to make local modifications
% 15Oct20: added option for regions: reg_align_btc_down to force regions to align with downsampling for btc statistics
% 03May21: added documentation for binarize
%
%   See also:  MLIS_TEST, MLIS_RUN, MLIS_SUB_RUN, FFDM_BTC_CALC, MLIS_ALGLIB_PILOT, MLIS_WARN, MLIS_REGS_ASK.
%
if (nargin<1)
    mlis_opts_def=[];
end
mlis_opts=mlis_opts_def;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters relevant to search algorithm and phase substition algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
mlis_opts=filldefault(mlis_opts,'target_mode','absolute');
% 'absolute': supplied target coords are the goal
% 'relative': supplied target coords are interpreted as an offset relative to starting point
mlis_opts=filldefault(mlis_opts,'distance_exponent',2); %how to compute distance, Inf means maximum
mlis_opts=filldefault(mlis_opts,'distance_success',0.05); %distance in btc coords to declare success
mlis_opts=filldefault(mlis_opts,'phmod_sfmin_nyq',0); %so that all SF's are included by default
mlis_opts=filldefault(mlis_opts,'phmod_sfmax_nyq',2); %so that all SF's are included by default
mlis_opts=filldefault(mlis_opts,'phmod_avoid_onaxis',0); %n avoids all frequencies that are less than n freq bins from an axis (n=0: no effect)
mlis_opts=filldefault(mlis_opts,'phmod_avoid_nyq',0); %n avoids all frequencies that are less than n freq bins from an axis + Nyquist
%
mlis_opts=filldefault(mlis_opts,'if_mirror',0); %1 to apply phase algorithms to an image that has been preprocessed by mirroring in each axis
mlis_opts=filldefault(mlis_opts,'if_mirror_btc',0); %1 to also mirror-flip the binary texture image, ignored if if_mirror=0
%
%statistics calculation options
%
mlis_opts=filldefault(mlis_opts,'btc_whitened',1); %modify statistics based on whitened map
mlis_opts=filldefault(mlis_opts,'binarize','median'); %'median','mean', 'hard' (cut at 0.5), 'none',...
%                                                      'median_nonan' [ignores NaN's], 'mean_nonan' [ignores NaN's], or a value for the quantile
mlis_opts=filldefault(mlis_opts,'btc_down',1); % downsampling for computation of btc statistics
%
mlis_opts=filldefault(mlis_opts,'if_cosbell',0); %windowing: 0 for none, 1 or 2 for cosine bell window (1: den=R-1, 2: den=R)
mlis_opts=filldefault(mlis_opts,'if_window_all',0); %0 to just use windowed map for whitening, 1 to use windowed map for phases (irrelevant if if_cosbell=0)
mlis_opts=filldefault(mlis_opts,'padfactor',1); %padding factor
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters relevant to search algorithm (mlis_run)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
mlis_opts=filldefault(mlis_opts,'niters_max',10); %maximum number of iterations
mlis_opts=filldefault(mlis_opts,'worse_maxp',0.1); %max probabily of using a new phase jitter if it worsens the cost function
mlis_opts=filldefault(mlis_opts,'worse_maxd',0.02); %max increase in distance to  allow 
% -- probability of backtracking shrinks linearly to zero as this is reached
mlis_opts=filldefault(mlis_opts,'save_select',2);
mlis_opts=filldefault(mlis_opts,'save_select_strings',{'only first','all','changes','improvement','best yet'});
%-1: only save initial step
% 0: save all
% 1: save if change
% 2: save if improvement over prev
% 3: save if best yet
mlis_opts.save_select_string=mlis_opts.save_select_strings{2+mlis_opts.save_select};
%
%phase jitter parameters
%
mlis_opts=filldefault(mlis_opts,'phjit_max',0.1); %maximum underlying phase step
mlis_opts=filldefault(mlis_opts,'phcum_max',0.3); %maximum cumulative phase change
mlis_opts=filldefault(mlis_opts,'phcum_lim','hard'); %options:
mlis_opts=filldefault(mlis_opts,'phcum_lim_opts',{'none','hard','rms','tanh'});
%   'none': cumulative phase shift is unbounded, phase shift used = cumulative shift
%   'hard': cumulative phase shift set to max if it exceeds phcum_max, phase shift used= cumulative phase shift
%   'rms' : rms cumulative phase shift cannot exceed max (all phases reduced proportionally if it does), phase shift used=cumulative phase shift
%   'tanh': cumulative phase shift is unbounded,  phase shift used = phcum_max*tanh(phcum/phcum_max)
mlis_opts=filldefault(mlis_opts,'phjit_num_max',Inf); %maximum number of phases that can be jittered
mlis_opts=filldefault(mlis_opts,'phjit_frac_max',1); %maximum fraction of phases that can be jittered
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters relevant to phase substition algorithm (mlis_run_sub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
mlis_opts=filldefault(mlis_opts,'phsub_vals',setdiff([-1:0.1:1],0)); %list of values of the  btc statistic to try
mlis_opts=filldefault(mlis_opts,'phsub_wts',[0.25:0.25:1]); %list of weights (1: use only the btc phase, 0: use only the original image weight
mlis_opts=filldefault(mlis_opts,'phsub_nex',4); %number of examples of each btc value and weight 
mlis_opts=filldefault(mlis_opts,'phsub_mix','uvec'); %how to mix phase of original image and btc phase
mlis_opts=filldefault(mlis_opts,'phsub_mix_opts',{'phase','uvec'});
%   'phase': blends as  linear combination of phases
%   'uvec': blends as  linear combination of unit vectors
mlis_opts=filldefault(mlis_opts,'phsub_nkeep',1); %number of value-weight combinations to keep
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters relevant to regional modification: names should all begin with reg_
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
mlis_opts=filldefault(mlis_opts,'reg_count',0); %0 to disable, 1 or more to enable
mlis_opts=filldefault(mlis_opts,'reg_size_mode','abs'); %how to specify region size: abs->in pixels, rel->as fraction of image
mlis_opts=filldefault(mlis_opts,'reg_size_mode_opts',{'abs','rel'});
mlis_opts=filldefault(mlis_opts,'reg_size',32); %region size, either in pixels or as fraction of image size
mlis_opts=filldefault(mlis_opts,'reg_placement_mode','grid'); %method of placing regions:  grid (with random jitter), or repelling
mlis_opts=filldefault(mlis_opts,'reg_placement_mode_opts',{'grid','repel'}); %placement mode options
mlis_opts=filldefault(mlis_opts,'reg_grid_arrange',0); %arrangement regularity in grid mode, in [0 1]:  0 is Poisson, 1 is maximally regular
mlis_opts=filldefault(mlis_opts,'reg_repel_exponent',0); %strength of repulsion (0: none, Inf: maximal)
mlis_opts=filldefault(mlis_opts,'reg_repel_tol',0.0001); %tolerance for ignoring diffrences in peak probabilitie
mlis_opts=filldefault(mlis_opts,'reg_add_weight',1.0); %weight fraction for adding regions
mlis_opts=filldefault(mlis_opts,'reg_add_mode','seq'); %how to add the regions: seq->sequentially, sim->all at once (by averaging)
mlis_opts=filldefault(mlis_opts,'reg_add_mode_opts',{'seq','sim'});
mlis_opts=filldefault(mlis_opts,'reg_blend_mode','hard'); %how to blend the regions:  hard, cos bell (den R-1), gaussian 2 SD
mlis_opts=filldefault(mlis_opts,'reg_blend_mode_opts',{'hard','cb1','cb2','gau'});
mlis_opts=filldefault(mlis_opts,'reg_align_btc_down',0); %1 to force region positions to align with the downsampling for calclation of image statistics (irrelevant if btc_down=1)
mlis_opts=filldefault(mlis_opts,'reg_verbose',0); %log level for region algorithms
return

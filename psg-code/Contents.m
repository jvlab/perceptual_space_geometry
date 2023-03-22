% Demos, etc for gr23 by J. Victor
%
% Texture demos
%   btc_extremes_demo: demonstrate extremes of theta-pair textures
%   spokes_layout_demo: demonstrate btc stimuli along spokes
%   spokes_setup_create: make a library of spoke specifications
%
% Perceptual space geometry experiments: setup
%   psg_cond_create: create a cell array of image names for a cond file
%   psg_cond_write: write a condition file
%   psg_defopts: set up default options
%   psg_sessconfig_make: create a set of sessions and trials
%   psg_session_stats: tally the statistics of a session (counts of how many trials each stimulus is used, etc)
%   psg_sess_perm: apply a permutation to stimulus numbers to create a new session
%   psg_setup_demo: demonstrate setup of a perceptual space geometry experiment
%   psg_showtrial: show a trial, from the cond (*.csv) and image files
%   psg_spec2filename: convert spec structure to file name
%   psg_spokes_setup: create stimuli and cond files for perceptual space geometry expt with btc stimuli on spokes
%   text_csvwrite: write a text csv file (e.g., for configuration file)
%
% Perceptual space geometry experiments: geometric analysis
%   psg_colors_legacy: get legacy colors
%   psg_consensus_demo: demonstrate Procrustes consensus and plotting
%   psg_coords_fix: fix extra stim_labels entries in a mat-file
%   psg_findrays: parse a set of stimulus coordinates into rays
%   psg_get_coordsets: read coordinates from psychophysical data or quadratic form model
%   psg_pcaoffset: pca after offset, and reconstruction by successive dimensions
%   psg_planecycle: analyze and order points in a plane
%   psg_plotangles: plot angles between rays
%   psg_plotcoords: plot psg coordinates
%   psg_procrustes_demo: compare multiple datasets via Procrustes method
%   psg_qformpred: predict perceptual space coords from quadratic form model of thresholds
%   psg_qformpred_demo: demonstrate predictions from quadratic form model of thresholds
%   psg_rayangles: compute angles between rays
%   psg_rayfit:  fit a coordinate structure to rays
%   psg_read_coorddata: read coordinates data inferred from a psg experiment
%   psg_typenames2colors: assign colors to array types, for plotting
%   psg_typenames2colors_test: test psg_typenames2colors
%   psg_visualize: plot several pages of visualization of psg coords
%   psg_visualize_demo: demonstrate basic visualization of psg coords
%
% Perceptual space geometry experiments: choice probability analysis
%   psg_ineq_apply: apply the inequality conditions of psg_ineq_logic to a eet of obsrvations, and do flips
%   psg_ineq_edgecount: counts the edges in output of psg_ineq_logic
%   psg_ineq_logic: sets up logic for excluded rank-choice-probabilities, for tests of symmetry, umi, addtree, etc.
%   psg_ineq_logic_demo: test psg_ineq_logic
%   psg_permutes_logic: sets up permutations for tests of symmetry, umi, addtree
%   psg_probs_check: compares versions of psg_umi_triplike
%   psg_quad_stats: calculate and display statistics of quadruplets relevant to testing additivity and addtree
%   psg_quad_stats_demo: demonstrate psg_quad_stats
%   psg_read_choicedata: read choice probability data from a psg experiment
%   psg_stats_tally: utility for to tally statistics, for psg_triad_stats, psg_umi_stats
%   psg_tent_stats: calculate and display tent statistics (condition for addtree)
%   psg_tent_stats_demo: demonstrate psg_tent_stats
%   psg_triad_stats: calculate and display triad statistics
%   psg_triplet_choices: extract triplets of choice probabilities from choices data file
%   psg_umi_stats: calculate and display statistics of trials relevant to testing ultrametric inequality
%   psg_umi_stats_demo: demontrate psg_umi_stats
%   psg_umi_triplike: analyze likelihoods that rank choice probabilities are consistent with ultrametric inequality and symmetry
%   psg_umi_triplike_demo: apply psg_umi_triplike to data
%   psg_umi_triplike_plot: plot detailed results of psg_umi_triplike_demo
%   psg_umi_triplike_plota: plot summary (asymptotic0 results of psg_umi_triplike_demo

%   Copyright (c) 2022, 2023 by J. Victor
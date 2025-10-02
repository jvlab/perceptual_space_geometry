% Demos, pilots, etc for projective space geometry by J. Victor
%
% Texture demos
%   btc_extremes_demo: demonstrate extremes of theta-pair textures
%   spokes_layout_demo: demonstrate btc stimuli along spokes
%   spokes_setup_create: make a library of spoke specifications
%
% Perceptual space geometry experiments: setups
%   faces_mpi_inventory: inventory of the contents of the MPI faces database
%   faces_mpi_inventory_demo: use faces_mpi_inventory
%   faces_mpi_get_setups: specify several face paradigm setups for MPI faces database
%   faces_mpi_setup_create: use faces_mpi_get_setups to select and structure face files
%   faces_mpi_psg_setup: generate the control files for a psg experiment with MPI faces
%   irgb_btcstats: compute binary texture stats at multiple scales from a grayscale image
%   irgb_gzranking_analyze: analyze Giesel-Zaidi ranking data
%   irgb_gzranking_getimages: one-time script to extract rgb values of material images
%   irgb_gzranking_read: read Giesel-Zaidi ranking data
%   irgb_modify: phase-scramble and whiten images
%   irgb_modify_demo: demonstrate irgb_modify on Giesel/Zaidi textiles
%   irgb_modify_dict: dictionary of image modifications
%   irgb_spec_defaults: set up default specifications for iid rgb stimuli
%   irgb_spec_make: make specifications for iid rgb stimuli
%   irgb_spec_modify: modif specifications for iid rgb stimuli
%   irgb_stim_make: make stimuli from specifications iid rgb stimuli
%   irgb_psg_imgs_setup: generate stimulus files files for a psg experiment with manipulated rgb images
%   irgb_psg_sess_setup: generate only control files for a psg experiment with manipulated rgb images
%   irgb_psg_setup: generate control files for a psg experiment with iid rgb stimuli
%   psg_cond_create: create a cell array of image names for a cond file
%   psg_cond_sess_split: split an existing condition file into sub-sessions
%   psg_cond_write: write a condition file
%   psg_defopts: set up default options
%   psg_saw_anon2id: rename and copy SAW data files from JN paper to de-anonymize subjects
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
%   find_psg_xform_test: test fitting projective transformations via persp_xform_find
%   psg_align_coordsets: align coordinate datasets with partially overlapping stimuli
%   psg_align_knit_demo: test psg_align_coordsets, psg_remnan_coordsets, procrustes_consensus with partial overlaps
%   psg_align_stats_demo: align coordinate datasets and find consensus, with partial overlaps
%   psg_align_vara_brief: utility to summarize selected outputs
%   psg_align_vara_demo: analyze variance explained in global and grouped consensus datasets
%   psg_align_vara_plot: plotting script for psg_align_vara_demo
%   psg_align_vara_task: customized version of psg_align_vara_demo for bulk analysis of task/texture experiments
%   psg_align_vara_util: plotting utility for psg_align_vara_demo
%   psg_btcmeta_symapply: apply symmetry operation to btc metadata
%   psg_btcremz: simplify coordinate descritors when a coordinate is set to zero
%   psg_btcsyms: determine symmetry operations on btc coordinates
%   psg_cent_recip:  compute centrality and reciprocity indices (Tversky and Hutchinson)
%   psg_cent_recip_demo:  demonstrate psg_cent_recip
%   psg_choices_fix: fix extra stimlist entries in a mat-file
%   psg_choicedata_fixmeta: fix metadata (column headers) in test files 
%   psg_colors_legacy: get legacy colors
%   psg_color2rgb: convert color letters to rgb triplets
%   psg_commonstims: restrict coords and metadata to stimuli in common to one or more datasets
%   psg_consensus_demo: demonstrate Procrustes consensus and plotting
%   psg_consensus_write_util: utility to write consensus data and metadata
%   psg_coord_pipe_demo: pipelining of transformations, including consensus and pca rotation
%   psg_coord_pipe_util: utility to create pipeline structure
%   psg_coorddata_parsename: determine domain type ('btc', etc), short file name, default setup file from full coordinate file name
%   psg_coords_fillin: fill in missing lower-dimensional models
%   psg_coords_fix: fix extra stim_labels entries in a mat-file
%   psg_curv_loglik_plot: plot log likelihoods from SAW's curvature analysis csv files
%   psg_dbclean_grouping: special-purpose script to delete grouping experiments with 10 sessions from psg_[raystats|llfits]_summ databases
%   psg_dimgrps_util: utilty for specifying nesting dimensions
%   psg_dimstat_task: statistics for dimensionality of perceptual space, task project
%   psg_dist_heatmaps_demo: create distance heatmaps from a coordinate file or model
%   psg_findray_setopts: set options for psg_findrays
%   psg_findrays: parse a set of stimulus coordinates into rays
%   psg_geomodels_apply: apply a geometric transformation of a generic model
%   psg_geomodels_apply_test: tests consistency of the suite of geometric transformations
%   psg_geomodels_apply_util: utility for psg_geomodels_apply_test
%   psg_geomodels_define: define geometric models, standard options, and nesting relationships
%   psg_geomodels_illus: illustrate geometric models with toy data
%   psg_geomodels_ndof: number of degrees of freedom in a geometric model
%   psg_geomodels_run: like pst_geomodels_test but across multiple dimensions and saves results
%   psg_geomodels_summ: summarize results from psg_geomodels_run
%   psg_geomodels_test: test multiple geometric models
%   psg_geo_affine: fit an affine model, with and without offset
%   psg_geo_general: general wrapper for fitting geometric models with standardized calling conventions
%   psg_geo_layouts_setup: set up example coordinates for simulating geometrical transformations
%   psg_geo_procrustes: fit procrustes model (rotation and translation), with and without scaling
%   psg_geo_projective: fit projective model
%   psg_geo_pw_va_test: tests psg_geo_pwaffine_va and psg_geo_pwprojective_va, also checking fitted transform
%   psg_geo_pwaffine: piecewise affine model fitting
%   psg_geo_pwaffine_va: optimize a piecewise affine model with known cutplanes
%   psg_geo_pwaffine_va_test: test psg_geo_pwaffine_va for a range of dimensions and number of cuts, focusing on the coord change
%   psg_geo_pwaffine_obj: objective function for psg_geo_pwaffine
%   psg_geo_pwaffine_test: test options for psg_geo_pwaffine
%   psg_geo_pwprojective_pset: set up projection params to ensure piecewise continuity
%   psg_geo_pwprojective_va: optimize a piecewise projective model with known cutplanes and projection params
%   psg_geo_transforms_setup: set up example illustrative geometrical transformations
%   psg_get_coordsets: read coordinates from psychophysical data or quadratic form model
%   psg_get_geotransforms: get a geometric transformation from a file
%   psg_get_transform: specify a coordinate transformation (from keyboard)
%   psg_getgps: utility to get groups of files
%   psg_isomap: isomap embedding from distances
%   psg_isomap_demo: demonstrate psg_isomap and plot
%   psg_legend_keep: select graphic targets of a legend
%   psg_legend_util: relabel the legends of a plot made by psg_visualize, psg_plotcoords
%   psg_llfits_dbplot: plots data from databases created by psg_llfits_summ
%   psg_llfits_summ: create a database of log likelihood params from coordinate files
%   psg_lljit: compute log likelihoods of a model with jittered coordinates
%   psg_lljit_crit: use psg_lljit to compute critical jitters to reduce log likeklhood by a given p-value
%   psg_lljit_crit_demo: test psg_lljit_crit, psg_lljit
%   psg_localopts: local options (file names, etc.) for psg coords and choices
%   psg_majaxes: analyze axes of an affine transformation
%   psg_majaxes_reorder: utility for psg_majaxes
%   psg_meta_restrict: restrict a metadata structure to specific stimuli
%   psg_natcoords: determine natural coordinates from empirical axis trajectories in a perceptual space
%   psg_natcoords_demo: demonstrate psg_natcoords
%   psg_noneuc_dbplot: plots data from databases created by psg_llfits_summ
%   psg_noneuc_summ: create a database of log likelihood params from coordinate files
%   psg_parse_filename: parse a file name to determine paradigm type, paradigm, subject id, file type
%   psg_pcaoffset: pca after offset, and reconstruction by successive dimensions
%   psg_planecycle: analyze and order points in a plane
%   psg_plotangles: plot angles between rays
%   psg_plotcoords: plot psg coordinates
%   psg_plotmults: plot multipliers (gains) on each ray
%   psg_procrustes_demo: compare multiple datasets via Procrustes method
%   psg_procrustes_regr_demo: compare Procrustes, regression,and projective transforms
%   psg_procrustes_regr_test: test comparisons of Procrustes, regression, and projective transforms
%   psg_procrustes_task: procrustes analysis for the task experiments
%   psg_pwaffine_apply: apply a piecewise affine model
%   psg_pwprojective_apply:  apply a piecwise projective model
%   psg_pwprojective_test: test psg_pwprojective_appply, psg_geo_pwprojective_pset
%   psg_qformpred: predict perceptual space coords from quadratic form model of thresholds
%   psg_qformpred_demo: demonstrate predictions from quadratic form model of thresholds
%   psg_qform2coord_proc: create a coordinate file from a quadratic form model
%   psg_rayang_summ: summarize ray angles
%   psg_rayangles: compute angles between rays
%   psg_rayfit:  fit a coordinate structure to rays
%   psg_raymults: compute multipliers (gains) along rays
%   psg_raystats: parametric bootstrap with jitters for ray angles and multipliers
%   psg_raystats_dbplot: summarizes raystats tables, makes plots
%   psg_raystats_dbplot_axes: adjust axis scales interactively for psg_raystats_dbplot
%   psg_raystats_dbplot_style: utility script to set plot style in psg_raystats_dbplot
%   psg_raystats_dbplot_ticks: utility script for tick labels for plot panels in psg_raystats_dbplot
%   psg_raystats_summ: compute multiplers, angles, and their statistics across multiple datasets, make tables
%   psg_read_coorddata: read coordinates data inferred from a psg experiment
%   psg_remnan_coordsets: remove NaN's from aligned datasets
%   psg_selaxes_demo: create a data and metadata file with selected axes
%   psg_spec2legend: create a nice legend entry from spec_labels or typenames
%   psg_task_loaddata: load the data for all task experiments
%   psg_typenames2colors: assign colors to array types, for plotting
%   psg_typenames2colors_test: test psg_typenames2colors
%   psg_visualize: plot several pages of visualization of psg coords
%   psg_visualize_demo: demonstrate basic visualization of psg coords
%   psg_write_coorddata: write a coordinate data file and metadata
%
% Perceptual space geometry experiments: choice probability analysis
%   btcsel_like_analtable: as in psg_like_analtable, but for analyses with selected subsets of stimuli
%   ord_char_simul_cphists: plot choice probability histograms from ord_char_simuil_demo
%   ord_char_simul_demo: analysis of simulations of ultrametric, addtree, and Euclidean spaces
%   ord_char_simul_dist: generate distances for ord_char_simul_demo
%   ord_char_simul_plot: plot results fom ord_char_simul_demo
%   ord_char_simul_rmlabels: script for removing labels from plots
%   psg_choicedata_makeeven: prune choice probability data so that every triad has an even number of trials
%   psg_choicedata_merge: merge choice probability data files
%   psg_colors_like: set up default colors and symbols
%   psg_conform: determine how to flip a response to conform a dataset to sym, umi, or addtree
%   psg_dirichlet_loglike: expected log likelihood of trials with underlying Dirichlet prior for choice probabilties
%   psg_ineq_apply: apply the inequality conditions of psg_ineq_logic to a set of obsrvations, and do flips
%   psg_ineq_edgecount: counts the edges in output of psg_ineq_logic
%   psg_ineq_logic: sets up logic for excluded rank-choice-probabilities, for tests of symmetry, umi, addtree, etc.
%   psg_ineq_logic_compare: compares inequality structures (partitions) and a priori log likelihoods based on tents and tetrahedra
%   psg_ineq_logic_demo: test psg_ineq_logic
%   psg_indq_logic_scission: implement logic table for scission criteria
%   psg_ineq_logic_tetra: bootstrap an inequality structure from a tent to a tetrahedron
%   psg_ineq_lookup: look up a triadic comparison in a table
%   psg_ineq_tetra_test: test psg_ineq_logic_tetra and tetrahedral options in psg_ineq_logic
%   psg_ineq_triads: create a table of triadic comparisons
%   psg_like_analtable: analyze table from consolidated outputs of psg_umi_triplike_plota via psg_umi_trip_tent_run
%   psg_like_maketable: create table from consolidated outputs of psg_umi_triplike_plota via psg_umi_trip_tent_run
%   psg_logic_v_demo: test psg_ineq_logic and correspondence to manuscript
%   psg_permutes_logic: sets up permutations for tests of symmetry, umi, addtree
%   psg_probs_check: compares versions of psg_umi_triplike
%   psg_quad_stats: calculate and display statistics of quadruplets relevant to testing additivity and addtree
%   psg_quad_stats_demo: demonstrate psg_quad_stats
%   psg_read_choicedata: read choice probability data from a psg experiment
%   psg_resample_conform: draw samples from a Dirichlet distribution that are consistent with a set of inequalities
%   psg_resample_conform_demo: production execution of psg_resample_conform
%   psg_resample_conform_test: test psg_resample_conform
%   psg_select_choicedata: select the choice data from a subset of tokens
%   psg_stats_tally: utility for to tally statistics, for psg_triad_stats, psg_umi_stats
%   psg_tent_choices: extract tents of choice probabilities
%   psg_tent_stats: calculate and display tent statistics (condition for addtree)
%   psg_tent_stats_demo: demonstrate psg_tent_stats
%   psg_tentlike_demo: do max likelihood calculation for addtree for rank choice data
%   psg_triad_stats: calculate and display triad statistics
%   psg_triplet_choices: extract triplets of choice probabilities from choices data file
%   psg_umi_stats: calculate and display statistics of trials relevant to testing ultrametric inequality
%   psg_umi_stats_demo: demontrate psg_umi_stats
%   psg_umi_triplike: analyze likelihoods that rank choice probabilities are consistent with ultrametric inequality and symmetry
%   psg_umi_triplike_demo: apply psg_umi_triplike to data
%   psg_umi_triplike_plot: plot detailed results of psg_umi_triplike_demo
%   psg_umi_triplike_plota: plot summary (asymptotic0 results of psg_umi_triplike_demo
%   psg_umi_trip_tent_run: script to automate running of psg_umi_triplike_demo and psg_tentlike_demo
%
%  Example figures
%   bright_customplot_demo: custom plots for brightness (Aguilar) paradigm, unselected and selected
%   bright_customplot_demo2: further customization of brightness plots for gr23
%   btcsel_customplot_demo: demonstrate custom plots by merging and selecting database tables
%   btcsel_customplot_demo2: custom plots for bc6 unselected and selected, by selecting database tables
%   btcsel_customplot_demo3: further customization for ICERM ppt
%   faces_customplot_demo: custom plots of faces analysis, after selection by gender or age
%   faces_customplot_demo2: further customization of plots of faces analysis for gr23
%   psg_vss24_task_fig_meth: show a psg space in standard fashion with model, as points, and as distances
%   psg_vss24_task_fig_pool: show pooled bc55 data
%   psg_vss23_methfig: show threshold model, psg space, data points, and heatmap of distances, bgca d123 (same as psg_vss24_task_fig_meth)
%   psg_vss25_methfig: show threshold model, psg space, data points, and heatmap of distances, dgea d123, d234
%
%  Utilities/compatibility
%
%   psg_save_mat:  save a mat-file, for Octave compatibilty
%   psg_strfind: for Octave compatibility with strfind
%
%   Copyright (c) 2022, 2023, 2024, 2025 by J. Victor

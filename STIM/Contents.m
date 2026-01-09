% Visual stimuli m-files by J. Victor
%
% Utility.
%   bmpshow     - show an array of uint8 (typically read as a .bmp)
%   border      - put a border on a matrix, enlarging it, or crops it
%   borderin    - put a border on a matrix, not enlarging it
%   covtest     - demonstrate speedup of cov compared with non-vectorized approach
%   fixation    - put a circular fixation point on a bitmap
%   fopensaf    - open a file safely (checking for over-writing)
%   getfilelist - get a list of files conforming to a particular pattern from any path
%   longestr    - find the longest run of matches in a one-dimensional array
%   pxlrep      - pixel-replicate (enlarge) a matrix. Inelegant but faster than reppxl.
%   pxlsampl    - pixel-sample (shrink) a matrix. No antialiasing, averaging, etc.
%   qform_plot_util - plot a quadratic form
%   reprow      - replicate a multidimensional array along a row
%   repblk      - replicate each element of a multidimensional array into a block 
%   reppxl      - pixel-replicate (enlarge) a matrix, based on repblk.
%   sizepos     - interactively get size and position of a box
%   timer_test  - test Matlab's timer, showing how to retrieve the calling object
%   timer_test_func - callback for Matlab's timer, showing how to retrieve the calling object
%   zpad        - pad non-negative integer with zeros (for making file names)
%
% Glider manipulation, histo, stats, and maxent demos
%   cstats_energy     - high-order correlations (energies) for cstats_metro, multiple gray levels and any glider
%   cstats_metro_demo - demonstrates cstats_metro
%   cstats_metro_fig - creates a figure from cstats_metro
%   cstats_metro     - metropolis algorithm for making ME stimuli with constrained histograms and correlations
%   donut_metro      - metropolis algorithm with permutations within matching donuts
%   donut_metro_demo - demonstrates donut_metro:  permutations within matching donuts
%   donutw_metro      - variant of donut_metro that allows for weighted choices of permutations
%   donutw_metro_demo - demonstrates donutw_metro:  permutations within matching donuts
%   entropy_mtc_demo  - entropy for simple deviations from a random multiple-gray-level texture
%   fcdecomp_demo    - comparison of decomposition into Fourier components
%   flfs_bc_findp  - determine mixture of horizontal and vertical sticks to achieve desired weights for two-point correlations
%   flfs_bc_findp_of  - objective function for flfs_bc_findp
%   flfs_bc_findp_demo - test flfs_bc_bcfactors
%   flfs_define      - define structure for falling leaves/falling sticks constructions
%   flfs_define_perm - define "trivial" structure for falling leaves/falling sticks constructions, in which every permutation is a valid configuration
%   flfs_eivec_test  - for debugging the eigenvector-finding procedure
%   flfs_enumerate   - enumerate distinct sequences of falling leaves/falling sticks
%   flfs_enumerate_demo   - demonstrates flfs_[define|seqreduce|enumerate|makepropag],
%   flfs_findeivec   - find the stable eigenvector of a propagation matrix
%   flfs_makepropag  - make propagation matrices for one category of falling leaves/falling sticks, by tracking the configuration of sequences   
%   flfs_makepropag_cells  - make propagation matrices for one category of falling leaves/falling sticks, by the configuration of cells   
%   flfs_seqreduce   - reduce a sequence of falling leaves/falling sticks to the ones that show
%   glider_addcoords - create a glider structure from a matrix
%   glider_adjustp01 - adjust a probability array to [0 1] up to a tolerance
%   glider_getbcp    - get the parameters for the configutations of blocks
%   glider_cf2ca     - change block counts to independent coordinates on absolute blocks
%   glider_ca2cf     - change independent coordinates on absolute blocks to block counts
%   glider_ca2cr     - change independent coordinates on absolute blocks to independent coordinates on relative blocks
%   glider_cr2ca     - change independent coordinates on relative blocks to independent coordinates on absolute blocks
%   glider_choose    - create a glider library and choose which glider(s) to use
%   glider_test      - test this package
%   glider_addsubs   - add a description of the subgliders to a glider
%   glider_findsubs  - find the subgliders of a single rank
%   glider_fitmp     - find the best-fitting relative block probabilities consistent with a 2-d Markov process
%   glider_fitmp_afp - apply "forward" Pickard constraints
%   glider_fitmp_con - constraint function for glider_fitmp
%   glider_fitmp_of  - objective function for glider_fitmp
%   glider_fitmp_pl  - generate pointer lists for Pickard constraints for glider_fitmp
%   glider_lkpind    - find a check's index in a glider
%   glider_mapubi    - count the different colorings of a glider in an image
%   glider_mapubi_red    - reduce block counts to a smaller glider
%   glider_mee_demo  - demo of maximum-entropy extension analysis of images
%   glider_metro     - run Metropolis algorithm for maximum-entropy extension analysis of images
%   glider_metro_cadj  - check and adjust block probabilities for Metropolis algorithm
%   glider_metro_def - set defaults for Metropolis algorithm for maximum-entropy extension analysis of images
%   glider_metro_demo  - demo of Metropolis algorithm for maximum-entropy extension analysis of images
%   glider_metro_show  - show output of Metropolis algorithm for maximum-entropy extension analysis of images
%   glider_mp1d      - create a 1-d Markov process with specified block probabilities, arb. gray level and init
%   glider_mp2d      - create a 2-d Markov process with specified block probabilities
%   glider_pickapply - apply Pickard probabilities for a 2x2 box
%   glider_pickard   - calculate Pickard probabilities for a 2x2 box
%   glider_pcvtmtx   - calculate matrices to convert between relative, absolute, and full probabilities
%   glider_tm1d      - determine transition matrix and stable point of a 1-d Markov process
%   glider_2x2extend - extend from subblocks to Pickard-consistent probabilities
%   metro_demo       - demonstrate Metropolis algorithm for non-Pickard and Pickard textures
%   metro_demo_beta  - demonstrate Metropolis algorithm for textures with nonzero beta
%   metro_demo_theta - demonstrate Metropolis algorithm for textures with nonzero theta
%   pglider2markov   - determine the Markov matrix from block probabilities on any glider
%   ubi_codel_map    - create a map of unique block indices for individual local image statistics
%   ubi_codel_map_test - test ubi_codel_map
%
% Make stimuli and threshold models.
%   ag_domain_demo - show alpha, gamma domain and find bounaries along rays
%   alloy2strength - compute strength of quadratic mechanisms from alloy data, for mtc_fitted_ellipse (Rizvi)
%   arden       - make an arden plate
%   atg_block_scramb_demo - make examples of blocked and scrambled alpha-theta-gamma textures
%   at_mrf      - make texture with fourth- and third-order correlations, specified continuously
%   atg_space   - show the space of possible alpha-theta-gamma textures
%   atg_space7  - show the space of possible alpha-theta-gamma textures, for Matlab 7
%   bin2sr      - make 2x2 general block probabilities via stochastic recoloring
%   bin2tern_111 - make 2x2 ternary block probabilities from binary probabilities, w:g:b 1:1:1
%   bin2tern_121 - make 2x2 ternary block probabilities from binary probabilities, w:g:b 1:2:1
%   blevdemo    - demonstrates binary, even/odd propagated, and 2fold symmetric textures, b nonneg
%   binoise     - make binary noise
%   btca_deflayouts - set up defaults for (measured, predicted) experiments
%   btca_makeconds - make a conds array for (measured, predicted) experiments
%   btc_alpharange - determine possible values of alpha consistent with lower-order stats
%   btc_alpharange_bc_demo - demonstrate btc_alpharange and make maps with specified beta-h and beta-v
%   btc_alpharange_bg_demo - demonstrate btc_alpharange and make maps with specified gamma and beta-h
%   btc_asymstats_demo - tabulate asymmetries in psychophysical data
%   btc_atg_alt_demo - consider alternate ways of generating textures specified by a,t,g
%   btc_augcoords  - augment a subset of coordinates to the full 10-coord system, indicating a "method"
%   btc_augcoords_mix  - augment a subset of coordinates for a mixture of textures
%   btc_augcoords_mults  - augment a subset of coordinates to the full 10-coord system, with one or more multipliers, extrapolate if necessary
%   btc_autocorr_orient_demo - show oriented power in btc maps
%   btc_autocorr_pspec_demo - show autocorrelation and power spectrum of btc maps
%   btc_auxopts   - set up auxiliary options for map generation
%   btc_axesneeded - find the axes needed to characterize a metric
%   btc_axis_map_demo -- show btc textures along an axis, save in a mat-file (for Balasubramanian/Prentice)
%   btc_axis_plane_demo -- show btc textures along an axis and in a plane
%   btc_confluence_stats -- determine number of confluent regions in a texture
%   btc_coorkinds - determine what kinds of coordinates are present from texture coord letter codes gbcdetuvwa
%   btc_corrs2vec - convert a correlation structure (with gamma, beta, theta, alpha) to a 10-vector
%   btc_define   - define the binary texture blocks (2x2 system)
%   btc_dimrf    - create Markov random fields textures on diagonal lattice (for beta_diags)
%   btc_disting_est - estimate number of distinguishable textures
%   btc_d4matrix - generate a permutation matrix that represents the action of D4
%   btc_edirs    - set up a structure of experimental directions in a plane for binary experiments
%   btc_edirs_test - tests btc_edirs
%   btc_eivecs_def - sets defaults for btc_eivecs_reaqd, btc_[sym|hvi]dirlib
%   btc_eivecs_mink_compare:  compare eigenvectors used for stimuli with Minkowski directions
%   btc_eivecs_read - read file of eigenvector statistics, produced by btc_eivecs_stats
%   btc_eivecs_stats - summarize statistics of the eigenvectors found by btc_soid_demo
%   btc_eivecs_stats_madj - as btc_eivecs_stats, but defaults to madj=1 and also computes dot-products with self excluded.
%   btc_eivecs_stats_comp - compare statistics in two files, e.g., one with and one without madj=1
%   btc_eivecs_sum  - calculate statistics of eigenvectors
%   btc_eivecs_sumx - calculate statistics of eigenvectors comparing pairs of subjects
%   btc_exptname - determine standard experiment name string from texture coord letter codes gbcdetuvwa
%   btc_gamut    - show 1-d gamut for each btc parameter
%   btc_gamut_lim - create borders of 2-d gamut for a pair of params
%   btc_gamut_twostats    - show 1-d gamut for pairs of btc parameters
%   btc_genstim  - make maps with coordinates specified in strips (replaces of GenStimFastExec for btc series)
%   btc_genstim_de - make maps with coordinates specified in strips for the DE pair via a no-seam algorithm
%   btc_genstim_de_demo - demonstrate btc_genstim_de, including noseam option
%   btc_genstim_noseam_demo  - demonstrate "trynoseam" option, make maps without seam between adjacent nonrandom blocks, all planes EXCEPT de
%   btc_genstim_noseam_quaddemo  - demonstrate "trynoseam" option, four quadrants with different parameters (possibly for MDB pilot)
%   btc_getstats_demo - demonstraton of calculation of image statistics 
%   btc_hflip    - horizontal flip for a string of code letters
%   btc_hvi      - horizontal/vertical interchange (on NW to SE axis) for a string of code letters
%   btc_hvidirlib - generate library of directions in horizontal-vertical inversion space of btc
%   btc_hvirev   - horizontal/vertical interchange (on NE to SW axis) for a string of code letters
%   btc_letcode2vec - convert a letter code structure (with gamma, beta, theta, alpha) to a 10-vector
%   btc_letcode_strip - strip a letter code of non-numeric fields
%   btc_makebigconds_show - show the maps for a bigconds setup
%   btc_makemaps  - make maps from a method, such as pickard
%   btc_makemaps_scale  - make maps from a method, such as pickard, local coords at arbitrary scale
%   btc_makemaps_scale_test  - test btc_makemaps_scale
%   btc_makeonemap - make one btc map (only), used for debugging btcm series
%   btc_makepickard - complete a set of texture parameters to one that satisfies Pickard conditions
%   btc_makepickard_demo - demo/tryouts for completing Pickard conditions from gamma, beta_horiz
%   btc_makesem_demo - make a texture in the semantic d portion of the space (g=0,d=e,t=u,v=w)
%   btc_makespec_demo - make textures in special directoins of the space, including sym and hvi, and comparing
%   btc_makehvi_demo - make a texture in the horizontal-vertical invert 2-d portion of the space
%   btc_makesym_demo - make a texture in the 5-d (symmetric) portion of the space
%   btc_map2counts - calculate counts of block probabilities in a map
%   btc_map2counts_demo - simple demo for btc_map2counts (iid textures)
%   btc_mat_symmetrize - symmetrize a matrix with respect to one or more of the dihedral symmetries
%   btc_meabb    - find maxent alpha for betah and betav
%   btc_meabb_asymp - examine asymptotic behavior of maxent alpha for betah and betav
%   btc_meabt    - find maxent alpha for (betah or betav) and theta
%   btc_metro_demo - show stats of Metropolis algorithm working on btc textures
%   btc_metro_demo_replot - replot from btc_metro_demo
%   btc_metro_mix_demo - show stats of Metropolis algorithm working on mixtures of btc textures
%   btc_metro_mix_demo_replot - replot from btc_metro_mix_demo
%   btc_mix2tex_demo - create textures by Metrolis mixing of textures made with btc_vec2tex_pickard
%   btc_mix2tex_tuvw_showsamp_demo - show a sample of a tuvw texture created with btc_mix2tex_demo
%   btc_nometro_demo - stats of maps made without Metropolis algorithm
%   btc_nonstat_demo - demonstrate nonstationary nature of naive construction of tu texture via 2-dMarkov
%   btc_nopick   - nonPickard map-making routine for thetas sharing an edge or beta and theta
%   btc_offset_plan - plan an experiment to measure thresholds offset from the origin
%   btc_pairdomain - show the domain of a pair of coordinates, as a continuously varying texture
%   btc_pairdomain_token - show the domain of a pair of coordinates, as a continuously varying texture, standard and token
%   btc_pairextend - determine equivalent pairs of coordinates by rotation and hv interchange
%   btc_pairrays - show maps based on a pair of coordinates on rays from the origin
%   btc_pairsneeded - find the planes needed to characterize a metric
%   btc_pairtable - determine which axis pairs are covered by a set of paradigmns, including symmetries
%   btc_pairverify - verify construction of maps based on pairs of coordinates
%   btc_plan    - which planes are needed to extend the (ag) model
%   btc_precomp_demo - look at adjusting parameters of non-Pickard textures, so that Markov probabilities "come out right"
%   btc_pred_demo - predict thresholds on an axis and for mixture from a set of experimental directions and a quadratic form 
%   btc_pred_emink - predict thresholds in eigenvector and Minkowski directions from a set of experimental directions and a quadratic form 
%   btc_pred_emink_data - set up data structure for btc_pred_emink
%   btc_quad_bcde_demo - show statistics for textures with four second-order correlations
%   btc_quad_bcde_make - make examples of textures with four second-order correlations
%   btc_quad_bgca_demo - show statistics for textures with mixtures of first, second, and fourth-order correlations
%   btc_quad_bgca_make - make examples of textures with mixtures of first, second, and fourth-order correlations
%   btc_quad_bgca_onec_make - make examples of textures with mixtures of first, second, and fourth-order correlations, one-component strategy
%   btc_quad_bgca_onec_test - test how perturbing gamma affects other texture statistics
%   btc_quad_tuvw_make - make examples of textures with four third-order correlations
%   btc_qform_customize - customize a quadratic form model to fit a subject's thresholds
%   btc_quad_merge - merge resource files made with btc_quad*make
%   btc_qform_customize_run - process a file of quadratic form models and surrogates via btc_qform_customize
%   btc_qform_customize_test - test btc_qform_customize
%   btc_rend:  render an arbitrary [0 1] map in btc coordinates
%   btc_rend_noedge:  render an arbitrary [0 1] map in btc coordinates without edge effects, one btc coordinate specified
%   btc_rend_noedge_demo:  demonstrate btc_rend_noedge
%   btc_rend2_noedge:  render an arbitrary [0 1] map in btc coordinates without edge effects, two btc coordinates specified (Pickard pair only)
%   btc_rend2_noedge_demo:  demonstrate btc_rend2_noedge, only Pickard planes
%   btc_rend2_noedge_demo2:  demonstrate btc_rend2_noedge, de (diagMRF) also
%   btc_rotcode - rotate a string of code letters
%   btc_semstrats - how to split a "semantic" directio into two Pickard components
%   btc_soid_arch2ds - convert an xls file into a structure readable by btc_soid_getdata
%   btc_soid_arch2ds_demo - demonstrate btc_soid_arch2ds, btc_soid_arch2ds_demo
%   btc_soid_archscan - scan an xls file for correct formatting
%   btc_soid_demo - apply routines for finding best-fitting ellipsoid to psychophysical data
%   btc_soid_stats - version of btc_soid_demo with additional stats re fits on axes vs. diagonals
%   btc_soid_summ_avg - summarize results of btc_soid_demo and produce a similarly-formatted average subject file
%   btc_soid_ebfix - fix the error bars of an edirs structure
%   btc_soid_find - find the multiplier for a tangent vector to achieve a desired value of a quadratic form
%   btc_soid_find_of - objective function for btc_soid_find
%   btc_soid_find_mix - find the multiplier for a mixture of tangent vectors to achieve a desired value of a quadratic form
%   btc_soid_find_mix_of - objective function for btc_soid_find_mix
%   btc_soid_fit  - find the best-fitting ellipsoid to a psychophysical threshold dataset
%   btc_soid_getdata - get data psychophysical data for a set of planes from a matlab file
%   btc_soid_limits  - set up sanity limits for psychophysical data
%   btc_soid_plot - plot experimental setup and data for finding best-fitting ellipsoid
%   btc_soid_plotqf - plot quadratic form corresponding to best-fitting ellipsoid
%   btc_soid_predict - predict thresholds from a set of experimental directions and a quadratic form
%   btc_soid_summ_read - read a file of fits of psychophysical data to quadratic model (of the sort written by btc_soid_demo)
%   btc_soid_sym - symmetrize a btc psychophsyical dataset, and uniformize axes that occur in multiple planes
%   btc_soid_symtest - test btc_soid_sym
%   btc_soid_test - test routines for finding best-fitting ellipsoid
%   btc_soid_xlsfmt  - set up formats for reading a spreadsheet of psychphysical threshold data
%   btc_soid_xlsread - read spreadsheet with threshold data
%   btc_soid_xlsread_demo - demonstrate btc_soid_xlsread
%   btc_soid_xlstst - test the format and contents of a psychophsycial data spreadsheet
%   btc_soidfg_anal - analyze models fitted by btc_soidfg_fit
%   btc_soidfg_anal_demo - test btc_soidfg_anal
%   btc_soidfg_anal_of - objective function for projection in btc_soidfg_anal
%   btc_soidfg_calcstats - calculate statistics of a figure-ground fit
%   btc_soidfg_define - define options structure for fits to psychophysical figure-ground data
%   btc_soidfg_demo - demonstrate quadratic fit to figure-ground data modules
%   btc_soidfg_find - find the threshold along a direction
%   btc_soidfg_find_of - objective function for btc_soidfg_find
%   btc_soidfg_fit - carry out quadratic fit to figure-ground data
%   btc_soidfg_fit_makex - set up a regression fit for quadratic figure-ground model
%   btc_soidfg_model - set up model for fits to psychophysical figure-ground data
%   btc_soidfg_model_proj - project a phenomenological model into the ideal-observer space
%   btc_soidfg_model_test - test btc_soidfg_model, btc_symclasses
%   btc_soidfg_modelranks - calculate model ranks for phenomenological soidfg models, before and after projection into ideal-observer space
%   btc_soidfg_name2lets - extract coordinate letters from a model parameter name
%   btc_soidfg_predict - predict a figure-ground threshold from a quadratic model
%   btc_soidfg_symmetrize - symmetrize psychophysical data based on sign
%   btc_soidfg_validate - validate thresholds and eror bars in a psychophysical dataset
%   btc_split_alpha - adjust alpha parameters to try to keep probabilities non-negative
%   btc_sum_eccen - summarize rmse of fits and characterize shapes of ellipses in a fitted model
%   btc_surv      - survey the available methods for coordinate augmentation
%   btc_symbasis - find a symmetrized basis for the btc coordinates
%   btc_symclasses - find classes of btc coordinates equiavalent under various symmetries
%   btc_symdirlib - generate library of directions in symmetric space of btc
%   btc_symstrat - implement one of several strategies for making symmetric textures
%   btc_symstrat_demo - demonstrates and compares strategies for making symmetric textures
%   btc_symstrats - find best strategy for making symmetric textures
%   btc_sym_proj_tensor - find projection and tensor for symmetrizing with respect to one or more of the dihedral symmetries
%   btc_test     - test btc-series
%   btc_token_edirs - get stimulus coordinates for token ellipse experiments
%   btc_token_xlstst - test reading of summary spreadsheet from token ellipse experiments
%   btc_token_xls2ds - create mat structures raedable via btc_soid_getdata from token ellipse experiments
%   btc_trule_explore - determine whether T-rules are self-consistent
%   btc_trule_explore2 - another approach to determine whether t-rules are self-consistent, 2 generations
%   btc_trule_explore3 - another approach to determine whether t-rules are self-consistent, 3 generations
%   btc_trule_getcorrs  - get correlations from a T-rule 2x2x2x2 set of probabilities (ABC top row, then D)
%   btc_trule_iterate - iterate a T rule
%   btc_trule_search1 - Nelder-Mead search for T rules, constraining probabilities within 2x2 block
%   btc_trule_search1of - objective function for btc_trule_search1
%   btc_vecnan2letcode - convert a 10-vector to a letter code structure (with g, b, c, d, e, t, u, v, w, a), omitting NaN's in the vector
%   btc_vec2corrs - convert a 10-vector to a correlation structure (with gamma, beta, theta, alpha)
%   btc_vec2letcode - convert a 10-vector to a letter code structure (with g, b, c, d, e, t, u, v, w, a)
%   btc_vec2tex_pickard - create a texture from a 10-vector that satisfies the Pickard conditions
%   btc_vflip    - vertical flip for a string of code letters
%   btc_whichsym - determine the symmetries of a set of block probabilities
%   btc_3x3anal - propagate a 2x2 rule into 3x3 via Markov, to allow checking for consistency
%   btcm_super  - supervisor and doc for "btcm" routines to run metropolis algorithm on btc stimuli (Jason Mintz)
%   btcstats - simple calculation of binary texture coordinates from a binary image
%   btcstats_demo - demonstrate btcstats
%   changebstats - change the statistics of a binary image
%   changebstats_gamma - change the gamma statistic of a binary image
%   clevdemo    - demonstrates binary and even/odd propagated textures, c neg or pos
%   clevdemo_atg  - demonstrates general atg textures, c neg or pos, via genmrfmg and getp2x2_atg
%   clevdemo_atbg  - demonstrates general atbg textures, c neg or pos, via genmrfm and getp2x2_atg, getp2x2_pickard1
%   clewdemo    - demonstrates binary, even/odd propagated, and 2fold symmetric textures, c neg or pos
%   convert_bigconds_whichdir - modify a BigConds file so that each correct answer lies along a ray, and that target flag assumes both values
%   convert_bigconds_whichdir_pool - modify a BigConds file so that selected rays are pooled, and that target flag assumes both values
%   decode_map  - decode a map created with ptbxcvt_btcdemo into its coordinates
%   decode_map_demo  - demonstrates decode_map
%   decode_map_maketables  - use decode_map to create configuration tables for graylevel experiments
%   dg_btc_pic  - explore properties of dichotomized Gaussians
%   d1strips    - make strips of 1-d correlated textures with controlled flips
%   d1flipqu    - make quads of 1-d correlated textures with controlled number of flips
%   d1maps      - make maps of 1-d correlated textures; no consideration of flips
%   diagrule_ptov - convert a "diagonal" rule from block probabilities to rule number as vector
%   diagrule_show - show a sample of a "diagonal" rule from rule number, arbitrary number of gray levels
%   diagrule_stov - convert a "diagonal" rule from rule number as scalar to rule as vector
%   diagrule_vtos - convert a "diagonal" rule from rule as vector to scalar rule number
%   dilumetro_test - demonstrate and test Metropolis donut routine for diluting texture params
%   dg_btc_pic  - relationship of dichotomized Gaussian and btc parameterizations, and Pickard conditions
%   eig_mink_define  - define block structure and contrast levels for eigenvector/Minkowski experiment
%   eig_mink_setup  - automated setup (adaptation of ptbxcvt_btcdemo) of blocks for eigenvector/Minkowski experiment
%   eo_glider   - make binary textures from arbitrary gliders, with propagated and sporadic errors
%   eobademo    - demonstrates eobanal and evenoddb
%   eobagind    - makes even/odd/random textures (prop/spor or Markov random fields) specified via indices
%   eobagmat    - makes even/odd/random textures (prop/spor) with eo and lum bias specified continuously
%   eobagmrf    - makes even/odd/random textures via Markov random fields with eo and lum bias specified continuously
%   eobanal     - analyzes parameters for evenoddb
%   eobbdemo    - demonstrates eobanal and evenoddb, with 4 separate fields with edges in common
%   eobdemo     - demonstrates evenoddb
%   eobmrf      - generates uniform even/odd/random textures many ways, demonstrating pitfalls if not isotripole
%   eobmagst    - makes a map with demonstrates eobmrfa, eobmrfm
%   eobmdemo    - demonstrates eobmrfa, eobmrfm
%   eobmagst    - makes even/odd/random textures (Markov random fields) with parameters specified in strips
%   eobmrfa     - analyzes parameters for Markov random field generation of eo and lum bias
%   eobmrfbp    - calculate block probabilities for a Markov random field texture
%   eobmrfm     - makes even/odd/random textures (Markov random fields)
%   eoflipga    - make a gallery of pairs of propagated textures with controlled number of interchanges
%   eoflipqu    - make a gallery of quads of propagated textures with controlled number of interchanges
%   eoreduce    - reduce a texture to a sparse texture by flipping against an even texture
%   eoworst     - look for "worst cases" for reduction by eoreduce
%   eqindemo    - samples of stimuli from information-theoretic-equalized classes
%   evenodd     - make even, odd, and random textures, with sporadic error
%   evenoddb    - make even/odd/random textures (prop/spor) with eo and lum bias, can specify first row, first col
%   evenoddp    - make even, odd, and random textures, with sporadic and propagated error
%   evenodd_demo- standalone demonstration of even/odd/random (prop/spor), slow and fast generation
%   genmrfm     - generate binary Markov random fields from general 2x2 block probabilities (not checking for consistency)
%   genmrfm_tee - generate binary Markov random fields from block probabilities in a T-shaped region (not checking for consistency)
%   genmrfm_tee_trap - generate binary Markov random fields from block probabilities in a T-shaped region, using first row only to initialize
%   genmrfm_1d  - generate 1-d binary Markov random fields; these are 1-d initializations for genmrfm
%   genmrfmg    - generate non-binary Markov random fields from general 2x2 block probabilities (not checking for consistency)
%   genmrfmg_tee  - generate non-binary Markov random fields from block probabilities in a T-shaped region, using first row to initialize
%   genmrfmg_1d - generate 1-d non-binary Markov random fields; these are 1-d initializations for genmrfmg
%   genmcig     - generate binary Markov random fields for a CIG texture
%   gen_eyeris_beta_prelim - generate maps for the Rucci EyeRIS system, cardinal betas only, preliminary version of gen_eyeris_maps
%   gen_eyeris_maps - generate maps for the Rucci EyeRIS system, cardinal betas only
%   getcorrs_allnames- create a structure specifying all 2x2 param names
%   getcorrs_diagrule- get correlations from specification of a "diagonal" rule (per Chubb)
%   getcorrs_farness - compare two sets of correlations at multiple scales
%   getcorrs_p2x2   - get correlations from 2x2 block probabilities of a binary rule
%   getcorrs_rescaleb   - get correlations after rescaling a binary texture
%   getcorrs_rescaleb_demo   - demonstrates getcorrs_rescaleb on textures in the library made by texturelib_get
%   getcorrs_scales  - get correlations at multiple scales, via getcorrs_rescaleb
%   getftps_p2x2    - get Fourier transform coordinates from 2x2 block probabilities
%   getmoms_p2x2   - get moments from 2x2 block probabilities, not necessariliy binary rule
%   getp2x2_abg     - get 2x2 block probabilities from first- second, fourth -order correlations
%   getp2x2_ag     - get 2x2 block probabilities from first- and fourth -order correlations
%   getp2x2_atg     - get 2x2 block probabilities from first- third, fourth -order correlations
%   getp2x2_cig - get 2x2 block probabilities from Champagnat Idier Goussard parameters (non-Pickard)
%   getp2x2_diagrule - get 2x2 block probabilities for "extreme" (diagonal) rules
%   getp2x2_corrs   - get 2x2 block probabilities from correlations
%   getp2x2_ftps    - get 2x2 block probabilities from Fourier transform coordinates
%   getp2x2_pabcde  - get 2x2 block probabilities from a structure of params for eobmrf series
%   getp2x2_pickard1- get 2x2 block probabilities from seven correl parameters in first Pickard component
%   getp2x2_tg      - get 2x2 block probabilities from first- and third-order correlations
%   getpabcde_ag    - get structure of params for eobmrf series from first- and fourth-order correlations
%   grating     - make a grating
%   gtc_define  - define the binary texture coordinates for any neighborhood shape
%   howfevod    - determine how many flips are possible for even-to-even, even-to-odd conversions
%   iso34sup_demo   - make superpositions of third- and fourth-order textures at different scales
%   mdbcirc_demo - demonstrate an experiment using mdbcirc_gen, mdbcirc_trial
%   mdbcirc_diag - diagram of an experiment using mdbcirc_gen, mdbcirc_trial
%   mdbcirc_gen  - generate stimuli with a border between two or three sinusoidally varying texture parameters
%   mdbcirc_getdata  - get data from one or more sets of trials from mdbcirc_demo
%   mdbcirc_getkinds  - get kinds of trials to run in mdbcirc_demo
%   mdbcirc_graybar   - add a gray bar to an image
%   mdbcirc_mefill - fill an (alpha, theta, gamma) triplet via maximum entropy
%   mdbcirc_show - show a sequence of stmimuli made by mdbcirc_gen
%   mdbcirc_test - test mdbcirc_gen, mdbcirc_show, mdbcirc_trial
%   mdbcirc_timer_func - timer function for mdbcirc_trial
%   mdbcirc_trial - run a trial based on a sequence of stimuli made by mdbcirc_gen, get results by waitforbuttonpress
%   mdbcirc_trial_ginput  - run a trial and get results via GINPUT
%   mdbcirc_verify    - verify properties of the mdbcirc stimuli by analyzing many samples
%   mee_decomp_demo - demonstrate decomposition of complex binary 2x2 textures into several one-parameter textures
%   mpentag     - entropy per unit area of the (alpha,gamma) or (alpha, gamma, theta) texture
%   mpentkld    - plots entropy per unit area and K-L divergences of the (alpha,gamma) texture
%   mpentkld2   - pretty plots of entropy per unit area and K-L divergences of the (alpha,gamma) texture
%   mpentkld_atg  - plots entropy per unit area and K-L divergences of the (alpha,gamma,theta) texture
%   mpgrent     - calculate generalized relative entropy (for characterizing MRF textures)
%   mpkldag     - K-L divergence between two (alpha,gamma) textures
%   mpkldagt    - table of K-L divergence between two (alpha,gamma) textures, for IO calculation in luis.doc
%   mrfen       - find maximum entropy mrf even/odd/biased textures, with isotripole, isodipole, and laxer constraints
%   mrfen_df    - set up default structure for mrfen and mrfen_of
%   mrfen_of    - objective function for mrfen
%   mrfent22    - find block entropy of a Markov random field texture
%   mtc_addbias - add bias to a coordinate group structure in specified directions
%   mtc_addbias_cvt - convert a probability vector to a coordinate group structure suitable for adding bias (needed for ng composite)
%   mtc_addbias_cvt_test - tests mtc_addbias_cvt
%   mtc_augcoords - augment coordinates and generate strategies for multiple-gray-level-coordinates (work in progress)
%   mtc_augcoords_bc_boundary_demo -- shows the fuzziness of the border in one of the bc planes for the consensus method
%   mtc_axis_plane_verify - tests mtc syntheses in all axes and planes
%   mtc_bc_de_test - test stimulus generation in bc and de planes
%   mtc_bc_de_flfs_tryout - try out syntheses in beta, beta plane esp for non-Pickards, including diabonals, using falling-leaves construction
%   mtc_bc_flfs_makemaps_test - test mtc_flfs_makemaps for bc direction
%   mtc_bc_flfs_setup - utility for mtc_bc_flfs_tryout, mtc_bc_de_flfs_tryout
%   mtc_bc_flfs_tryout - try out syntheses in beta, beta plane esp for non-Pickards, using falling-leaves construction
%   mtc_bc_me_search - search for starting conditions in the bc plane, special cases for ng=2,3
%   mtc_bc_me_tryout - try out various approaches for the bc plane
%   mtc_bc_me_tryout2 - demonstrate an example of params that are Pickard in both directions but lead to different transition matrices
%   mtc_bc_me_tryout3 - try out approaches for the bc plane, focusing on maxent routine when both diags are Pickard
%   mtc_bd_me_tryout - try out various approaches for the bd plane
%   mtc_bd_thetaopts_demo - demonstrate the two "theta" options for augmenting coords in beta_card, beta-diag, 
%   mtc_bt_nopick_tryout - test methods for Pickard and non-Pickard  beta, theta plane
%   mtc_btc_model_demo -  combined modeling and plotting of gray-level and binary thresholds
%   mtc_cgs_check - check non-negativity, normalization, in a multiple-gray-level coordinates,and fix when possible
%   mtc_cgss2probs - convert multiple-gray-level coordinates to probabilities of each configuration, using sparse (intuitive) method
%   mtc_cgss_reduce - project multiple-gray-level coordinate groups to a sparse subset
%   mtc_cgst2probs - convert multiple-gray-level coordinates to probabilities of each configuration, using totient method
%   mtc_cgst_reduce - project multiple-gray-level coordinates via the totient method
%   mtc_cgs2cgsstruct - convert a matrix of values to a coordinate group structure (named values)
%   mtc_cgs_remove - remove coordinates of one or more orders from a multiple-gray-level coordinate array
%   mtc_cgs_rmclass - remove a coordinate class from a multiple-gray-level coordinate array
%   mtc_cgs2btc_xform - create a matrix to convert multiple-gray-level coordinates to btc coordinates, taking care of reordering and sign conventions
%   mtc_cgs2probs_xform - create a matrix to convert multiple-gray-level coordinates to probabilities of each configuration, using sparse (intuitive) method
%   mtc_cgs2trico - convert coordinate group values to trico probability values
%   mtc_cgsstruct_check - check non-negativity, normalization, and names in a cgs structure,and fix when possible
%   mtc_cgsstruct_clean - clean a cgs structure by removing unbiased coordinate groups
%   mtc_cgsstruct_fill -  fill in the unbiased coordinate groups in a cgs structure
%   mtc_cgsstruct_merge - merge one coordinate group structure into another
%   mtc_cgsstruct_remove - remove coordinate groups of one or more orders from a cgs structure
%   mtc_cgsstruct_rename- rename the fields of a coordinate group structure when the template is changed
%   mtc_cgsstruct_reverse - determine the coordinate structure resulting from reversing the check positions
%   mtc_cgsstruct_rmclass - remove coordinate class from a cgs structure
%   mtc_cgsstruct_project - project coordinate group structure to a unique set of values
%   mtc_cgsstruct_me_demo - test mtc_cgsstruct2gcs_me for various gray levels and gliders
%   mtc_cgsstruct_me_test - test mtc_cgsstruct2gcs_me for ng=2, comparing against analytic values
%   mtc_cgsstruct2cgs - convert a coordinate group structure (named values) to a matrix of values
%   mtc_cgsstruct2cgs_me - find the maxent coordinate group structure for given coords of order nchecks-1 or less
%   mtc_cgsstruct2cgs_me_of - objective function for mtc_cgsstruct2cgs_me
%   mtc_choosecoords - interactively choose coordinate specifications, typically within a plane
%   mtc_countparams - count and summarize the free parameters in the coordinate groups
%   mtc_coord_rename - rename coordinates or coordinate groups for a new template
%   mtc_define  - set up multiple-gray-level coordinates
%   mtc_define_test  - test mtc_define
%   mtc_defopts - set up options for mtc package (tolerances, logging)
%   mtc_design_show - show a BigConds array, both in terms of raw params and triangle (alloy) transformation
%   mtc_design2conds - make a BigConds array from a design description
%   mtc_edata_mergeplanes - merge psychophysical data from standard directions and interleaved
%   mtc_edirs    - set up a structure of experimental directions in a plane for multigray-level experiments
%   mtc_edirs_merge   - merge two planes of experimental data (e.g., standard and interlaved) multigray-level experiments
%   mtc_fitted_ellipse - find the best-fitting ellipse to alloy-plot data, and determine mechanism strengths (Rizvi)
%   mtc_flfs_makemaps - make multilevel maps via falling leaves/falling sticks construction
%   mtc_freeparams - find the free parameters in a coordinate group structure
%   mtc_gbc2let - convert a gamma or cardinal beta coordinate axis to a letter code (for mtc_plans.xlsx)
%   mtc_gbc2let_demo - demonstrate mtc_gbc2let, mtc_let2gbc: make a table of conversion of gamma or cardinal beta coordinate axes to a letter code
%   mtc_genstim: upgrade of btc_genstim, make maps with patches specified by cgs structures
%   mtc_genstim_de - upgrade of btc_genstim_de, make maps with coordinates specified in strips for the DE pair via a no-seam algorithm
%   mtc_genstim_noseam_demo - demonstrate "trynoseam" option, make maps without seam between adjacent nonrandom blocks
%   mtc_genstim_test - test mtc_genstim
%   mtc_getcorrs_p2x2 - check entropies and normalizations for a set of 2x2 probabilities
%   mtc_getdirvec - get a direction, either a pure direction or arbitrary bias
%   mtc_let2gbc - convert a letter code to a gamma or cardinal beta coordinate axis (for mtc_plans.xlsx)
%   mtc_make_coord_table - utility for mtc_cgs2probs, mtc_probs2cgs
%   mtc_makeconds - make a BigConds array for plane experiments with axes specified by mtc
%   mtc_makeconds_demo - demonstrate mtc_makeconds, and make sample maps
%   mtc_makeconds_test - test mtc_makeconds
%   mtc_maxaug - find the maximal deviation from the unbiased point for which a method exists
%   mtc_maxaug_sensitiv_demo - demonstrate sensitivity of maxaug to rounding error
%   mtc_method_choose - interactively choose from a list of methods for augmentation
%   mtc_mtcinfo - create an "mtcinfo" structure to allow compatibility between btc and mtc series
%   mtc_onecg_demo - demonstrates synthesis of multiple-gray-level textures in a single coordinate group
%   mtc_parsename - parse the name of a coordinate or coordinate group, into a list of checks and a list of multipliers
%   mtc_plane_demo - demonstrate coordinate planes for multiple-gray-level coordinates
%   mtc_probs2cgs - convert probabilities of each configuration to multiple-gray-level coordinates
%   mtc_probs2cgs_xform - create a matrix to convert probabilities of each configuration to multiple-gray-level coordinates
%   mtc_ramp_demo - demonstrate maps with pairwise ramp correlation structure
%   mtc_ramp_frac - utility to generate probability vector for ramp and streaks, fractional step sizes to allow for arbitrary repeat period
%   mtc_ramp_frac_demo - demonstrate maps with pairwise ramp correlation structure and streaks, fractional step sizes
%   mtc_ramp_frac_eucl - Euclidean distances for maximum-contrast ramps, with fractional step sizes
%   mtc_ramp_seam_demo - demonstrate maps with pairwise ramp correlation structure and seams vs noseams
%   mtc_ramp_streak_demo - demonstrate maps with pairwise ramp correlation structure and streaks
%   mtc_rays_demo - demontrate maps along one or more rays, intended for multiple gray levels
%   mtc_samps_demo - make sample maps and save as png
%   mtc_selectmeths - select methods from the output of mtc_augcoords
%   mtc_simpmodels_demo -- isodiscrimiation contours for simple discrimination models
%   mtc_simple_plot_demo:  simplified version of mtc_soid_plot_demo to demonstrate access to psychophysical thresholds
%   mtc_soid_arch2ds_demo - demonstrate btc_soid_arch2ds, btc_soid_arch2ds_demo for gray-level data
%   mtc_soid_avg - compute several types of averages across subjects
%   mtc_soid_plot_demo:  plot merged psychophysical data
%   mtc_soid_xlsfmt  - set up formats for reading a spreadsheet of gray-level psychphysical threshold data 
%   mtc_soid_xlsrun - read gray-level psychophyscial data spreadsheets and create an omnibus mat file of all fits
%   mtc_soid_xlsrun_summ - merge, summarize, and plot output of mtc_soid_xlsrun
%   mtc_soid_xlsrun_summ_34 - merge, summarize, and plot output of mtc_soid_xlsrun, use planes of {1,2,3}, {1,2,4} order to find common scaling across subjects
%   mtc_soid_xlstst - test the format and contents of a gray-level psychophsycial data spreadsheet
%   mtc_thresh_check -- check the fitted threshold values,to verify coordinate transforms, BigConds, etc.
%   mtc_trule_531 - iterate a t-shaped rule in a triangular region, to check stability
%   mtc_tt_nopick_tryout - test methods for Pickard and non-Pickard theta, theta planes
%   mtc_twocplane_demo - demonstrate maps in a plane of two coordinates, for designing experiments
%   mtc_vec2cgsstruct_test - convert a probability vector into a structure, dealing with the nontrivialities of a composite ng
%   mtc_unbiased - create the unbiased probability vectors
%   mtc_xalloy_mtx - set up a matrix for display of raw coords as alloy plot
%   mtc_xform_demo - demonstrate mtc_cgs2probs_xform and mtc_probs2cgs_xform
%   mxe22prb    - find the 2x2 blocks with maximum entropy, given prescribed nearest-nbr flips
%   mxe22sub    - required for mxe22prb
%   randflip    - apply the same random flip (and rotation) to all maps in a stack
%   rfp_btc_demo - demonstrate rfp_btc_make
%   rfp_btc_explore - demonstrate rfp_btc_make and rfp_curve, one map at a time, interactive
%   rfp_btc_make - render radial frequency patterns with btc textures
%   rfp_curve   - make a radial frequency pattern as a curve
%   rfp_defopts - define options for radial frequency pattern generation
%   rfp_findr   -  find the radius at given orientations along an rfp_curve
%   rfp_mask    - make a radial frequency pattern as a mask
%   rfp_synth_demo - demonstrate rfp_curve, rfp_mask, rfp_findr
%   sqtoken     - make a square token
%   st_glider   - make binary textures from arbitrary spatiotemporal gliders, with propagated and sporadic errors
%   symmdemo    - demonstrates various twofold symmetry classes
%   tbln_demo   - demo of truncated band-limited noises
%   tbln_harvest- show results (entropies, correlations) from tbln_demo
%   tda_btc_examples - make examples of btc stimuli used in topological data analysis paper (Rodrigues, Purpura, Guidolin, Desroches) 
%   terndemo    - demonstrates terntex
%   terntex     - make some ternary textures
%   texect      - count number of possible stimuli given specific number of errors
%   texedemo    - demonstrates texect
%   texres_entbounds - determine bounds on entropy of a texture in texres
%   texres_getlist - get a list of existing texture resource (texres_*.mat) files
%   texres_getmap - get a map from a texture resource
%   texres_make - make texture resource (texres_*.mat) files
%   texres_make_fix - fix texture resource by eliminating multipliers with no maps
%   texres_make_refine - refine texture resource (texres_*.mat) files, focusing on on maps that are poorly mixed
%   texres_quality - use JS divergence (entropy of mixing) to quantify quality of texture resources
%   texres_setup - set up the definitions of the texture resources to make
%   texres_show - show selections from texres
%   texturelib_anal - analyze contents of texture library
%   texturelib_desc_sipi - get descriptors of SIPI texture library
%   texturelib_desc_pvt - get descriptors of private texture library
%   texturelib_get - get a texture from the library (specifies histogram equalization, cropping, etc.)
%   texturelib_mee_demo - read some texture and make 2x2 MEE approximations via donut_metro
%   texturelib_mee_monitor - monitor progress of texturelib_mee_demo
%   texturelib_meew_demo - read some texture and make 2x2 MEE approximations via donutw_metro
%   texturelib_meew_monitor - monitor progress of texturelib_meew_demo
%   texturelib_read - read a natural texture from the library
%   tg_anal     - analyze parameter space for first- and third-order correlations, related to getp2x2_tg
%   tg_mrf      - make texture with first- and third-order correlations, specified continuously
%   tg_domain_demo - show theta, gamma domain and find bounaries along rays
%   trico_axes_demo - demonstrate axes of ternary 2x2 space
%   trico_cvt   - convert among ternary probability specifications
%   trico_fc2probs - convert Fourier coefficients to pvals sfor ternary probability specifications
%   trico_generate_demo - generates textures from ternary probability specifications
%   trico_probs2fc - convert pvals to Fourier coefficients for ternary probability specifications
%   trieo       - make a triangle/even hybrid texture (Chubb texture)
%   twomdemo    - demonstrates stimuli that flip between two symmetry classes
%   usetoken    - use a set of tokens to elaborate a texture
%   usetdemo    - demonstrate usetoken
%   vismoncal   - visual calibration of monitor
%   vismoncal_demo   - demonstrates calibration of monitor
%   vismoncal_meas   - make a measurement for visual calibration of monitor
%   vismoncal_rgb_demo   - demonstrates calibration of monitor, r, g, and b separate
%   vismoncal_keypress_func   - visual calibration of monitor keypress callback
%   vismoncal_timer_func      - visual calibration of monitor timer callback
%
% Manipulation of local image stats
%   mlis_alglib - set up library of algorithms for mlis routines
%   mlis_alglib_compare - compare modified patches made b mlis_alglib_pilot, several algs on one page
%   mlis_alglib_filename_aug - utility to add path information to stimulus files made by mlis_alglib_pilot
%   mlis_alglib_local -analysis and display of local (pixel-wise) modifications made by mlis_alglib_pilot
%   mlis_alglib_pilot - pilot algorithms on multiple images, multiple btc targets; generate database and stmiulus patches
%   mlis_alglib_pilot_init - initialization for mlis_alglib_pilot
%   mlis_alglib_plot - plot results of mlis_alglib_pilot
%   mlis_alglib_reg_demo - demonstrate region-by-region image modifications
%   mlis_alglib_reg_show - show results (maps and statsitics) of region-by-region image modifications
%   mlis_alglib_sdbcheck - check integrity of sdb (stimulus database) structure 
%   mlis_alglib_sdbwrite - utility to write sdb structure
%   mlis_alglib_summ - summarize sdb (stimulus database) from mlis_alglib_pilot
%   mlis_alglib_summ_plot - plot summaries of sdb (stimulus database) from mlis_alglib_pilot, run after mlis_alglib_summ
%   mlis_alglib_vecreorg - utility for mlis_alglib_pilot
%   mlis_btcstats - binarize and calculate statistics from a gray-level map (no whitening)
%   mlis_define - define options
%   mlis_local_show2jh - utility for mlis_alglib_local to display two joint histograms
%   mlis_mirror - augment a map by mirroring in both axes
%   mlis_regs_ask - get aoptions for sub-region manipulation
%   mlis_regs_locs_setup - set up locations of sub-regions for image manipulation
%   mlis_regs_locs_demo - demonstrate mlis_regs_locs_setup, mlis_regs_weights_setup
%   mlis_reg_weights_setup - set up weights for mixing sub-regions for image manipulation
%   mlis_run - run the phase-jitter (search) algorithm
%   mlis_run_alg - run general map modification, phase substitution followed by phase-jitter (search)
%   mlis_run_b  - run a proof-of-principle phase-substitution algorithm
%   mlis_run_calcdist - calculate distance of btc params to target
%   mlis_run_restoremv - restore a map's mean and variance
%   mlis_run_setup - set up whitening filter, cosine bell, phase mask for mlis_run and mlis_run_sub
%   mlis_run_sub  - run a phase substitution algorithm
%   mlis_sdbseq_acmd - set up a sequence of files from an sdb database the ACMD (Algorithm Comparison for Manipulation of Digital Mammogram) study
%   mlis_sdbseq_design - choose a subset of patches for stmiuli, satisfying constraints on non-reuse of sources
%   mlis_sdbseq_usource - set up a list of unique sources fro stimuli based on metadata generated by mlis_sdbseq_acmd
%   mlis_test   - test mlis_run (phase searching)
%   mlis_test_b - test mlis_run_b
%   mlis_test_sub - test mlis_run_sub (phase substitution)
%   mlis_test_suse - test mlis_run_sub followed by mlis_run
%   mlis_warn - check compatibility of an options structure
%   mlis_whitening_filter - create whitening filter for a map, using options in mlis_define
%   mlis_window_setup - create standard windows (flat, cosbell, Gaussian)
%   whiten_window_demo - simple demo of whitening with cosine bell
%
% Response modeling.
%   sm_clim     - estimate confidence limits on fitted parameters
%   sm_couer    - estimate counting error
%   sm_deff     - set up default parameter fitting options
%   sm_defo     - set up default stimulus options
%   sm_defp     - set up default model params
%   sm_demo     - create a simulated data set
%   sm_demor     - create a simulated data set, set up stimulus and modeling options, show resids
%   sm_fcor     - fraction correct, for an array of target corr vals vs one distractor corr val
%   sm_fcorp    - plotting routine for sm_fcor and sm_fitm
%   sm_fitm     - fit a model, calls sm_fcor
%   sm_fitmr    - find a residual for model fitting
%   sm_hess     - find Hessian at the minimum
%   sm_mgof     - test goodness of fit (following sm_summ)
%   sm_mtda     - histogram of maximum number of targets detected several arrays (calls sm_tdsa)
%   sm_nicep    - make nice plots of data and fit
%   sm_summ     - summarize and consolidate multiple analyses (category 2, 3, 5)
%   sm_summ2324 - summarize and consolidate multiple analyses (category 23, 24, 24neg)
%   sm_tdsa     - histogram of number of targets detected within summation area (calls sm_tesa)
%   sm_tep      - histogram of number of targets that are present
%   sm_tesa     - histogram of number of targets present within summation area (calls sm_tep)
%       Note additional documentation in sm_doc.txt
%   sm_where    - finds where an array has pixels that satisfy a constraint
%   tdpmodel    - simple model with local processing and global limit
%   tdpmtarg    - version of tdpmodel for finding critical values as function of btarg
%   tdpmdist    - version of tdpmodel for finding critical values as function of bdist
%   tdpmdemo    - demonstrates tdpmodel
%   weibull_gfit - log-likelihood fit of one or more datasets with individual Weibulls or a common b
%
% Make a set of stimuli for talks.
%   makebino    - make binary noises, series of contrasts and check sizes
%   makeevod    - make even, odd, random textures
%   makegrat    - make gratings, series of orientations and spatial frequencies
%
% Demos for movies.
%   achidemo    - demo of drifting gratings and plaids for achiasmatic stimulus
%   aperdemo    - demo of drifting gratings in various apertures (barber pole illusion)
%   envexpan    - expanding filter envelope over Gaussian noise
%   envrotat    - rotating filter envelope over Gaussian noise
%   flbademo    - flickering bars with a graph of contrast vs. time
%   nnfmdemo    - demo of non-Fourier motion drifting gratings
%   plaidemo    - demo of component gratings and a drifting plaid
%   verddemo    - moving of drifting vernier, two bars
%   verfdemo    - demo of flickering vernier
%   vergdemo    - demo of drifting vernier and phase offsets with gratings
%
% Visual texture generation and control.
%   matshow     - show a matrix of numbers as an image
%   texmake     - make a sample of a specified texture
%   texshow     - show a sample of a specified texture
%   grayadj     - quick and dirty estimate of gray level between black and white
%
% File generation for VSG.
%   cal_bwg_map_demo - create calibration maps for black/white/gray
%   chgcmaps    - change colormap values (look-up tables) in a group of .bmp files
%   ctexchk     - check where a column texture changes
%   ctexfix     - scramble where a column texture changes
%   cuerast     - convert cue and target locations to raster values
%   cuerast_jitt- convert cue and target locations to raster values, with added jitter
%   getmapnums_vsm  - get map numbers from each trial
%   gvsgbwcl    - make B and W textures, graded correlation, for single-stimulus expts
%   gvsgfi      - file generation for J. Tsai's VSG program -- *.vsa, *.vst, *.vsm and *.bmp maps
%   gvsgfisc    - show parameter names and captions
%   gvsgfram    - make frame masks, with various styles
%   gvsggchk    - generic check of parameters
%   gvsgdef     - set up default or initial values of parameters
%   gvsgd1qu    - make quads of 1-d correlated textures with controlled flips and imperfections
%   gvsgeoqu    - make quads of even/odd/random textures with controlled flips and imperfections
%   gvsgeorp    - report on output of gvsgeoqu to ensure that it is close to what is desired
%   gvsgevod    - make even and odd textures with controlled number of flips
%   gvsgfram_test   - test cuerast_jittand gvsgfram
%   gvsgjitt_test   - test cuerast_jitt
%   gvsgmaps    - make maps
%   gvsgmdif    - check whether a set of maps are all different
%   gvsgneed    - define a structure of needed parameters
%   gvsgnew     - get a new set of parameters
%   gvsgobjc    - check object gallery
%   gvsgobjm    - make arrays from object gallery
%   gvsgroun    - round, randomly
%   gvsgsel     - select a set of parameters for a trial
%   gvsgshow    - show information about parameters
%   gvsgsycl    - make 2-fold symmetric textures, graded correlation, for single-stim expts
%   gvsgsyct    - make balanced or unbalanced textures in columns, for use with gvsg2fsy
%   gvsgszct    - make balanced or unbalanced textures in columns, for use with gvsg2fsz
%   gvsg2fin    - convert textures with twofold symmetry to textures with inversion
%   gvsg2fsy    - make balanced maps with twofold symmetry
%   gvsg2fsz    - make balanced maps with two independent twofold symmetries
%   gvsg2fva    - validate number of flips in textures made with gvsgsyct
%   hrb_annulus - make annulus masks for hrb_merge and hrb_map
%   hrb_check   - check options for hrb_merge
%   hrb_getnex  - get options for hrb_merge
%   hrb_map     - make maps for "herringbone" stimuli
%   hrb_map_demo- demonstrate hrb_map for several parameter values and profile types
%   hrb_map_survey_demo     - demonstrate hrb_map for several check sizes and strip sizes
%   hrb_merge   - read vsa, vst files, and make maps for "herringbone" stimuli (two bands of textures)
%   hrbxls_craikplot - plot hrb craik (missing fundamental) data as Weibull functions
%   hrbxls_demo - demonstrates reading and plotting of hrb data
%   hrbxls_maskingplot - plot hrb masking (circular mask) data as Weibull functions
%   hrbxls_merge - merge datasets read from a spreadsheet
%   hrbxls_scalingplot - plot hrb scaling data as Weibull functions
%   hrbxls_survplot - plot hrb survey data as contour maps
%   hrbxls_read - read and analyze a spreadsheet with hrb data
%   hrb_unitprofile - make unit profiles for "herringbone" stimuli (modulators for two bands of textures)
%   makebigcondp_btc - utility routines to get params for makebigconds_btc
%   makebigconds_atten - make a BigConds arrays for a btc adaptation/attention experiment
%   makebigconds_btc - make a BigConds array for 10-coordinate stimuli (for specifying a set of stimuli)
%   makebigconds_demo - make a BigConds array (for specifying a set of stimuli)
%   mpamake_map  - turn a map into micropatterns, optionally scrambled
%   ptbxparse   - parses an array specifying conditions, for Psychophysics Toolbox software of Charlie Chubb
%   ptbxcvt_atglabs - get short and long alpha, theta, gamma labels from texture mode
%   ptbxcvt_batch - convert a batch of data files for Psych Toolbox software, for analysis 
%   ptbxcvt_btcatest - short script to test ptbxcvt_btcdemo
%   ptbxcvt_btccheck - check map layout and texture parameter arrays for conversion for ptbxcvt_btcgetnex
%   ptbxcvt_btcdemo - demonstrate conversion from Psychophysics Toolbox software of Charlie Chubb for 10-d space and texture mode, extends ptbxcvt_demo
%   ptbxcvt_btcgetnex - get map layout and texture parameters for conversion from Psychophysics Toolbox files, 10-param extension
%   ptbxcvt_btclabs - get short and long labels from texture mode, 10-parameter set
%   ptbxcvt_btcmap - extension of ptbxcvt_map to 10-parameter set and textures, called by ptbxcvt_btcmerge and calls btc_genstim
%   ptbxcvt_btcmap_samps  - samples of maps created with ptbxcvt_btcmap with target in each position, 10 params
%   ptbxcvt_btcmerge - merge and convert from Psychophysics Toolbox files to VSG , for 10-parameter set
%   ptbxcvt_check - check map layout and texture parameter arrays for conversion from Psychophysics Toolbox software
%   ptbxcvt_condense - condense one-trial-per-line, 8-column data file to a one line per condition, 5-column file
%   ptbxcvt_demo - demonstrate conversion from Psychophysics Toolbox software of Charlie Chubb
%   ptbxcvt_demo_safe - obsolete but working version to convert Psychophysics Toolbox files to VSG 
%   ptbxcvt_getnex - get map layout and texture parameters for conversion from Psychophysics Toolbox files
%   ptbxcvt_makeray - make a component of a BigConds array that specifies a ray
%   ptbxcvt_merge - merge and convert from Psychophysics Toolbox files to VSG 
%   ptbxcvt_map - create maps for conversion from Psychophysics Toolbox software to VSG, using GenStimFastExec
%   ptbxcvt_map_demo   - demonstrate ptbxcvt_map and cigc_par, with 2x2 correction options
%   ptbxcvt_map_samps  - samples of maps created with ptbxcvt_map with target in each position
%   ptbxcvt_mapsq_demo - demonstrate ptbxcvt_map and cigc_par, square targets and 2x2 correction options
%   ptbxcvt_mtcdemo - extends ptbxcvt_btcdemo to gray levels
%   ptbxcvt_mtcmap  - extension of ptbxcvt_btcmap to gray levels, called by ptbxcvt_mtcmerge and calls mtc_genstim
%   ptbxcvt_mtcmerge - merge and convert from Psychophysics Toolbox files to VSG , for 10-parameter set with gray-level extension
%   ptbxcvt_tpanal - analyzes effect of target position and target structured vs. random
%   ptbxcvt_vdat - read output files from VSG software and convert to inputs for Chubb analysis
%   rep_line_reg - replace the regions and lines of a bitmap with regions of alternating colors
%   sens_fromefit - get the Weibull sensitivity parameter 'a' from an ellipse fit
%   sens_fromefit_test - test sensfromefit
%   shrink_map  - shrink a set of bmp files (reduce check size) and embed in a gray surround
%   shrink_mask - create masks for the shrink_map maps
%   token_map_demo   - replacess black-and-white checks in a map with tokens, and also makes mask maps via min, max, xor
%   token_map   - subroutine for token_map_demo
%   vsaread     - read and parse vsa file, after they are created by gvsgfi
%   vsawrite    - write vsa file (answers and parameters) independent of gvsgfi
%   vsmparse    - parse a vsm file, already read in by vsmvstread
%   vsmvstmod   - modify vsm and vst files, after they are created by gvsgfi
%   vsmvstmod_cue_post_s2   - variant of vsmvstmod
%   vsmvstread  - read vsm and vst files, after they are created by gvsgfi
%   vsmwrite    - write vsm file (map sequences) independent of gvsgfi
%   vstwrite    - write vst file (timing sequences) independent of gvsgfi
%   vswrite_demo- demonstrates vsawrite, vsmwrite, vstwrite
%
% Analysis of psychophysical data files from VSG.
%   cianal      - classification image analysis (data collected from maps made by ptbxcvt)
%   cianal_cv      - canonical variates (and related) calculations for classification image analysis
%   cianal_cvplot  - plot output of cianal_cv
%   cianal_cvset   - create matrices for canonical variates (and related) calculations
%   cianal_demo - demonstrates generalized canonical variates: cianal_cv, cianal_prep, cianal, cianal_plot
%   cianal_down    - downsample derived images for cianal
%   cianal_gcvrr_demo - demonstrates generalized canonical variates and regularized regression, cianal_cv, cianal_rr
%   cianal_getopts - get options for cianal
%   cianal_pen     - apply a ridge term and a nonsmoothness penalty
%   cianal_plot - plot output of cianal
%   cianal_prep - setup the maps and data for cianal
%   cianal_rr      - regularized regression (and related) calculations for classification image analysis
%   cianal_xv_demo  - demonstrate cross-validation of revcor maps
%   cianal_xv       - cross-validation of revcor maps
%   cigc_cosyne_demo   - demonstrates cigc_par, and also flips, for CoSyNe 2006 poster
%   cigc_demo   - demonstrates cigc_ell, and others
%   cigc_ell    - determine classification image parameters via empiric log likelihood
%   cigc_emp    - determine classification image parameters via empiric probability of blocks
%   cigc_llr    - determine classification image parameters via log likelihood ratio for 2x2 blocks
%   cigc_par    - determine classification image parameters via parity within a template
%   cigc_sum    - determine classification image parameters via sum within a template
%   cispr       - simple penalized regression and cross-validation via drop-one for revcor, rr
%   cispr_demo  - demonstrates cispr
%   cixv        - cross-validation calculation for revcor, gcv, rr
%   cixv_demo   - general cross-validation demo
%   compmaps    - extracts *.bmp files from an archive, compares them, and creates sparse matrices
%   ellipsefit_expand - expand parameter vector from active ellipse parameters to all parameters
%   ellipsefit_reduce - expand parameter vector from all ellipse parameters to active parameters
%   errpata_demo  - demonstrate ptbxcvt_errpata
%   extamaps    - extracts *.bmp files from an archive into arrays of s1 and the answer in s2
%   discdemo    - script that illustrates the use of many of these functions
%   discmaps    - shows average maps for right and wrong answers, and their differences
%   gcv_smqf    - make quadratic form for nonsmoothness penalty for generalized canonical variates calculations
%   hrbcianal   - classification image analysis (data collected with hrb* paradigms)
%   hrbcianal_demo - demonstrate classification image analysis for data collected with hrb* paradigms
%   hrbcianal_plot - plot reverse correlation classification image analysis (data collected with hrb* paradigms)
%   hrbcianal_prep - prepare map and response data for classification image analysis for data collected with hrb* paradigms
%   mapubi      - generate a map of unique block indices
%   mdbcirc_anal - analyze results of an experiment using mdbcirc_gen, mdbcirc_trial
%   mdbcirc_fill - fill in (regularize) results of an mdbcirc experiment
%   mdbcirc_mdl_demo - demonstrates modeling of mdb responses
%   mdbcirc_mdl_fit  - fit a model to mdb responses
%   mdbcirc_mdl_get  - get constraint flags for modeling of mdb responses
%   mdbcirc_mdl_nllold - calculate negative log likelihood of a set of mdb responses via AllModels_JV 
%   mdbcirc_mdl_nll   - calculate negative log likelihood of a set of mdb responses, streamlined
%   mdbcirc_mdl_of    - objective function for mdbcirc_mdl_fit
%   mdbcirc_mdl_pnll  - plot negative log likelihood of a set of mdb responses
%   mdbcirc_read - read and analyze results of an experiment using mdbcirc_gen, mdbcirc_trial
%   oppok_demo  - demonstrates reanalysis allowing opposite errors to count as correct
%   paramaps    - parameterizes *.bmp files from an archive into arrays of s1 and the answer in s2
%   ptbx_effwab - find effective Weibull A, B in from a generic ellipse model
%   ptbxcvt_ellipplot - plot contrasts at criterion fraction correct for a set of pie slices
%   ptbxcvt_ellipplot_c - plot contrasts at criterion fraction correct for each ray, and global model
%   ptbxcvt_ellipplot_mtc - plot contrasts at criterion fraction correct for a set of pie slices, after alloy transformation
%   ptbxcvt_ellipplot_mtc_zoom - like ptbxcvt_ellipplot_mtc, but allows for zooming out to see higher thresholds
%   ptbxcvt_errpata - analyze error patterns in a "details" output of data acquired with ptbxcvt files
%   ptbxcvt_findrays - group the two texture parameters into rays, and parameterize them by contrast
%   ptbxcvt_fixedb     - import a fixed-b fit into an unconstrained fit
%   ptbxcvt_fixedb_gen - import a fixed-b or general ellipse fit into an unconstrained fit
%   ptbxcvt_getci  - get confidence limits via trimmed means and resampled values
%   ptbxcvt_oppok  - convert a "details" array into a "BigConds" array, with opposite considered correct
%   ptbxcvt_plot - plot a conds file
%   ptbxcvt_invpieslice - invert Weibull function within a pie slice
%   ptbxcvt_reorder - reorder a 4-column summary file 
%   ptbxcvt_weibeplot - plot ellipse threshholds based on cardinal-axis Weibull fits with uniform beta, minkowski params
%   ptbxcvt_weibrplot - Weibull plots on each cardinal or diagonal axis
%   ptbxcvt_weibsplot - Weibull plots on multiple rays from Chubb Results structure
%   ptbxcvt_weibull - invoke Chubb software to do Weibull fits
%   ptbxcvt_weibull_boot - invoke Chubb software to do bootstrap resamples on Weibull fits
%   readvdat    - read a psychophysical output file (written by J. Tsai's programs) from VSG.
%   readvsm     - read a vsm file to determine which maps are presented on which trials
%   weibull_demo - demonstrates Weibull analysis and plotting of data from Psychophysical Toolbox or VSG
%   weibull_gen_demo - demonstrates Weibull analysis and plotting of data for general ellipse models
%   weibull_multi_demo - demonstrates Weibull analysis for multiple datasets linked by the b-param
%
% Graphical user interface: even, odd, random textures.
%   texuictl    - user interface for textures with hooks to add cues
%   texuicbk    - callback for above
%   texeoctl    - control for even, odd, and random textures
%   texeocbk    - callback for above
%
% Graphical user interface: pilot attention experiments.
%   attuictl    - user interface for textures and cues
%   attuicbk    - callback for above
%   atteoctl    - control for even, odd, and random textures
%   atteocbk    - callback for above
%   attcuctl    - user interface for setting up cues
%   atttactl    - user interface for setting up target
%   attfrzv     - freeze graphical interface values as a structure
%   combtexm    - combine texture and a cue or target mask
%   cuebal      - balances the cues within the enabled locations
%   cuemask     - creates masks for cues
%   cuerast     - turns cue coordinates (as floating values) to raster locations
%   cuesc2vc    - turns cue coordinates as strings to cue coords as floating values
%   cuescoor    - extracts cue coordinates as strings
%   tarmask     - creates masks for targets
%
% Graphical user interface: insight experiments
%
%   insight     - runs the insight experiment
%
% Graphical user interface: utilities.
%   handv2vv    - extract a vector of values from a vector of handles
%   popedimk    - make a popup box from an edit box
%   popedirv    - revise a popup box already made from an edit box
%   pop2vsc     - extract value, string, and choice from a popup box
%   myquest     - customization of QUESTDLG for my purposes

%   Copyright (c) 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026 by J. Victor



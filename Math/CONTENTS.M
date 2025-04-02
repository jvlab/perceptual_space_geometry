% Miscellaneous m-files by J. Victor
%
% Utility.
%   betaincinv_fzero - inverse incomplete beta function using fzero
%   bitrev      - integers in bit-reversed order, any base
%   bootsa      - bootstrap statistics for an array of structures
%   compstruct  - compare two structures (modified from "Compare.m" of Gigis)
%   cootodsu    - convert two sets of coordinates to a matrix of distances (L1, sum)
%   cootodmx    - convert two sets of coordinates to a matrix of distances (Linf, max)
%   cootodsq    - convert two sets of coordinates to a matrix of squared distances
%   dilog       - dilogarithm function, from Didier Clamond, 2006 (downloaded from Mathworks)
%   dunique     - find the unique points in a space, from distances (not from coordinates)
%   edss        - find closest squared distances within a set; efficient memory use
%   ellipcor    - samples from a series with covariance lying on an elliptical shell (also see gnormcor)
%   eqhisted    - find bins boundaries for histogram equalization, new strategy
%   eqhistedu   - find bins boundaries for histogram equalization, original strategy
%   findleastc  - find the least values in a list, specified by count
%   findleastf  - find the least values in a list, specified by fraction
%   fig_movie_demo - demonstrate showing a movie by replacing figs
%   figs_to_ps   - print all open figures to output.ps
%   figs_to_ps2  - print all open figures to output.ps, Ver 2019b and above
%   figs_to_ps_bestfit - print all open figures to output.ps, best fit, Ver 2019b and above
%   figs_to_ps_landscape   - print all open figures to output.ps in landscape mode
%   figs_to_ps_landscape2  - print all open figures to output.ps in landscape mode, ver 2019b and above
%   filldefault - fill in default values into an option structure
%   find_max_rect - find the maximum rectangle in a bitmap
%   flipalign   - flip vectors to maximize alignment
%   fzipexec    - look in common places to find the wzunzip path
%   gammaincinv_fzero - inverse incomplete gamma function using fzero, JV's substitute for matlab's builtin
%   get_colormap_size - get the size of a colormap of a graphic object
%   get_maxsize_bytes - get maximum memory available to matlab (pc and some unix systems)
%   get_key_vals - get numerical values from a string with text keys
%   getinp      - get an input and check ranges, with a default
%   getinpq     - structure-based utility to get an input and check ranges, with a default
%   if_posdefqf - determine whether a set of quadratic forms is positive definite
%   if_triangle - determine whether a putative distance matrix satisfies triangle inequality, and max power to conform
%   int2nary    - convert integers to n-ary arrays
%   isemtystruct - determine whether a varialbe is empty or an amepty structure
%   legend_strings  - ports some of the functionality of previous versions of legend into Matlab V7
%   lognz       - log if positive, zero if nonpositive
%   logz        - log if positive, zero if nonpositive (same as lognz, for compat with David Goldberg)
%   jack        - jackknife statistics for arrays
%   jack_demo   - demonstrates jack and jacksa
%   jacksa      - jackknife statistics for an array of structures
%   jacksa_demo - demonstrates jacksa
%   mean_s      - calculate means of all fields in a structure
%   mindex_inv  - multi-indexing designed for loops
%   mindex_make - create a multi-index
%   mindex_make_check - create a multi-index or check a pre-supplied one
%   minfew      - find the first few smallest values
%   mixedmod    - successive reduction mod an arbitrary sequence of integers
%   mobius      - mobius function
%   nary2int    - convert n-ary sequences to integer
%   nchoosekb   - binomial coefficients, with no complaints for negative, large, or out-of-range args
%   nchoosekpad - array of binomial coefficients, padded, with protection like nchoosekb
%   nchoosek2seq - determine lexicographic position of a sequence of increasing integers
%   nchoosek2seq_2v - vectorized special case of nchoose2kseq for pairs of increasing integers
%   nchoosek2seq_3v - vectorized special case of nchoose2kseq for triplets of increasing integers
%   nchoosek2seq_3vr - vectorized special case of nchoosek2seq, alternate ordering scheme
%   nicesubp    - find an number of rows and columns for a nice subplot array
%   phacolor    - create a colormap suitable for graphing phase
%   pickrand    - pick a random element from an array
%   primroot    - find and test primitive roots (prime or composite modulus)
%   psi_real    - digamma (psi) function for large real arguments
%   ranpaswd    - generate a random password
%   regress_linwt - weighted linear regression
%   phase_reduced - reduced phase (between -pi and pi)
%   pingdemo    - demonstrate pingtest
%   pingtest    - test the network connection via a ping
%   randperm_lim - random permutations with limited shift of each entry
%   randperm_partial - random partial permutations, fast if very few numbers are drawn
%   randround   - round, randomly
%   randseq     - random sequence of integers from a limited range, limiting nearest-neighbor repeats
%   quantize    - represent an array of floats as integer multiples of a float
%   seq2nchoosek - convert lexicographic position to a sequence of increasing integers
%   series61    - some series related to 61
%   setfields   - set multiple fields of a structure
%   shading_patchonly - modify the shading of a 3-d surface, but leave the edges unchanged
%   show_ss     - show scalar values in a structure of a structure
%   singmult    - find singletons and multiple values in non-negative column
%   sizetrim    - size, with trailing 1's deleted
%   std_s       - calculate standard devs of all fields in a structure
%   surf_augvec - augmented versionof surf that also works if data are a vector
%   tablecol2char - convert a table column to a character array, regardless of type or missing data
%   tags_findedges - find edges (dimension 1) between differently-tagged voxels in a volume (multidimensional)
%   tags_findedges_test - test tags_findedges
%   tags_findfaces - find faces (codimension 1) between differently-tagged voxels in a volume (multidimensional)
%   tree_count_demo - count the trees generated from two predecessors, and the shortest path for the one-dimensonal optimal segmentation algorithm
%   tensor_axb_cxd_demo - demonstrate solution of AXB+CXD=Y via inversion of Kronecker tensor product matrices
%   totient     - Euler phi (totient) function and a generalization of it
%   trigamma_demo - trigamma function and estimates of it
%   warplist    - warp a sequence of reals into an equally-spaced sequence from 0 to 1, taking averages when there are repetitions
%   uvdecor     - decorrelate a set of columns, preserving volume
%
% Distributions and statistical tests.
%   MCMCProb     - Dirichlet-prior sampling of a count distribution (from Kennel,Shlens, Abarbanel, Chichilnisky, 2005)
%   auc_gau      - area under ROC curve for Gaussians, one-dimensional, not nec equal variance
%   auc_gau_cues - use auc_gau to show how mean and variance differences combine to influence discriminability
%   auc_gau_demo - demonstrate auc_gau
%   auc_roc_cut  - receiver operating curve, and area under curve, via single cutpoints applied to empirical PDFS
%   auc_var_demo - area under ROC curve for isotropic Gaussian distributions that only differ in variance
%   auc_varc_demo - area under ROC curve for isotropic Gaussian distributions that differ in variance but in contraasting directions
%   circularcorr - circular correlation coefficient and significance via permutation test
%   circularcorrcoef - circular correlation coefficient (from F. Mechler)
%   circularcovar - circular covariance (from F. Mechler)
%   circstats   - basic circular statistics
%   cohenkappa  - Cohen's kappa statistic for inter-rater reliability
%   dirichlet_sample - sample from a multidimensional Dirichlet distribution (downloaded and modified by Mike Repucci)
%   empconfi    - empirical confidence limits from samples
%   esr_demo    - demonstrates t-circ-like approach for checking additivity of VEP components
%   findsegflips - find surrogate-data shuffles consistent with maintaining trials in segments
%   fisherdisc  - Fisher discriminant, with cross-validation and shuffle test
%   fisherdisc_classify  - classify new data based on classifiers built from fisherdisc
%   fisherdisc_def  - set up default options for fisherdisc
%   fisherdisc_demo - test and demo routine for fisherdisc, real EEG spectral data
%   fisherdisc_do   - utility routine for fisherdisc
%   fisherdisc_outclassify  - similar to fisherdisc_classify, but uses halfway and mapbayes, and different results format
%   fisherdisc_test - test and demo routine for fisherdisc, artificial data
%   fdr         - false discovery rate correction
%   gauexpmom   - definite integrals of Gaussian * exponential* power
%   gemvdist    - generate multivariate Gaussian, uniform, or exponential distributions
%   genedist    - generate Gaussian, uniform, or exponential distributions
%   hartigans   - Hartigan dip statistic with bootstrap evaluation of significance (from F. Mechler)
%   hartigan_test - Same as hartigans, but plots into next window
%   hartigans_noplot   - Same as hartigans, no plot (from F. Mechler)
%   hartiganu   - Same as hartigans, but with plotting options (from F. Mechler)
%   hartigan    - Hartigan dip statistic (from F. Mechler)
%   hidden_maxent_demo -- maximum-entropy distribution with a hidden variable
%   hmm_demo    - simple hidden Markov process demo, showing that it is not simple Markov
%   ig2_illus   - graphics for information geometry, two neurons
%   ksstat      - calculate Kolmogorov-Smirnov statistic
%   ksstatb     - calculate Kolmogorov-Smirnov statistic and assess significance via bootstrap
%   loglik_beta - log likelihood of p-values drawn from a beta distribution
%   loglik_beta_demo - demonstrate loglike_beta, and estimates of beta-distribution parameter in symmetric case
%   loglik_beta_demo2 - as in loglik_beta, but also includes point mass
%   loglik_beta_deriv - derivative of log likelihood of p-values drawn from a beta distribution
%   loglik_beta_except - identifies exceptions for log likelihood calculations with beta-function
%   loglik_beta_test - test relative normalization of discrete and continuous part in loglike_beta
%   loglik_bias_demo -  estimates of log likelihoods and biases from Bernoulli data
%   median_diff_cdf - calculate cumulative probability for difference in medians or quantiles
%   median_diff_cdf_locut - as in median_diff_cdf, but only low estimates at cutpoints, and much faster
%   median_diff_cdf_hicut - as in median_diff_cdf, but only high estimates at cutpoints, and much faster
%   median_diff_demo - demonstrate tests for differences in medians
%   median_diff_pdf - calculate probability density for difference in medians or quantiles
%   median_diff_setup - setup for nonparametric tests for differences in medians
%   median_diff_tail - calculate tails of cumulative distribution for difference in medians
%   median_diff_test - test performance for tests for differences in medians
%   merefrac    - maximum-entropy spike count distributions, with hard refractory period
%   me3mix_demo - demonstrate strains for mixtures of 3-neuron distributions
%   me3pgen     - generate probabilities for 3-neuron maxent models
%   me3strain   - strain index for 3-neuron maxent models
%   me3strain_condit_survey - demonstrate strain index for conditional data and unconditional data
%   me3strain_data_anal - compare strains across regions
%   me3strain_data_anal2 - compare strains across regions, alternate calc for pairwise interaction strength
%   me3strain_data_boxplot - demonstrate strain index for data and make boxplot
%   me3strain_data_count - compare strains across regions and counts significances, plots as function of total spike rate
%   me3strain_data_lock - compare lockout and naive strain and neg log likelihood
%   me3strain_data_lock2 - compare lockout and naive strain and U=(neg log likelihood pair/indep)
%   me3strain_data_survey - demonstrate strain index for data
%   me3strain_demo    - Monte Carlo simulation of distribution of strain index for 3-neuron maxent models
%   me3strain_demo2   - Monte Carlo simulations, with multiple presets
%   me3strain_demo3   - Monte Carlo simulations, with multiple presets, alternate plotting
%   me3strain_dr_demo - me3strain_demo with Dirichlet resampling
%   mvnlps            - integral of a multivariate normal Gaussian over an ellipse, Copyright (C) 2013, Alan Genz
%   mvnlps_inv        - inverts mvnlps, i.e., finds radius of an ellipsoid that yields a given probability
%   mvnlps_inv_of     - objective function for mvnlps_inv
%   mvnlps_test       - tests mvnlps
%   odi_demo    - demonstrate orientation difference index for von Mises distributions
%   opt_box_iter - iteration step for multidimensional optimal-segmentation algorithm
%   opt_box_run - omnibus routine to run the optimal segmentation algorithm, multidimensional version
%   opt_box_setup - setup for multidimensional optimal segmentation algorithms
%   opt_box_setup_inds - setup of indices (not data-dependent) for multidimensional optimal segmentation algorithms
%   opt_box_setup_sizes - computes sizes of index arrays in opt_box_setup_inds but does not set them up
%   opt_box_setup_test - test opt_box_setup
%   opt_box_showcuts2d - show 2-d projections of cuts from optimal segmentation algorithm
%   opt_box_showcuts3d - show 3-d projections of cuts from optimal segmentation algorithm
%   opt_box_traceback - recover subdivisions for multidimensinal segmentation algorithm
%   opt_box_t1    - make tables to intialize multidimensional optimal segmentation algorithm
%   opt_covatags2d - compute andshow covariance descriptors of coding regions from general segmentation algorithm
%   opt_cutpoint_ratio  - toy model of optimal cutpoint for information/coding entropy, linear case
%   opt_cutpoint_ratio2 - toy model of optimal cutpoint for information/coding entropy, quadratic case
%   opt_segend_demo - demonstrates finding a segmentation that optimizes information, comparing segment and tethered versions
%   opt_segend_illus - illustrate the tethered routine
%   opt_segend_illus_draw - draw a segment for opt_segend_illus
%   opt_segend_iter - iteration step for algorithm of opt_segment_demo
%   opt_segend_run - omnibus routine to run the one-dimensional optimal-segmentation algorithm tethered to a single end
%   opt_segend_setup - setup for one-dimensional optimal segmentation algorithm tethered to a single end
%   opt_segend_traceback - recover subdivisions for segmentation algorithm, tethered routine
%   opt_segend_t1 - make tables for opt_segment_demo: single-segment partitions tethered to a single end
%   opt_segend_t2 - make tables for opt_segment_demo: two-segment partitions tethered to a single end
%   opt_segment_defopt - set up options structure for finding a segmentation that optimizes information
%   opt_segment_demo - demonstrate algorithm for finding a segmentation that optimizes information
%   opt_segment_iter - iteration step for algorithm of opt_segment_demo
%   opt_segment_kndemo - demonstrate an arithmetic function relevant to opt_segment_demo and opt_segment_iter
%   opt_segment_primes - find primes that are alternatives for the dynamic programming algorithm
%   opt_segment_run - omnibus routine to run the optimal-segmentation algorithm
%   opt_segment_setup - setup for algorithms of opt_segment_demo
%   opt_segment_traceback - recover subdivisions for segmentation algorithm
%   opt_segment_t1 - make tables to initialize segmentation algorithm: single-segment partitions
%   opt_segment_t2 - make tables to initialize segmentation algorithm: two-segment partitions
%   opt_showtags2d - show tags associated with coding regions from general segmentation algorithm
%   opt_tag_defopts - set default options for opt_tag_iter
%   opt_tag_isconnnected - determine whether a coding region is connected, via isconnected or bwlabeln
%   opt_tag_init - de novo initialization for tag-based optimization of coding regions
%   opt_tag_iter - iteration step for tag-based optimization of coding regions
%   opt_tag_run -  run multiple iterations of tag-based optimization of coding regions
%   opt_tag_setup - set up arrays for repeated calls to opt_tag_iter
%   orderstat_dist - distributions for orders statistics
%   orderstat_mom  -  moments for orders statistics
%   orderstat_de -differential entropy (natural base) for orders statistics
%   pbetabayes_compare - Bayesian comparison of probabiliies estimated from finite samples, Dirichlet prior
%   poicount    - Poisson counting distribution
%   poismemoms  - moments of some maxent distributions with Gaussian input, Poisson output
%   poismemoms_demo - demonstrates poimemoms
%   quantile_diff_demo - demonstrate tests for differences in a quantile
%   quantile_diff_setup - setup for nonparametric tests for differences in a quantile
%   quantile_diff_test - test performance for tests for differences in a quantile
%   randpermvar_demo - demonstrates relationship of mse distance to a random permutation, and variance
%   sampcg_demo - sampling from a Gaussian with linear constraints
%   samp_covar_util - estimate covariance matrix from a sample, and plot iso-density ellipses
%   tcirc       - generic t-circ for one population, arguments as real and imaginary parts
%   t2circ      - generic t-circ for one population, arguments as complex, per Drew Hudson
%   t2circ2     - generic t-circ for two populations, arguments as complex, per Drew Hudson
%   umc_gmc_demo - compares uniform and Gaussian Monte Carlo (from David Goldberg)
%
% Contingency tables.
%   fishexac    - Fisher exact test for 2 x 2 tables
%   mwaytres    - create surrogate multidimensional contingency tables with constrained marginals
%   tbstrain    - calculate strain in a 2 x 2 x n contingency table
%   tbtbnres    - create surrogate 2 x 2 x n contingency tables with constrained marginals
%
% Combinations and partitions: enumerations.
%   allcombs    - arrange indistinguishable things into distinguishable boxes
%   lxcm        - find next combination in lexicographic order
%   lxcmdemo    - demo of lxcm
%   multi_boot_groups - bootstrap resamles, among groups
%   multi_shuff_enum - make an array of all ways of ordering counts(1) 1's, counts(2) 2's, etc.
%   multi_shuff_enum_novec - nonvectorized version of multi_shuff_enum, slower but demonstrating logic
%   multi_shuff_groups - permutations that shuffle among groups, optionally preserving a tag
%   multi_shuff_groups_test - tests multi_shuff_groups
%   multi_shuff_mixent - mixing statistics and entropy of shuffles from multi_shuff_groups
%   nallcomb    - number of arrangements of indistinguishable things into distinguishable boxes
%   nlallcmb    - natural log of number of arrangements of indistinguishable things into distinguishable boxes
%   pairings    - enumeration of pairings of distinguishable items into pairs, indepenent of order
%   part        - partitions of a number, written to a file (from Jan Van lent)
%   partit      - partitions of a number
%   xval_configs_make - make configurations for "folds" in a cross-validation
%   xval_confmtx_make - make a confusion matrix, with separate out-of-sample and in-sample coordinates
%
% Partitions and symmetric functions: manipulations.
%   partdemo    - demo of partition routines
%   partdifc    - find the category-count partition
%   partdual    - dual partition
%   partmpys    - multiply two symmetric functions
%   partmunc    - find multinomial coefficient corresponding to a partition
%   partsgec    - parity and number of elements in the symmetric group whose cycle structure is given by the partition
%   partsymc    - number of terms in a symmetric function, with exponents drawn from a partition
%  Note:  see partdoc.txt for additional documentation.
%
% Partitions: format conversions.
%   partcomp    - reformat partition from complete to compact (inverse of partunco)
%   partpack    - reformat partition from full to compact (inverse of partunpa)
%   partunco    - reformat partition from compact to complete (inverse of partcomp)
%   partunpa    - reformat partition from compact to full (inverse of partpack)
%
% Information calculations, binned.
%   biasentd    - demo and comparison of histent, jackent, bootent
%   bayest      - Bayesian estimates of probabilities and their moments
%   binless1    - script to create entropy illustrations for binless.doc
%   binless1cs  - as in binless1, but also includes Chao-Shen estimator
%   binless2    - script to create information illustrations for binless.doc
%   binless3    - script to create information illustrations focusing on many trials, few samples
%   binpoien    - calculate binned entropy of a Poisson process, for a given event rate per bin
%   bootent     - histogram entropy estimates and bootstrap correction
%   bootents    - entropy and bootstrap entropy of a single sample histogram
%   chaoshen    - Chao-Shen entropy estimate
%   channel_model_calcs - plots a function that arises in optimizing gains
%   channel_model_calcs2 - as in channel_model_calcs, but looks over a wider range
%   control_info_gw - a function that arises in the generalized water-filling problem (see control_info*.doc)
%   entdiscpois - entropy of a Poisson process, discretized
%   entdiscpois_demo - demonstrate entropy of a Poisson process, discretized, with various approximations
%   histbent    - histogram entropy estimates: naive, Bayesian, and std. dev.
%   histbi      - general bias estimator (t-p or jackknife) for histogram entropy
%   histent     - histogram entropy estimates: T-P bias calculation 
%   histentd    - demo of histogram entropy asymptotics, # 1 
%   histente    - demo of histogram entropy asymptotics, # 2
%   histentf    - demo of histogram entropy asymptotics: # 3 (used for Bell ms figures)
%   histentg    - demo of histogram entropy asymptotics: # 4
%   histentj    - demo of histogram entropy asymptotics, # 5 (includes jackknife, like # 3)
%   histinfo    - plugin entropy estimate in a histogram (from probabilities)
%   histinfo_nz - plugin entropy estimate in a histogram (from probabilities), if all entries guaranteed >0
%   histcsbi    - Chao-Shen bias for plugin entropy estimate from a histogram (from count)
%   histjabi    - jackknife bias for plugin entropy estimate from a histogram (from count)
%   histtpbi    - Treves-Panzeri bias for plugin entropy estimate from a histogram (from count)
%   jackent     - histogram entropy estimates and jackknife correction
%   kl_markov_demo - demo of Kullback-Leibler divergences for Markov chains
%   legecoef    - coefficients of Legendre polynomials
%   legecosc    - Legendre polynomial of cos(theta), in terms of cos(n*theta) -- recur and nonrecur
%   legeeval    - values of Legendre polynomials
%   legetest    - script to test legecoef, legeval and compare with builtin Legendre
%   mabound     - calculate Ma bound and Strong "Ma bound"
%   mabound_demo  - demonstrate mabound
%   mix_renyi_pwr - demonstrates mixing inequality for power-law measure of concentration, for renyi_bayes_ms
%   nnasympt_local   - subroutine for nnasympt1
%   nnasympt1   - asymptotics of nearest neighbor estimates of entropy: local mean log nearest neighbor
%   nnasympt2   - numerical coefficient in asymptotics of nearest neighbor estimates of entropy
%   nnasympt3   - asymptotics of nearest neighbor estimates of entropy: global av of mean log nearest neighbor
%   nnasympt4   - another try at asymptotics of nearest neighbor estimates of entropy
%   nnmv        - empiric mean and variance for minimum of correlated Gaussian random variables
%   nn2pt       - two-point limit for log nearest neighbor distance, asymptotic in number of dimensions
%   nn3pt       - mean, variance, and covariance of distances between three points drawn from a Gaussian
%   nnkpt       - k-point limit for log nearest neighbor distance, asymptotic in number of dimensions
%   nnrat_demo  - distribution of ratio of nearest neighbor and next-nearest-neighbor
%   nnrat2_demo  - distribution of ratio of multi-nearest-neighbor distances
%   nn_dim_ent  - nearest neighbor distance distribution in n dimensions
%   poiedemo    - demo of poiest, compares naive, T-P, Poisson 
%   poifdemo    - demo of poiest, compares naive, T-P, Bayes, Poisson
%   poient      - compares poiest with naive, T-P, and Bayesian, shows variances
%   poisdemo    - graphical demo of poiest and its asymptotics
%   poiest      - histogram entropy estimates via assumed Poisson distribution
%   poixinfo    - transinformation from Poisson processes
%   rbias2_demo  - asymptotic bias of coding measures, 2 alternatives, superimposed plots
%   rbias2var_demo  - asymptotic bias and variance of coding measures, 2 alternatives, superimposed plots
%   rbias2a_demo - asymptotic bias of coding measures, 2 alternatives, individual plots
%   rbias3_demo  - asymptotic bias of coding measures, 3 alternatives
%   rbiasvar3a_demo - asymptotic bias and variance of coding measures, 3 alternatives, nice ternary plots
%   rbiasvar3b_demo - asymptotic bias and variance of coding measures, 3 alternatives, unlabelled ternary plots
%   renyi       - renyi entropy
%   renyi1_demo - renyi entropy of histograms
%   renyi2_demo - relative renyi entropy (conditional sense) and table info
%   renyi3_demo - extremal renyi entropy (conditional sense) distributions via Lagrange multipliers
%   renyi4_demo - as in renyi3, but better behavior near Shannon case
%   renyi5_demo - relative renyi entropy (K-L sense) and table info
%   renyi_dpi   - search for exceptions to Data Processing Inequality for various Renyi table entropies
%   renyi_mix   - search for exceptions to mixing (convexity) property for various Renyi table entropies
%   renyie_mixr - search for exceptions to mixing (convexity) property of exponential version, with row sums constrained
%   renyi_mixr  - search for exceptions to mixing (convexity) property, with row sums constrained
%   renyikl     - Kullback-Leibler divergence in the Renyi sense
%   renyikl1_demo    - demo of Kullback-Leibler divergence in the Renyi sense
%   renyiklt    - calculates three Renyi table entropies in the K-L sense
%   renyiklt_demo -  tests whether renyiklt satisfies data processing inequality
%   renyit      - calculates three Renyi table entropies in the conditional sense
%   renyite     - calculates exponential version of Renyi table entropies
%   renyit_demo -tests whether renyit satisfies data processing inequality
%   renyitklt_demo - compares renyit and renyiklt
%   tblxbi      - plugin transmitted information estimate, jackknife or Treves-Panzeri (from count)
%   tblxinfo    - plugin transmitted information estimate in a table (from probabilities)
%   tblxinfo_count  - plugin transmitted information estimate in a table (from counts, uses look-up table for speed)
%   tblxinfo_test - test tblxinfo, tblxcinfo_count
%   tblxtpbi    - Treves-Panzeri bias for plugin information estimate in a table (from count)
%   tgicalc     - calculate transmitted generalized information (for renyi_bayes_ms)
%   tgi_demo    - generalized transmitted information for some 2x2 tables (for renyi_bayes_ms)
%   tgi_bias_demo  - generalized transmitted information for 2x2 tables with bias estimates (for renyi_bayes_ms)
%   tgi_motiv_demo  - generalized transmitted information for majority rule (related to renyi_bayes_ms and "error correcting codes"
%   tgi_motiv_demo  - generalized transmitted information for some 2xcontin tables (for renyi_bayes_ms)
%   tgi_motiv_makey - utility for tgi_motiv_demo
%   tgi_surf    - create a surface of transmitted info values for 2x2 tables
%   tgi_3inp_demo   - generalized transmitted information for some 3x3 tables (for renyi_bayes_ms)
%
% Information calculations, unbinned, especially transform strategy.
%   difidemo    - demonstrates difinmin and difinmin_grad
%   difinmin    - as in difinmin, but special calculations for 1-d cases
%   difinmin_old     - differential information , discretely sampled continuous stim space; discrete resposes
%   difinmin_grad    - as in difinmin_old, but calculates gradient also
%   difjdemo    - further tests of difinmin and difinmin_grad
%   difkdemo    - further tests of difinmin, comparing special 1-d calc with general
%   emvdemo     - demonstrates entbinmv, entminmv
%   entest      - entropy ests of specific dists via standard approach and via a binless transform
%   entestbw    - entropy ests of specific dists via standard approach (parametric in bin width) and via a binless transform
%   entestmi    - entropy ests of specific dists via a binless transform and nearest neighbor
%   entestd     - demo of entest and entestmi via entestbw
%   entbin      - entropy estimates of sets of data via binning
%   entbinmv    - entropy estimates of multivariate data via a binning
%   enthist     - histogram entropy (binning already done)
%   entillus    - illustrates change of variables for entropy calculation
%   entint      - entropy estimates of data via transform and intervals
%   entmin      - entropy estimates of data via transform and nearest neighbor
%   entmin_info - mutual information estimate using entmin, 1-d data
%   entminmv    - entropy estimates of multivariate data via transform and nearest neighbor, with jackknife
%   infsptr     - script demo of xinminmv, xinminho, xinbinmv, xinbinho
%   infsptrf    - function version of infsptr, with features for analyzing real data
%   infsubdv    - subdivides a set of spike trains of different classes into approximately equal blocks
%   infs2tls    - converts spike train structure used by infsptr to what is read by tlsread
%   memstrat    - information equivalent of a memory task
%   multicat_demo  - demonstrates ways of analyzing multivariate-to-categorical data in a 2-pixel space
%   nn1d        - nearest neighbor in 1-d
%   nn1d_sorted - nearest neighbor in 1-d, points already sorted
%   nnz1d       - nearest neighbor at nonzero distance in 1-d
%   nsphvol     - volume of an n-sphere
%   tls2infs    - converts spike train structure created by tlsread to what is needed by infsptr
%   unifxfrm    - demonstrates an attempt to find a uniformizing transformation
%   varent_demo - relationship of histogram variance and entropy
%   warptls     - warp a spike train structure to have uniformly spaced times
%   xinbinmv    - calculation of transinformation from multivariate embedding, binned
%   xinbinho    - xinbinmv, but assumes all spike trains embedded in same dimension
%   xincldmv    - calculation of transinformation from geometric-means clustering algorithm, from from dists
%   xincldho    - xincldmv, but assumes all spike trains embedded in same dimension
%   xinclumv    - calculation of transinformation from geometric-means clustering algorithm, from coords
%   xincluho    - xinclumv, but assumes all spike trains embedded in same dimension
%   xinminmv    - calculation of transinformation from spike trains, using unbinned strategy
%   xinminho    - xinminmv, but assumes all spike trains embedded in same dimension
%
% Spike metrics.
%   aronov_test - tests labdist_faster (from T. Kreuz)
%   aranov_para_test - tests labdist_faster_qkpara (from T. Kreuz)
%   clusterc    - wrapper around clusterg, works with non-continguous classes
%   clusterg    - safer version of distclust_fast, works with singleton or empty classes
%   distclust_fast    - calculates confusion matrix from matrix of distances (from Danny Reich)
%   distclust_row     - row-by-row version of distclust_fast, allowing matching to an external set of standards
%   distclust_test    - test distclust_fast, distclust_row, clusterg
%   distclust_fast_demo - compares distclust_fast with metricclust
%   labdist_fast      - calculate spike time distance for labeled trains (from Dmitriy Aronov)
%   labdist_faster    - calculate spike time distance for labeled trains, faster (from Thomas  Kreuz)
%   labdist_faster_qkpara    - extension of labdist_faster to parallel q, k (from Thomas  Kreuz)
%   spkdallq    - calculates spike time metrics for all values of cost by finding critical values of q
%   spkdallqab  - Alex Backer's modification of spkdallq
%   spkdallq_dist        - does full calculation (recursion and final stage) for parallel algorithm, all values of
%   spkdallq_final       - does final stage of parallel algorithm for unlabeled trains, all values
%   spkdallq_recur       - does recursion for parallel algorithm for unlabeled trains, all values of q
%   spkdallqk_dist       - does full calculation (recursion and final stage) for labeled trains, all values of q and k
%   spkdallqk_final      - does final stage of algorithm for labeled trains, all values of q and k
%   spkdallqk_recur      - does recursion for algorithm for labeled trains, all values of q and k
%   spkdallqk_recur      - does recursion for algorithm for labeled trains, all values of q and k
%   spkdallqk_recur_ctc       - version of spkdallqk_recur_zollweg for calling from spkdallqk_recur
%   spkdallqk_recur_zollweg   - optimized version of spkdallqk_recur from John Zollweg of Cornell Theory Center
%   spkdallqk_showcrit   - shows critical values of k and q
%   spkdallqk_trunc      - reduce the results of a recursion for the full train, to results for truncated trains
%   spkdallqk_truncseq   - convert a full recursion for the full train, to results for sequentially truncated trains
%   spkd        - spike time metric between two trains (from Danny Reich)
%   spkd_qpara  - spike time metric between two trains, parallel in q (From Thomas Kreuz)
%   spkd_slide  - spike time metric, allowing for free translatin of one spike train w.r.t. the other
%   spkdl       - spike time metric between two lists of trains (m and c versions) (from Danny Reich)
%   spkdshow    - show spike metric as a function of cost
%   traceminpath  - show the minimal cost path (from Alex Backer)
%   trtodist    - create a distance matrix (spike time distance) from a set of spike trains
%   tstodist    - create a distance matrix (spike time distance) from two sets of spike trains
%
% Computational geometry.
%   cascade_pois   - make a cascaded (multidimensional) Poisson process
%   cascade_pois_demo   - demonstrate cascade_pois
%   circumci    - find a circumcircle
%   circumcircle2d - find a circumcircle using matlab's triangulation and checking for collinearity
%   fibspiral  - finds approximately equally-spaced points on the surface of a hypersphere
%   fit_conic   - fit a conic to points
%   fit_conic_frac - find the minimum-volume ellipse bounding a fraction of points
%   fit_conic_stats - statistics of variance estimation by several methods
%   fit_conic_unpack - utility for fit_conic and fit_conic_unpack
%   hsphere_sample - sample points on the surface of a hyperbolic sphere via several methods
%   hypangle_xr  - included angle in the Poincare disk model; distances via hypdist_xr and hyperbolic law of cosines
%   hypcompose_xr - compose translations in hyperbolic space as distances and angles, and compute parabolic transformations for Poincare disk model
%   hypcompose_xr_tryout - try several ways of composing fractional linear transformations
%   hypdist     - hyperbolic distance in the Klein-Beltrami-Poincare disk model of hyperbolic geometry (geodesics are straight lines)
%   hypdist_uhp - hyperbolic distance in the Poincare upper half plane model, via log of cross-ratio
%   hypdist_xr  - hyperbolic distance in the Poincare disk model, via log of cross-ratio
%   hyperbolic_contin_demo - demonstrates a strategy for approaching a continuous limit
%   hyperbolic_geom_demo - demonstrates equally-spaced circles on hyperbolic plane, and equal sectors of them
%   hyperbolic_tessel - construct triangular tessellation of hyperbolic plane
%   hyperbolic_tessel_demo - demonstrate hyperbolic_tessel
%   hyperbolic_tessel_demo2 - demonstrate hyperbolic_tessel for a pair of superimposed tessellations
%   hyperbolic_tessel_demo3 - demonstrate color symmetry of the [7,14] hyperbolic tessellations
%   hyperbolic_tessel_demo4 - demonstrate a generic three-fold color symmetry even-order hyperbolic tessellations
%   hyperbolic_tree_asym - embedding of asymmetric binary tree in the Poincare disk via application of frac linear transforms
%   hyperbolic_tree_demo - embedding of a binary tree in the Poincare disk via application of frac linear transforms
%   hyperbolic_tree_demo_alt - alternate method for embedding of a tree in the Poincare disk
%   hyperbolic_tree_gen_demo - modular, generalizable version of hyperbolic tree_trans_demo
%   hyperbolic_tree_gen_run - better modular, generalizable version of hyperbolic tree_trans_demo
%   hyperbolic_tree_trans_demo - embedding a tree in Poincare disk via transformations
%   hyperbolic_walkjump_contlim - continuum limit of walk-with-jump trains
%   hyperbolic_walkjump_convhull - show convex hull of some extreme walk-with-jump trains
%   hyperbolic_walkjump_demo - demonstrate walk-with-jump construction of J. DeBlois
%   hyperbolic_walkjump_dists - distances for Poisson walk-with-jump construction of J. DeBlois
%   hyperbolic_walkjump_dists2 - hyperbolic_walkjump_dists but adds special spike trains and better plots of distributions
%   hyperbolic_walkjump_dists3 - hyperbolic_walkjump_dists but also shows distance to mean endpoint
%   hyperbolic_walkjump_lim - demonstrate walk-with-jump construction of J. DeBlois, limiting cases
%   hyperbolic_walkjump_lim2 - as in hyperbolic_walkjump_lim, but adding random trains
%   hyperbolic_walkjump_lim3 - as in hyperbolic_walkjump_lim2, but rendering spike trains with different spike counts via colors
%   hyperbolic_walkjump_split - calculation of probability that two spike trains cannot cross in the future, for walk with jumps
%   hyperbolic_walkjump_trees - draws trees with limiting cases of spike size and speed
%   hyperbolic_walkjump_w - plots the "W(s,phi)" function from HyperbolicGeom_WWJ_notes.docx
%   hypersph_tet_vol - volume of a hyperspherical tetrahedron
%   hypersph_tet_vol_demo - test hypersph_tet_vol
%   hyp_centroid - hyperbolic centroid in Poincare disk model
%   hyp_flt_apply - apply a fractional linear transformation
%   hyp_flt_find - find a fractional linear transformation that maps three points to three other points
%   hyp_geodesic - find the coords of a geodesic between two points in Poincare disk
%   hyp_interp_demo - demonstrate interpolated spike trains
%   hyp_loid2poincare - convert a point on the hyperboloid model to a point in the Poincare disk
%   hyp_poincare2loid - convert a point in the Poincare disk to a point on the hyperboloid model
%   hyp_quad_tri_demo - trigonometry to demonstrate a coincidence of spike trains
%   hyp_reldist_demo - locus of points in H2 at a fixed relative distance or distance ratio from two points
%   hyp_tree_getexclude - get list of strings to exclude
%   hyp_tree_getparams - get parameters for a tree embedding
%   hyp_tree_make: make a tree in the Poincare disk
%   hyp_tree_show: show a tree in the Poincare disk
%   hyp_tree_words: create the words corresponding to a tree
%   mcdsym_cstep - concentration sep of Rousseeuw-van Driessen algorithm used in mcdsym_fastrvd
%   mcdsym_cvn - calculate non-normalized covariance, for mcdsym_fastrvd
%   mcdsym_fastrvd - determine ellipse of minimum covariance determinant containing a set of points
%   mcdsym_setdef - sets defaults for mcdsym_fastrvd
%   minkcirc_demo - circumference of a circle in Minkowski space
%   persp_xform_find - find a perspective (projective) transformation to align two sets of points in n dimensions (n=2: homography)
%   persp_ssqdif_fit - obective function for persp_xform_find
%   persp_ssqdif - compute deviation between a perspective model and target
%   persp_fit - fit a perspective transformatoin with viewpoint known
%   psubspace   - find projection onto a subspace, and closest point in a hyperplane to the origin
%   stellicos   - stellated icosahedra
%   thin_pois   - make a multidimensional thinned Poisson process
%   thin_pois_demo   - demonstrate thin_pois
%   three_line_problem - plot the surface of points from which three lines appear concurrent
%   vornnsum    - statistics of Voronoi regions and nearest-neighbor distances
%   vorstat     - statistics of Voronoi regions of random arrays of points
%   vorstats    - plots results of vorstat
%   wwj_eigs_demo - eigenvalues of operators for walk-with-jumps
%
% Orthogonal synthesis functions and demos.
%   funkern     - synthesis of the Funk kernel from spherical harmonics (Legendre polynomials)
%   synodemo    - demonstrates orthogonal synthesis (Fourier, Hermite,...)
%   synodemo_geom    - like synodemo, but more shapes, and also quasi-analytic calculation
%   synorth     - performs orthogonal synthesis (Fourier, Hermite,...)
%   synorthf    - creates functions for orthogonal synthesis (Fourier, Hermite,...)
%
% Signal processing utilities and demos.
%   afalwt      - iterative affine aligmnent, 1 pc, with missing data (or weights)
%   afalwt_init - initial part of affine alignment, 1 pc, with missing data (i.e., data filled in by mean)
%   afalwt_iter - iterative stage of affine alignment, 1 pc, with missing data
%   afalwt_test - tests affine alignment, 1 pc, possibly with missing data
%   afalwtm      - iterative affine aligmnent, 1 pc, with missing data (or weights)
%   afalwtmg     - greedy iterative affine aligmnent, 1 pc, with missing data (or weights)
%   afalwtmg_test - tests greedy and non-greedy affine alignment, multiple pcs, possibly with missing data
%   afalwtm_init - initial part of affine alignment, multiple pcs, with missing data (i.e., data filled in by mean)
%   afalwtm_iter - iterative stage of affine alignment, multiple pcs, with missing data
%   afalwtm_test - tests affine alignment, multiple pcs, possibly with missing data
%   afalwtm_wca  - weighted centering of affine term for multiple-pc affine alignment, missing data
%   covxform_demo - demonstrate transformations of a covariance matrix
%   dichot_gau_alpha - demonstrate analytic calculation of 4-point statistic alpha for dichotomized Gaussian
%   dichot_gau_quads - demonstrate analytic calculation of quad probabilities for dichotomized Gaussian
%   dichot_gau_triplets - demonstrate analytic calculation of triplet probabilities for dichotomized Gaussian
%   dpri2fc     - convert d-prime to fraction correct
%   dpri2fc_test- tests dpri2fc
%   eyemov_diagram - figure showing effects of eye movements on 1/f and grating
%   eyemov_ptlstab - test calculations for partially stabilized Brownian model
%   eyemov_sac_demo - demo of effects of ramp saccades on input spectrum 
%   eyemov_sac_msparams - get main sequence params for saccades
%   eyemov_sac_qktvar - calculate variable part of q(k,t), for saccades
%   eyemov_spec_cv_test - test calculations of spectra for constant velocity, gamma-process reversal
%   eyemov_spec_fbm -  spectra (q(k,t), q(k,w)) for fractional Brownian motion
%   eyemov_spec_v2oup_demo -  a special function that arises in analysis of the OUP model for eye movements
%   fc2dpri     - convert fraction correct to d-prime
%   fc2info     - convert fraction correct to information
%   fc2tgi      - convert fraction correct to generalized information
%   gaborpxlwts - calculate weights of a Gabor in square pixels
%   gaborpxlwts_demo - demonstrate gaborpxlwts
%   gaborpxlwts_plot - plot results of gaborpxlwts
%   gau         - Gaussian, parametrized by mean and standard dev
%   gaubox      - integral of Gaussian over a multidimensional box, orthogonal to the axes
%   gauint      - cumulative Gaussian, parametrized by mean and standard dev
%   gauf        - folded Gaussian, parametrized by mean and standard dev
%   gaumixinv   - inverse of cumulative of a Gaussian mixture
%   gaumixinv_of- objective function for gaumixinv
%   gaussian_btc_demo - binary image statistics of Gaussian noises
%   info2fc     - convert information to fraction correct; iterative and slow
%   info2dpri   - convert information to d-prime, much faster than fc2info(dpri2fc(.))
%   maxddist    - calculate distribution of maximum of several draws from a discrete distribution
%   multicat_rf1  - pieces for figure 1 (simulated stim/resps) for revised multivariate-categorical paper
%   multicat_rf2  - pieces for figure 2 (2-pixel scatterplots) for revised multivariate-categorical paper
%   prewhiten_demo - demonstrate prewhitening of synthetic blobby images, and calculate image stats
%   prewhiten_mri_demo - demonstrate prewhitening of sample MRI images, and calculate image stats
%   qseqdemo    - demonstrates q-sequences (per Jewett)
%   rocrdemo    - demonstrates ROC curves and related
%   rocrfn      - function used by rocrdemo
%   rfci_demo   - demos several ways (revcor, difference, mseq, GLM) of calculating RF's
%   rfcip_demo  - same as rfci_demo, but 3-d plots and other defaults and labeling
%   rfcj_demo   - demos several ways (revcor, mseq, GLM, Fisher) of calculating RF's
%   rrsim_demo  - regularized regression demo, related to cianal*
%   saccfem_demo - simulations of saccadic eye movement followed by fixational eye movement
%   sos_demo     - musical demo of sum of sinusoids
%   spec_renewdec_demo - power spectrum of process that is decorrelated after one step of a renewal process
%   tgi2fc      - convert generalized information to fraction correct 
%   udasym_demo - demos the asymmetry between detecting increases and decreases in Poisson distributions
%   udasym_simple - simple plots of asymmetry between detecting increases and decreases in Poisson distributions
%   udasyma_demo - an attempted lower bound for KL distance for decreases in Poisson distributions
%   udasymn_demo- adds noise to the above demo, and focuses on the asymptotics of very slow rates
%
% Graphs and their tours, random matrices.
%   connmemw    - find a connectivity matrix for a memory wheel tour
%   gdi_cyclereduce - reduce the cycle structure by exclusive or with singleton overlaps
%   gdi_explore - explore inequalities of distances on graphs
%   gdi_explore_bip_demo - special case of gdi_explore_demo, for bipartitie partitions, general number of nodes
%   gdi_explore_demo - demonstrate gdi_explore
%   gdi_explore_N3bip_demo - special case of gdi_explore_demo, for bipartitie partitions, N(half-nodes)=3 only 
%   gdi_explore_N3spec_demo - special case of gdi_explore_demo, for Sylvester subsets of exhaustive pairings, N(half-nodes)=3 only 
%   gdi_explore_plot_orig - for gdi_explore: plot a graph, pruned graph, and pairings
%   gdi_explore_plot_pruned - for gdi_explore: plot a pruned graph, pairings, partial orderings
%   gdi_ptlorder - partial-ordering utility for gdi_explore
%   gdi_ptlorder_augment - augment a partial ordering of routes by an ad hoc inequality
%   gdi_ptlorder_augment_demo -- demonstrate gdi_ptlorder_augment
%   gdi_ptlorder_routes - partial-ordering of routes on graphs
%   isconnected - determines whether multidimensional region is simply connected
%   isgraphc    - determines whether a matrix corresponds to a connected graph
%   randmtx_eigmom_demo - eigenvalues of random matrices, or with same-sign-in-column constraint
%   shufsegp    - shuffle a tour, as a list of roads, preserving successive pairs, using two pairs
%   shufsegt    - shuffle a tour, as a list of roads, preserving successive pairs, using triples
%   simexp_design - design a graph for the similarity experiment
%   simexp_design_demo - demonstrate simexp_design
%   simexp_design_stats - compute statistics of design generated by simexp_design
%   tourdirg    - find an all-edge round-trip tour of a directed graph
%   toureul     - find an all-edge round-trip tour ("Eulerian") of a symmetric graph
%   tourmemc    - combine and check a standard or multidimensional memory wheel tour
%   tourmemm    - find a multidimensional memory wheel tour
%   tourmemw    - find a memory wheel tour
%
% Bell polynomials.
%   bellds1     - find ptl(b)/ptl(s) at s=1 where b(x,s) is the sth Bell polynomial
%   bellds1d    - demo of asymptotics for bellds1
%   bellnext    - iterative construction of the Bell polynomials
%
% Hermite polynomials and simultaneous band-space-limited functions.
%   bdband      - calculate spatial profile corresponding to a bandlimiting filter
%   bddemo      - demonstrates bdband, bdprofil
%   bddemo_pert - demonstrates bdband, bdprofil with perturbing factor in exponent -- inelegant
%   bdnzdemo    - demonstrates behavior of Wigner kernel (from bddemo) near zero
%   bdprofil    - calculate space-limiting profile (in space) or band-limiting profile (in frequency domain)
%   bdprofsq    - square root of bdprofil and derivs of square roots
%   bdwcline    - calculate Wigner transform of symmetrized b- and s- (self-adjoint) limiting kernel and finds contour lines
%   bdwclinf    - objective function for fzero in bdwcline
%   bdwclinp    - plot Wigner transform of symmetrized (self-adjoint) b- and s- limiting kernel and draws contour lines
%   bdwidth     - interconvert 1/e-width and standard deviation width of a window function
%   bdwigner    - calculate Wigner transform of symmetrized (self-adjoint) b- and s- limiting kernel
%   bdzsqwig    - calculate Wigner transform of "square root" of symmetrized (self-adjoint) b- and s-limiting kernel
%   bsgaup      - parameters for Hermite band- and space-limited functions
%   bsgaurad    - radial component of Hermite band- and space-limited functions in polar separation
%   bsgauval    - values of Hermite band- and space-limited functions, in Cartesian separation
%   bsgdemo     - demonstrates bsgcalc
%   bsgcalc     - full complex calculation of Hermite b- and s- limited functions from polar separation
%   bsgccalc    - calculation of Hermite b- and s- limited functions in Cartesian separation
%   bsgcdemo    - demonstrates bsgccalc
%   bshpcalc    - calculation of orthogonal Hermite b- and s- limited functions in Cartesian separation
%   bsdfcalc    - calculation of differences of adjacent Hermite b- and s- limited functions in polar separation
%   bspicalc    - calculation of Hermite b- and s- limited functions in polar separation with spiral option
%   bspidemo    - demonstrates bspicalc
%   bsms_fig    - create a figure for bsms manuscript of polar and cartesian separations
%   bsms_wig    - create a figure for bsms manuscript of Wigner transforms
%   herdctop    - diagram for heroctop
%   hermcgen    - coefficients of generalized Hermite polynomials
%   hermcoef    - coefficients of Hermite polynomials
%   hermctop    - convert from cartesian to polar separation of generalized Hermite polynomials
%   hermdemo    - demonstrates the Hermite routines
%   hermeval    - values of Hermite polynomials, including simple asymptotic
%   heroctop    - same as hermctop, but assumes orthonormal real polynomials, not monic complex-valued
%   hgeneval    - values of generalized Hermite polynomials
%   hmu1dotx    - dot-products of mu=1 polar Hermite functions with x
%   ica_herm_demo - surrogate entropy measures based on Hermite polynomial expansion and projection on Gausian (per B. Knight)
%   ica_herm_gauproj - projection onto a Gaussian, utilty for ica_herm_demo
%   ica_hermite_notes - asymptotics for integrals of triples of Hermite functions, to estimate variance of Hermite coef estimates
%   ica_hermite_empir - simlar to ica_hermite_notes, but empirical data are standardized prior to estimating Hermite function coefs
%   tripint     - graphs and integrals of triple products of Hermite polynomials and functions, Slepian functions
%
% Special functions.
%   betalnv     - logarithm of beta function
%   ddlgamma    - second derivative of logarithm of gamma for real arguments
%   dlgamma     - psi (derivative of logarithm of gamma) for real arguments
%   dpssdemo    - examples of discrete prolate spheroidal sequences and wavefunctions
%   dpsswdemo   - shows Fourier windows associated with discrete prolate spheroidal sequences
%   eulermas    - Euler-Mascheroni constant and related
%   gammaiz     - incomplete gamma for complex arguments
%   gammalnz    - logarithm of gamma for complex arguments
%   gammaz      - gamma for complex arguments, extension of Matlab's gamma
%   gammazz     - gamma for complex arguments, another extension of Matlab's gamma
%   igdens      - density of inverse Gaussian distribution
%   igmom       - moments of inverse Gaussian distribution
%   nestedmean_proj - projection matrix to remove mean of subset in a full factorial design
%   nestedmean_ssdecomp - sum-of-squares decomposition in a full factorial design
%   psi         - psi (derivative of logarithm of gamma) for complex arguments
%   psi_jv      - rnamed version of psi to prevent namespace conflicts in v2019b
%   raywalkd    - probability density for distance moved in a Rayleigh random walk
%   radproj     - projection of a circularly symmetric function onto a line
%   sinc        - sinc (sin(x)/x) function
%   sincdemo    - demonstrates various windows (tapers) for fft
%   slitconv    - n-fold convolution of a slit
%   zrpsfda     - zeroth order radial prolate spheroidal function, discrete approximation
%
% Orthogonal matrices, matrix functions and group representations.
%   axistorot   - 3-d rotation matrix from axis and rotation angle
%   ccvorth     - finds orthogonal real vectors that have same span as z and conj(z)
%   cosdih      - cosine of dihedral angle between two hyperplanes
%   extorthb    - extend a vector to an orthogonal, or orthonormal, basis
%   extorthb_gen- extend one or more orthonormal vectors to an orthogonal, or orthonormal, basis
%   grmscmdt    - Gram-Schmidt procedure
%   infrot      - infinitesimal generator of a rotation
%   kron_demo   - two ways of transforming a matrix by pre- and post-multipling
%   liebrack    - Lie bracket of two square matrices
%   outprod     - outer product of two arrays, can be multidimensional (see also kron)
%   outprod_contract - outer product of two arrays, contracted over one or more dimensions, can be multidimensional
%   planesph    - polar projection from (n-1)-dim plane to surface of n-dim sphere, centered anywhere
%   planesph_inv- inverse of planesph
%   randorth    - random orthogonal matrix, rotations uniformly distributed in SO(d)
%   randorthu   - random orthogonal matrix, rotations uniformly distributed in SO(d)
%   randorthu_gen  - random orthogonal matrix with one or more directions fixed
%   rangdemo    - distribution of random angles, including generalized dihedral angles, in n-space
%   rot33to4    - maps two 3-d rotations to one 4-d rotation via the wedge product
%   rot4to33    - maps one 4-d rotation to two 3-d rotations via the wedge product
%   rotademo    - demo of rotation routines
%   rottoaxis   - axis and rotation angle for a 3-d rotation matrix
%   so3d5rep    - find the 5-dimensional representation of a 3-dimensional rotation matrix
%   so4twoang   - find the two canonical rotations of a 4-dimensional rotation matrix
%   so5twoang   - find the two canonical rotations of a 5-dimensional rotation matrix
%   sym_tensor_rep - find the symmetrized tensor product representation of a matrix (need not be orthogonal)
%   wedgerep    - wedge representation
%   xdet        - find the cross-determinant of two 2x2 matrices
%
% Canonical forms of matrices.
%   critlmut    - calculates criteria for whether a matrix is near upper-triangular
%   matsosmi    - iteratively rotate a matrix to minimize selected elements
%   matsos      - sum of squares of selected matrix elements, with optional transformation
%   mlindemo    - explore invariants of multilinear forms, to determine when a set of matrices is triangularizable
%   mliodemo    - uses gceroots and gcenull to determine when a set of matrices is triangularizable
%   randmeis    - random matrix with specified number of complex eigenvalues
%   rdr2sl      - inverse of sl2rdr
%   rdr2sl_demo - demonstrate rdr2sl and linear independence of triple products
%   sl2rdr      - write a 2x2 matrix, determinant 1, as a product of a rotation, a stretch, and a rotation
%   tri2all     - find rotation angle to minimize left lower-triangular element of a 2x2 matrix under adjoint transformation
%   tri2best    - apply rotation to minimize left lower-triangular element of a 2x2 matrix under adjoint transformation
%   triedemo    - demons/explores aspects of triangularized forms of matrices #1
%   trifdemo    - demons/explores aspects of triangularized forms of matrices #2
%
% Multilinear autoregression.
%   gcegetc     - conforms a coefficient array to shape [r r k]
%   gcenull     - nullspace for one root of generalized characteristic equation
%   gceroots    - roots of generalized characteristic equation
%   gnormcor    - samples from a multivariate Gaussian with arbitrary covariance
%   gnormcor_lc - samples from a multivariate Gaussian with arbitrary covariance and linear constraints
%
% Multidimensional scaling.
%   domds       - do a multidimensional scaling from a distance matrix
%   hdiskmds    - do a multidimensional scaling with points in a particular geometry
%   hdiskmex    - find minimum "fixing" exponent for the disk geometry
%   procrustes_compat - rename transformation fields of procrustes_consensus to that of procrustes.m
%   procrustes_consensus - iterative computation of consensus of several datasets by orthogonal Procrustes analysis
%   procrustes_conesnsus_ptl_test - test procrustes_consensus for datasets with only partial overlap
%   procrustes_conesnsus_test - test procrustes_consensus
%   procrustes_ref_devs - compute deviation of individual datasets from a reference or consensus
%   procrustes_varignon_demo - check Varignon theorem for consensus of consensus
%   procrustes_varignon_util - consensus of two datasets via standard Procrustes and midpoint
%   riemds      - MDS of points in a Riemannian manifold or Minkowski space with a "fixing" exponent
%   riemds_a    - analytic calculation of eigenvalues for riemds, initial method
%   riemds_b    - analytic calculation of eigenvalues for riemds, telescoping
%   riemds_c    - analytic calculation of eigenvalues for riemds, with series for S(a,n,b) summed
%   riemds_calc - does some of the calcs in riemds via a group-theoretic approach
%   riemds_d    - analytic calculation of eigenvalues for riemds, analytic conjecture for all a
%   riemds_e    - calculation of a sum for riemds in hypergeometric form
%   riemds_dmp2cos - calculates D<m><p/2>, orthogonal polynomials related to Gegenbauer/ultrasphericals
%   riemds_dorth- demonstrate orthogonality of the "D" polynomials, related to Gegenbauer/ultraspherical
%   riemds_fpa  - calculate the function f<p>(a), used in the "D" polynomials
%   riemds_hunt - iterate riemds to find highest "fixing" exponent
%   riemds_ica  - does integrals for riemds generalized eigenvalue calc
%   riemds_int  - does integrals in riemds_calc
%   riemds_mink - does integrals related to riemds for the hyperbolic case
%   riemds_sanb - calculation of a sum for riemds
%   riemds_sanbr- demonstrate a recursion relationship for riemds
%   riemds_tanb - calculation of a sum for riemds
%   riemds_tanbr- demonstrate a recursion relationship for riemds
%   riemds_sph_catalan -does integrals related to mds for the sphere-surface case related to Catalan numbers
%   riemds_subdim_show - displays results of several riemds runs, saved in riemds_subdim.mat
%   riemds_vmpc - calculate the "V" function, used in the "D" polynomials
%   sinidemo    - some definite integrals involving sines
%   umtpdemo    - modified ultraspherical polynomials, as trigonometric polynomials
%
% Renewal processes.
%   altcov      -  demo of auto and cross-covariances of alternating renewal process
%   altcovs     -  a gallery of altcov plots, different scales, rate normalization
%   altcovsa    -  a gallery of altcov plots, emphasizing the asymmetry
%   fit_gamma - fit parameters of a gamma process
%   fit_gamma_demo: demonstrate fit_gamma
%   gamma_process - generate samples of a gamma process
%   gamma_process_entropy - entropy of a gamma process
%   gamrate     -  z-distrib renewal process, specified by order and rate
%   gamrenew    -  z-distrib renewal process, specified by order and time const
%   stibeh_analyze  - analyze stimulation-behavior data
%   stibeh_demo  - demonstrate stimulation-behavior fitting package
%   stibeh_demo2 - same as stibeh_demo, but only fits constant model and ramp
%   stibeh_fit   - fit a model to stimulation-behavior data -- fmin with various starting points
%   stibeh_fmin  - carry out fmincon to model to stimulation-behavior data
%   stibeh_group  - group a behavioral dataset, preparing for analysis for presence of runs
%   stibeh_mersti- merge nearby stimulation periods in a stimulation-behavior dataset
%   stibeh_nll   - calculate negative log likelihood
%   stibeh_plot  - plot results of stimulation-behavior fitting package
%   stibeh_probs - calculate probabilities from a model for stimulation-behavior fitting package
%   stibeh_report - report analysis of stimulus-behavior fitting package
%   stibeh_runs_demo - demonstrate analysis for presence of runs
%   stibeh_runs  - analyze for presence of runs
%   stibeh_rvals - calculate convolution of stimulation train for stimulus-behavior fitting package
%   stibeh_synth - synthesize datasets for stimulation-behavior fitting package
%   stibeh_timlim  - time-limit a stimulation-behavior dataset
%   stibeh_timrev  - time-reverse a stimulation-behavior dataset
%   stibeh_xtoparams - convert a vector of parameters to a structure
%
% Plotting.
%   circle      - draw a circle (used for circstatsplot, from Ferenc Mechler)
%   circles     - draw circles with specified centers (from Ferenc Mechler)
%   circstatsplot - draw a circular histogram (from Ferenc Mechler)
%   conetest    - shows how to use coneplot and assign value-dependent colors to the cones
%   decisrule0 - diagram of a decision rule (see renyi_bayes_ms1*.doc) 
%   decisrule1 - diagram of range of acheivable 2x2 decision rules (see renyi_bayes_ms1*.doc) 
%   decisrule_line - utility for decisrule0
%   decisrule_patch - utility for decisrule0
%   img_gen_cartoon - cartoon diagrams of several image-generation methods (natural, xray, scintigram)
%   myerrorbar  - error bar plotting, wider error bars
%   myerrorbar_xy  - error bar plotting with x and y error bars; adjustable width error bars
%   mypolar     - polar plotting with more flexibility than Matlab's
%   stars_demo  - draw some stars, range of points, range of skips
%   stars_demo2 - draw some stars, fixed number of points, range of skips
%   surf2stl    - convert surfaces to stl files
%   surf2stl_demo - demonstrate surf2stl with a monkey saddle
%   torus       - draw a torus, analogous to Matlab's SPHERE.
%   torus_angc  - draw a torus from angular coordinates
%   torus_demo  - demonstrate torus draing and plotting
%   wedge       - draw a wedge (used for circstatsplot, from Ferenc Mechler)
%
% Dynamical systems -- ODE's.
%   qualode     - make a nice plot of trajectories of an ODE.  Any of the functions
%                 below can serve as fodder
%   coullet     - third-order ODE of Coullet, Tresser, and Arneodo
%   couplvdp    - slow-time approximation for two coupled van der Pol oscillators
%   glider      - idealized glider
%   lotke       - Lotke-Volterra predator-prey system
%   monkeysa    - monkey saddle vector field
%   rotpen      - rotating pendulum
%   sample3     - sample third-order ODE
%   vanderpol   - van der Pol oscillator
%
%  Information calculations and adaptation
%   ita_default - set up default values for information calculations
%   ita_firingrate - firing rate from tuning params
%   ita_firingrate_demo - demonstrate ita_firingrate for von Mises
%   ita_firingrate_test - demonstrate ita_firingrate for all types of curves
%   ita_mi_direct - computer mutual information by direct calculation of stimulus-response distribution
%   ita_mi_direct2 - compute mutual information by direct calculation of stimulus-response distribution, two neurons
%   ita_mi_directn - compute mutual information by direct calculation of stimulus-response distribution, n neurons
%   ita_mi_montecarlo - compute mutual information by Monte Carlo method, n neurons
%   ita_mi_test - test ita_mi_direct, ita_mi_montecarlo
%   ita_pcount - probability of a number of spikes (Poisson) from tuning (von Mises params) and response duration 
%   ita_tuning_default - set default tuning parameters
%   poiss_info_calcs_compare - time-dependence of information for Poisson neuron, von Mises or triangle tuning
%   poiss_info_calcs_compar2 - time-dependence of information for Poisson neuron, roof tuning
%   poiss_info_calcs_demo - asymptotic estimates of information for Poisson neuron, von Mises or triangle tuning
%   poiss_info_calcs_test - test asymptotic estimates of information for Poisson neuron, von Mises or triangle tuning
%
%  Toy control simulations
%   contoym_asn - active sensing noise
%   contoym_asnexp - expected value of active sensing noise
%   contoym_default - set up default params
%   contoym_sim_demo  - simulation demonstration
%   contoym_sim_demo2 - simulation demonstration, graphics tweaked for presentation
%   contoym_sim_do   - simulate one or more control runs, with a single parameter set
%   contoym_sim_lag_demo  - simulation demonstration with time lags
%   contoym_sim_test - tests simple controller simulations
%
% Protocerebral Bridge simulation, and related fly navigation
%   pcb_pcasim_demo - simulation of PCA analysis of protocerebral bridge activity
%   heading_ambig_demo - demonstrate heading ambiguity
%
% Ultrametric
%   ultra_spike_demo - demonstrations related to hybrid spike metric/ultrametric
%   ultra_spike_single - calculate hybrid spike metric/ultrametric for single-spike trains
%
% Flow statistics
%   flow_spec_ftdemo -- verify Fourier transformation of diffusion-with-advection 
%   flws_a0make - create time-domain descriptor of flow effects on first-order statistics
%   flws_b0make - create time-domain descriptor of flow effects on second-order statistics
%   flws_bdiffmake - average time-domain descriptor over mean transport times
%   flws_b0demo - demonstrate flws_b0make and related
%   flws_bwigmake - Fourier transform of bdiff
%   flws_bwig_demo - demonstrate a family of bwigs
%   flws_defopts - define options
%   flws_getparams - get parameters foor source and sensor
%   flws_normalizations - analyze normalizations, trajectory details
%   flws_ptsim_a0b0test - test flws_ptsim_a0, flws_ptsim_b0
%   flws_ptsim_a0 - flws_a0make for particle tracking 
%   flws_ptsim_b0 - flws_b0make for particle tracking
%   flws_ptsim_sparse_test - test a sparse-matrix version of flws_ptsim_[a0|b0]

%   Copyright (c) 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025 by J. Victor



function opts_use=psg_defopts(opts)
%opts_use=psg_defopts(opts) sets up default options for perceptual space geometry experiments
%
% opts: optional overrides, may be omitted
%   also may be an array created by psg_sessconfig_make, in which case
%   opts.cond_nstims, opts.cond_ncompares, and opts.cond_nsess are taken from sess
%
% opts_use: full options structure
%
%  17Nov22: add setseq: options for stimulus sets across sessions
%  20Nov22: add typeno_zpad
%  31Dec22: add option of using a session array as input
%  20Feb23: add opts.if_eachsess
%  31May23: add cond_nstims_toreplace
%
%  See also  FILLDEFAULT, PSG_SESSONFIG_MAKE, PSG_COND_WRITE, PSG_SETUP_DEMO, PSG_COND_CREATE,
%   PSG_COND_WRITE, PSG_COND_CREATE.
%
if (nargin<1)
    opts=struct;
elseif ~isstruct(opts)
    sess=opts;
    opts=struct;
    opts.cond_nstims=length(unique(sess(:)));
    opts.cond_ncompares=size(sess,2)-1;
    opts.cond_nsess=size(sess,3);
end
%
%logging and calculation
opts=filldefault(opts,'if_log',0); %1 to log
opts=filldefault(opts,'if_cumulative',0); %whether to compute statistics for each session cumulatively
opts=filldefault(opts,'if_eachsess',1); %whether to compute statistics separately for each session
%
%session parameters
opts=filldefault(opts,'cond_nstims',25);
opts=filldefault(opts,'cond_nstims_toreplace',0); %number of stimuli to replace by lower-numbered stims
opts=filldefault(opts,'cond_ncompares',8);
opts=filldefault(opts,'cond_novlp',2);
opts=filldefault(opts,'cond_nsess',10);
opts=filldefault(opts,'refseq',1); %randomization method for reference stimuli shared across contexts
opts=filldefault(opts,'refseq_labels',{'random','ordered','frozen randomization'});
opts=filldefault(opts,'setseq',1); %how stimulus sets are changed across sessions
opts=filldefault(opts,'setseq_labels',{'unique','same_sets','same_sets_same_positions','same_positions'});
%
%options related to stimulus example reuse, mostly for psg_cond_create
opts=filldefault(opts,'example_infix_mode',1); %whether to use different examples of each stimulus across sessions, or within sessions, or not at all
opts=filldefault(opts,'example_infix_labels',{'different examples across all sessions','different examples within session','single example','single example, no infix'});
opts=filldefault(opts,'example_infix_string','_'); %separator between stimulus name and example number
opts=filldefault(opts,'example_infix_zpad',3); %number of digits to zero-pad 
%
%file name creation options
opts=filldefault(opts,'sess_zpad',2); %zero-padding in session name
opts=filldefault(opts,'stim_filetype','png'); %file type
opts=filldefault(opts,'typeno_zpad',2); %padding for type number (unique ID for stimulus, 1 to cond_nstims)
%
opts_use=opts;
return

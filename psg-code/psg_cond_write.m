function [msgs,opts_used]=psg_cond_write(filename,textcell,opts)
% [msgs,opts_used]=psg_cond_write(filename,textcell,opts) writes a configuration file for the psg experiments
%
% filename: file name (and path), .csv appended if not given
% textcell: cell array of stimulus file names, without the .png extension
% opts: options
%  opts.if_log: 1 to log messages (defaults to 0)
%  opts.cond_write_prefix: a text string to precede every element of textcell
%
% msgs: error messages, empty if successful
% opts_used: options used
%
% 01Jun23: add opts.cond_write_prefix
%
%   See also: PSG_SETUP_DEMO, PSG_SESSCONFIG_MAKE, PSG_DEFOPTS, CSVWRITE, FOPEN, FPRINTF, FCLOSE.
%
if (nargin<3)
    opts=struct;
end
opts=psg_defopts(opts);
opts=filldefault(opts,'cond_write_prefix',[]);
opts_used=opts;
msgs=[];
%
ncompares=opts.cond_ncompares;
ncols=size(textcell,2);
if ncols~=(ncompares+1) %one reference + ncompares comparisons
    error(sprintf('Number of columns in textcell (%2.0f) incompatible.  Should be one more than number of comparison stimuli (%2.0f)',...
        ncols,ncompares));
end
ntrials=size(textcell,1);
cond_data=cell(size(textcell,1)+1,ncompares+1);
%create header row
cond_data{1,1}='ref';
for icompare=1:ncompares
    cond_data{1,1+icompare}=cat(2,'stim',sprintf('%1.0f',icompare));
end
%copy in data
for itrial=1:ntrials
    for icol=1:ncols
        cond_data{1+itrial,icol}=cat(2,opts.cond_write_prefix,textcell{itrial,icol});
    end
end
%
msgs=text_csvwrite(filename,cond_data,opts);
return

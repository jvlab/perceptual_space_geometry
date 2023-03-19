function [msgs,opts_used]=text_csvwrite(filename,textcell,opts)
% [msgs,opts_used]=text_csvwrite(filename,textcell,opts) writes a csv file
% from an cell array of text
%
% note that matlab's csvwrite only handles numeric values
%
% filename: file name (and path), .csv appended if not given
% textcell: cell array of data to be written, should be text only (i.e., numeric quantities already converted to ascii)
% opts: options
%  opts.if_log: 1 to log messages
%
% msgs: error messages, empty if successful
% opts_used: options used
%
%   See also: PSG_COND_WRITE, CSVWRITE, FOPEN, FPRINTF, FCLOSE.
%

if (nargin<=2)
    opts=[];
end
opts=filldefault(opts,'if_log',0);
%
opts_used=opts;
msgs=[];
filename=cat(2,filename,'.csv');
filename=strrep(filename,'.csv.csv','.csv');
%
[fid,fopen_msg]=fopen(filename,'wt'); %open for writing in text mode
if (fid<0)
    msgs=strvcat(msgs,fopen_msg);
    if (opts.if_log)
        disp(sprintf('error creating %s.',filename));
        disp(msgs);
    end
    return
end%
msgs=strvcat(msgs,sprintf('file %s open.',filename));
nlines=size(textcell,1);
for iline=1:nlines
    line_out=[];
    for icol=1:size(textcell,2)
        line_out=cat(2,line_out,textcell{iline,icol},',');
    end
    line_out=line_out(1:end-1);
    fprintf(fid,'%s\n',line_out);
end
msgs=strvcat(msgs,sprintf('file %s: %5.0f lines written.',filename,nlines));
%
status=fclose(fid);
if (status<0)
    fclose_msg=sprintf('error closing %s',filename);
    msgs=strvcat(msgs,fclose_msg);
    if (opts.if_log)
        disp(sprintf('error creating %s.',filename));
        disp(msgs);
    end
    return
end
msgs=strvcat(msgs,sprintf('file %s closed.',filename));
if opts.if_log
    disp(msgs);
end
return


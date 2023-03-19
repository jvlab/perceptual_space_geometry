function psg_stat_tally(caption,opts,s,name_zero,name_more,name_max)
% psg_stat_tally(caption,opts,s,name_zero,name_more,name_max) is a utility to display tallies of statistics
%
% caption: title string to display
% opts: opts strucure with cond_nsess, if_cumulative, if_eachsess
% s: field with data
% name_zero: data field name for a tally of 0
% name_more: data field name with a tally > 0
% name_max: data field with maximum of tallies
%
%  See also: PSG_UMI_STATS, PSG_TRIAD_STATS, PSG_QUAD_STATS, PSG_TENT_STATS.
%
disp(caption);
disp(cat(2,'         ',sprintf(' %5.0f',[0:s.combined.(name_max)])));
if (opts.if_eachsess)
    for isess=1:opts.cond_nsess
        disp(cat(2,sprintf('sess %4.0f',isess),sprintf(' %5.0f',s.eachsess.(name_zero)(1,isess),s.eachsess.(name_more){1,isess})));
    end
end
if (opts.if_cumulative)
    for isess=1:opts.cond_nsess
        disp(cat(2,sprintf('sess 1-%2.0f',isess),sprintf(' %5.0f',s.cumulative.(name_zero)(1,isess),s.cumulative.(name_more){1,isess})));
    end
end
disp(cat(2,' combined',sprintf(' %5.0f',s.combined.(name_zero),s.combined.(name_more))));
return


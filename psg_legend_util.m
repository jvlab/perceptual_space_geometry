function [legend_mode_used,strings_used]=psg_legend_util(legend_mode,strings)
% [legend_mode_used,strings_used]=psg_legend_util(legend_mode,strings): utility to change the labels in a legend created by psg_visualize
%
% if no arguments, then gets inputs from console
% legend_mode: 1 to label based on rays, 2 based on datasets, 3 for connections
% strings: cell array of new labels, must have same length as the current legend entries
%
% legend_mode_used: legend mode used
% strings_used: strings_used
%
% See also:  PSG_VISUALIZE, PSG_PLOTCOORDS, PSG_CONSENSUS_DEMO, PSG_LEGEND_KEEP.
%
if (nargin==2)
    if_auto=1;
else
    if_auto=0;
end
h_legend=get(gca,'Legend');
if isempty(h_legend)
    disp('current axis has no legend')
    legend_mode_used=0;
    strings_used=cell(0);
    return
else
    hc=get(gca,'Children');
    hc_keeps=psg_legend_keep(hc);
    hc_types=fieldnames(hc_keeps);
    %
    legend_strings=get(h_legend,'String');
    ifn_list=[];
    for ifn=1:length(hc_types)
        fn=hc_types{ifn};
        if length(legend_strings)==length(hc_keeps.(fn))
            ifn_list(end+1)=ifn;
            if ~if_auto
                disp(sprintf(' mode %1.0f ->%s',ifn,fn))
                for k=1:length(hc_keeps.(fn))
                    disp(sprintf('legend label %2.0f: %s',k,legend_strings{k}))
                end
            end
        end
    end
    if if_auto
        ifn_choice=legend_mode;
    else
        if length(ifn_list)>1
            ifn_choice=getinp('mode to use','d',[min(ifn_list) max(ifn_list)],min(ifn_list));
        else
            ifn_choice=ifn_list;
        end
        fn=hc_types{ifn_choice};
    end
    legend_mode_used=ifn_choice;
    if if_auto
        legend_strings=strings;
    else
        for k=1:length(hc_keeps.(fn))
            legend_strings{k}=getinp(sprintf('string to replace label %s',legend_strings{k}),'s',[],legend_strings{k});
        end
    end
    set(h_legend,'String',legend_strings);
    strings_used=legend_strings;
    return
end



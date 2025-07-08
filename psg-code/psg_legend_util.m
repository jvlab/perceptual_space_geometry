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
% See also:  PSG_VISUALIZE, PSG_PLOTCOORDS, PSG_CONSENSUS_DEMO.
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
    tags=cell(length(hc),1);
    for ich=1:length(hc)
        tags{ich}=get(hc(ich),'Tag');
    end
    hc_ds1=find(contains(tags,'ds 1'));
    hc_rays=find(contains(tags,'ray'));
    hc_sign=find(contains(tags,'signed'));
    hc_conn=find(contains(tags,'connection'));
    hc_c1=find(contains(tags,'connect   1'));
    hc_p1=find(contains(tags,'point   1'));
    hc_s1=find(contains(tags,'set  1'));
    hc_replot=find(contains(tags,'replot')); %when connections between datasets are present and if_norays=1
    hc_ds_any=find(contains(tags,'ds ')); %any dataset
    %
    %several options for where the labels might be
    hc_keep=struct;
    hc_keep.rays=intersect(intersect(hc_ds1,hc_rays),hc_sign); %label based on rays
    hc_keep.norays=intersect(setdiff(hc_ds_any,hc_replot),hc_s1);
    hc_keep_conn=intersect(intersect(hc_c1,hc_p1),hc_conn);
    hc_keep.connect_only=intersect(hc_keep_conn,hc_conn);
    %
    hc_types=fieldnames(hc_keep);
    %
    legend_strings=get(h_legend,'String');
    ifn_list=[];
    for ifn=1:length(hc_types)
        fn=hc_types{ifn};
        if length(legend_strings)==length(hc_keep.(fn))
            ifn_list(end+1)=ifn;
            if ~if_auto
                disp(sprintf(' mode %1.0f ->%s',ifn,fn))
                for k=1:length(hc_keep.(fn))
                    disp(sprintf('legend label %2.0f: %s',k,legend_strings{k}))
                end
            end
        end
    end
    if if_auto
        ifn_choice=legend_mode;
    else
        if length(ifn_list)>1
            ifn_choice=getinp('mode to use','d',[min(ifn_choice) max(ifn_choice)],min(ifn_choice));
        else
            ifn_choice=ifn_list;
        end
        fn=hc_types{ifn_choice};
    end
    legend_mode_used=ifn_choice;
    if if_auto
        legend_strings=strings;
    else
        for k=1:length(hc_keep.(fn))
            legend_strings{k}=getinp(sprintf('string to replace label %1.0f',k),'s',[],legend_strings{k});
        end
    end
    set(h_legend,'String',legend_strings);
    strings_used=legend_strings;
    return
end



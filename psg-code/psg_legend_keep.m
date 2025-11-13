function hc_keeps=psg_legend_keep(hc)
% hc_keeps=psg_legend_keep(hc) is a utility to determine a subset of
% children to be used in a legend, based on their tags
%
% Once a field of hc_keeps is selected, it can be used to create a legend with
% only those graphical elements, e.g., legend(hc_keeps.rays)
%
% hc: the children of an axis, typically created by psg_visualize
% hc_keeps: a structure
%   hc_keeps.rays: handles of the children to be kept in which the rays are labeled
%   hc_keeps.norays: handles of the children to be kept in which datasets are plotted as individual points
%   hc_keeps.connect_only: handles of the children to be kept in only connections between datasets are plotted
%
% 12Nov25: add hc_keeps.ds
%
% See also:  PSG_VISUALIZE, PSG_PLOTCOORDS, PSG_CONSENSUS_DEMO, PSG_LEGEND_KEEP.

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
hc_keeps=struct;
hc_keeps.rays=intersect(intersect(hc_ds1,hc_rays),hc_sign); %label based on rays
hc_keeps.norays=intersect(setdiff(hc_ds_any,hc_replot),hc_s1);
hc_keep_conn=intersect(intersect(hc_c1,hc_p1),hc_conn);
hc_keeps.connect_only=intersect(hc_keep_conn,hc_conn);
hc_keeps.ds=hc_ds_any;
return
end

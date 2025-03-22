function [ngps,gps,gp_list,nsets_gp,nsets_gp_max]=psg_getgps(sets,ngps_req)
% [ngps,gps,gp_list,nsets_gp,nsets_gp_max]=psg_getgps(sets,ngps_req) is a utility to get partition a list of datasets into groups
%
% sets: cell array, sets{iset}.label has the dataset labe, typically derived from a file name
% ngps_req: number of groups requested (may be omitted)
%
% ngps: number of groups (equal to ngps_req, if supplied
% gps: array of length length(sets), with group assignment
% gp_list: cell array of size ngps, gp_list{k} points to the sets that are in group k
% nsets_gp: array of length ngps, nsets_gp(k)=length(gp_list{k})
% nsets_gp_max: maximum of nsets_gp
%
%   See also: PSG_ALIGN_VARA_DEMO, PSG_ALIGN_VARA_TASK, HLID_GEOM_TRANSFORM_STATS
if (nargin<=1)
    ngps_req=[];
end
nsets=length(sets);
if_ok=0;
while (if_ok==0)
    if isempty(ngps_req)
        ngps=getinp('number of groups','d',[2 nsets]);
    else
        ngps=ngps_req;
    end
    gps=zeros(1,nsets);
    gp_list=cell(1,ngps);
    nsets_gp=zeros(1,ngps); %number of datasets in each group
    sets_avail=[1:nsets];
    for igp=1:ngps
        if ~isempty(sets_avail)
            if igp<ngps
                gp_list{igp}=getinp(sprintf('datasets for group %1.0f',igp),'d',[min(sets_avail) max(sets_avail)],sets_avail);
                gp_list{igp}=intersect(gp_list{igp},sets_avail);
                gps(gp_list{igp})=igp;
                sets_avail=setdiff(sets_avail,gp_list{igp});
            else
                gp_list{igp}=sets_avail;
                gps(sets_avail)=igp;
                sets_avail=[];
            end
        end       
    end
    %are all groups used, and is each stim assigned?
    for igp=1:ngps       
        nsets_gp(igp)=length(gp_list{igp});
        disp(sprintf('group %1.0f',igp))
        for iset=gp_list{igp}
            disp(sprintf(' dataset %2.0f: %s',iset,sets{iset}.label));
        end
    end
    if_ok=1;
    if ~isempty(sets_avail)
        disp('not all datasets assigned');
        if_ok=0;
    end
    if max(gps)<ngps
        disp('not all groups used');
        if_ok=0;
    end
    if (if_ok==1)
        if_ok=getinp('1 if ok','d',[0 1]);
    end
end
nsets_gp_max=max(nsets_gp);
return
end

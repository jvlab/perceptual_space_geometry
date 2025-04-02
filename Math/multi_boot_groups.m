function [boots,gp_info,opts_used]=multi_boot_groups(gps,opts)
% [boots,gp_info,opts_used]=multi_boot_groups(gps,opts) creates lists of bootstrap resamplings within groups
%
% Similar arguments as multi_shuff_groups, but in contrast:
%   * since these are bootstraps, the groups are not mixed, i.e., members
%     of one group can only be used in that group, but each can be used more than once
%   * does not consider an auxiliary tag (this would amount to subgroups)
%   * is not recursive
% 
% gps: a row vector of length n, that indicates the assignment of each
%   element to [1:ngps].  Not all of [1:ngps] must be present.
% opts: options
%  if_log: 1 to log, defaults to 0
%  if_ask: 1 to ask if_exhaust, nboots, defaults to 0
%  if_exhaust: 1 for exhaustive list, defaults to 0
%  if_sortrows: 1 to sort boot list by rows, defaults to 0
%  exhaust_max: maximum number for exhaustive list (defaults to 10^9)
%  nboots: number of bootstraps, defaults to 1000
%  nboots_max: maximum number of boots, defaults to 10^6
%
% boots: the bootstraps; each row is a selection of [1:n], possibly with repetition, and boots(is,:) goes into gps(:)
% gp_info: group information
%   gp_info.ngps: number of groups
%   gp_info.gps_unique: unique group numbers, in order
%   gp_info.gp_list: cell(1,ngps); which elements in each group, in order of gps_unique
%   gp_info.nsets_gp: (1,ngps); number of datasets in each group, in order of gps_unique
% opts_used: options used
% 
%  See also:  MULTI_SHUFF_GROUPS, HLID_GEOM_TRANSFORM_STATS.
%
if nargin<2
    opts=struct;
end
gp_info=struct;
opts=filldefault(opts,'if_log',0);
opts=filldefault(opts,'if_ask',0);
opts=filldefault(opts,'if_exhaust',0);
opts=filldefault(opts,'if_nowarn',0);
opts=filldefault(opts,'if_sortrows',0);
opts=filldefault(opts,'exhaust_max',10^9);
opts=filldefault(opts,'nboots',1000);
opts=filldefault(opts,'nboots_max',10^6);
%
%
n=length(gps);
gps_unique=unique(gps,'stable'); %unique group IDs, but in same order as in gps
ngps=length(gps_unique);
gp_list=cell(1,ngps);
nsets_gp=zeros(1,ngps);
%
boots=zeros(0,n);
%
exhaust=1;
for k=1:ngps
    gp_list{k}=find(gps==gps_unique(k));
    nsets_gp(k)=length(gp_list{k});
    exhaust=exhaust*nsets_gp(k).^nsets_gp(k);
end
gp_info.ngps=ngps;
gp_info.gps_unique=gps_unique;
gp_info.gp_list=gp_list;
gp_info.nsets_gp=nsets_gp;
%
exhmsg='exhaustive list too long; will use random strategy';
%
if opts.if_log | abs(opts.if_ask)
    disp(sprintf('list will have %10.0f bootstraps; max allowed for exhaustive list: %10.0f',exhaust,opts.exhaust_max));
end
%ask whether exhaustive if size permits, otherwise force random (non-exhaustive)
if opts.if_ask
    if exhaust>opts.exhaust_max
        disp(exhmsg);
        opts.if_exhaust=0;
    else
        opts.if_exhaust=getinp('1 for exhaustive bootstraps, 0 for random','d',[0 1]);
    end
    if opts.if_exhaust==0
        opts.nboots=getinp('number of bootstraps','d',[0 opts.nboots_max],opts.nboots);
    end
else %not interactive, make sure that size is OK and if not, log or warn
    if exhaust>opts.exhaust_max & (opts.if_exhaust==1)
        opts.if_exhaust=0;
        if opts.if_log
            disp(exhmsg);
        else
            if opts.if_nowarn==0
                warning(exhmsg);
            end
        end
    end
end
if opts.if_exhaust
    opts.nboots=exhaust;
end
%
%generate the bootstraps, firsst as pointers
%
if opts.if_exhaust
    %first create list of pointers
    s=zeros(0,0);
    for igp=1:ngps
        if nsets_gp(igp)>0
                sadd=multi_boot_all(nsets_gp(igp));
            if isempty(s)
                s=sadd;
            else
                sadd_rep=sadd(ceil([1:size(s,1)*size(sadd,1)]/size(s,1)),:);
                s=[repmat(s,size(sadd,1),1),sadd_rep];
            end
        end        
    end
else %random choices
    s=zeros(opts.nboots,n);
    sofar=0;
    for igp=1:ngps
        if nsets_gp(igp)>0
            s(:,sofar+[1:nsets_gp(igp)])=ceil(rand(opts.nboots,nsets_gp(igp))*nsets_gp(igp));
            sofar=sofar+nsets_gp(igp);
        end
    end
end
%go from s (pointers) to set numbers based on group membership
boots=zeros(size(s));
sofar=0;
for igp=1:ngps
    if nsets_gp(igp)>0
        boots(:,gp_list{igp})=gp_list{igp}(s(:,sofar+[1:nsets_gp(igp)]));
        sofar=sofar+nsets_gp(igp);
    end
end
%
%sort if requested
%
if opts.if_sortrows
    boots=sortrows(boots);
end 
opts_used=opts;
return
end

function s=multi_boot_all(n)
% s=multi_boot_all(n) creates an [n^n n] array of all possible selections
% of n elements from [1:n]
if n<=1
    s=n;
    return
else
    s=[1:n]';
    for k=2:n
        snew=ceil([1:n^k]/(n^(k-1)));
        s=[repmat(s,n,1),snew'];
    end
end
return
end

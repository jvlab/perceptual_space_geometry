function [shuffs,gp_info,opts_used]=multi_shuff_groups(gps,opts)
% [shuffs,gp_info,opts_used]=multi_shuff_groups(gps,opts) creates lists of shuffles
% between groups, that optionally:
%   *is exhaustive
%   *preserve an auxiliary tag
%   *is reduced by symmetry so that groups of the same size, and that contain the same number
%    of items for each tag, are considered identical.  That is, for every set of groups that have the same 
%    number of items with each tag, the groups should be "called out" in gps(shuff(is,:)) from lowest to highest
% 
%%%%%
%to do:
% exhaustive untagged 
% random tagged and untagged, might be able to do together
%%%%
% 
% Implementation of auxiliary tags relies on recursion, but recursion only one step deep
%
% gps: a row vector of length n, that indicates the assignment of each
%   element to [1:ngps].  All of [1:ngps] must be present.
% opts: options
%  if_log: 1 to log, defaults to 0
%  if_ask: 1 to ask if_reduce, if_exhaust, nshuffs, defaults to 0
%  if_reduce: 1 to reduce by symmetrizing, defaults to 0
%  if_exhaust: 1 to do exhaustive list, defaults to 0 (for random list)
%  if_justcount: just count up number of shuffles, do not create them, defaults to 0
%  if_nowarn: 1 to suppress warnings, defaults to 0
%  nshuffs: number of shuffles, defaults to 1000
%  exhaust_raw_max: maximum number for exhaustive list, before reduction (defaults to 10^9)
%  exhaust_reduced_max: maximum number for exhuastive list, after reduction (defaults to 10^6)
%  nshuffs_max: maximum number of shuffles, defaults to min(exhaust_raw_max,exhaust_reduced_max)
%  tags: array of size [1:n], containing values [1:ntags].
%   shuffles will only exchange values with the same tags.
%   If empty or all values equal, this is ignored.
%
% shuffs: the shuffles; each row is a perm of [1:n], and shuffles(is,:) goes into gps(:)
% gp_info: group information
%   gp_info.ngps: number of groups
%   gp_info.gps_unique: unique group numbers, in order
%   gp_info.gp_list: cell(1,ngps); which elements in each group, in order of gps_unique
%   gp_info.nsets_gp: (1,ngps); number of datasets in each group, in order of gps_unique
%   gp_info.exhaust_raw: number of shuffles if exhaustive, prior to reduction
%   gp_info.exhaust_reduced: number of shuffles if exhaustive, if reduced
%   gp_info.tags{itag}: gp_info for subsets with each tag
% a: a list of the shuffles, each row contains 
% opts_used: options used.
% 
%  See also:  NCHOOSEK, FILLDEFAULT, MULTI_SHUFF_ENUM, PSG_ALIGN_VARA_DEMO.
%
if nargin<2
    opts=struct;
end
gp_info=struct;
opts=filldefault(opts,'if_log',0);
opts=filldefault(opts,'if_ask',0);
opts=filldefault(opts,'if_reduce',0);
opts=filldefault(opts,'if_exhaust',0);
opts=filldefault(opts,'if_justcount',0);
opts=filldefault(opts,'if_nowarn',0);
opts=filldefault(opts,'nshuffs',1000);
opts=filldefault(opts,'exhaust_raw_max',10^9);
opts=filldefault(opts,'exhaust_reduced_max',10^6);
opts=filldefault(opts,'nshuffs_max',min(opts.exhaust_raw_max,opts.exhaust_reduced_max));
opts=filldefault(opts,'tags',ones(1,length(gps)));
%
if opts.if_ask
    opts.if_reduce=getinp('1 to reduce shuffles by considering groups of same size (and tag counts) to be equivalent','d',[0 1]);
end
if length(opts.tags)~=length(gps)
    if opts.if_nowarn==0
        warning(sprintf('length(tags), %3.0f, not equal to length(gps), %3.0f; tags ignored',length(opts.tags),length(gps)));
    end
    opts.tags=ones(1,length(gps));
end
%
if_tagged=double(~(min(opts.tags)==max(opts.tags)));
%
n=length(gps);
gps_unique=unique(gps);
ngps=length(gps_unique);
gp_list=cell(1,ngps);
nsets_gp=zeros(1,ngps);
%
shuffs=zeros(0,n);
%
for k=1:ngps
    gp_list{k}=find(gps==gps_unique(k));
    nsets_gp(k)=length(gp_list{k});
end
exhaust_raw=1;
exhaust_reduced=1;
%
%if if_tagged, recursive call to get counts
%
if if_tagged
    opts_recur=opts;
    opts_recur.if_log=0;
    opts_recur.if_ask=0;
    opts_recur.if_reduce=0; %only reduce at end
    opts_recur.if_exhaust=0; %irrelevant since if_justcount=1
    opts_recur.if_justcount=1;
    opts_recur.if_nowarn=1;
    tags=unique(opts.tags);
    ntags=length(tags);
    gp_info.tags=cell(1,ntags);
    gps_tag=cell(1,ntags);
    for itag=1:ntags %treat each tag as a separate subset
        gps_tag{itag}=gps(opts.tags==tags(itag));
        if ~isempty(gps_tag{itag})
            if opts.if_log
                disp(sprintf('multi_shuff_groups called recursively for counting, tag %1.0f [applies to %4.0f items]',tags(itag),length(gps_tag{itag})));
            end
            opts_recur.tags=ones(1,length(gps_tag{itag})); %will prevent further recursion
            [shuffs,gp_info_tag]=multi_shuff_groups(gps_tag{itag},opts_recur);
            exhaust_raw=exhaust_raw*gp_info_tag.exhaust_raw;
            gp_info.tags{itag}.gp_info=gp_info_tag;
            gp_info.tags{itag}.tag=tags(itag);
        end
    end
    exhaust_reduced=exhaust_raw;
    %determine whather any groups have the same number of each tag
    gp_profile=zeros(ngps,ntags);
    for igp=1:ngps
        for itag=1:ntags
            gp_profile(igp,itag)=sum((gps==gps_unique(igp)) & opts.tags==tags(itag));
        end
    end
    gp_profile_unique=unique(gp_profile,'rows');
    nu=size(gp_profile_unique,1);
    for iu=1:nu
        nk=sum(all(repmat(gp_profile_unique(iu,:),size(gp_profile,1),1)==gp_profile,2));
        exhaust_reduced=exhaust_reduced/factorial(nk);      
    end
else
    for k=2:ngps
        if k>1 %work iteratively to avoid quotients of large integers
            newmult=nchoosek(sum(nsets_gp(1:k)),nsets_gp(k));
            exhaust_raw=exhaust_raw*newmult;
            exhaust_reduced=exhaust_reduced*newmult/sum(double(nsets_gp(1:k)==nsets_gp(k))); %divide by multiplicity
        end
    end
    if sum(nsets_gp)~=n
        if opts.nowarn==0
            warning('some group numbers are not integers in the range [1 ngps]');
        end
    end
end
gp_info.ngps=ngps;
gp_info.gps_unique=gps_unique;
gp_info.gp_list=gp_list;
gp_info.nsets_gp=nsets_gp;
gp_info.exhaust_raw=exhaust_raw;
gp_info.exhaust_reduced=exhaust_reduced;
gp_info.if_tagged=if_tagged;
%
exhmsg='exhaustive list too long; will use random strategy';
if ~opts.if_justcount
    %choose reduced vs. non-reduced
    if opts.if_reduce
        exhaust=exhaust_reduced;
        reduced_string='reduced';
        exhaust_max=opts.exhaust_reduced_max;
    else
        exhaust=exhaust_raw;
        reduced_string='not reduced';
        exhaust_max=opts.exhaust_raw_max;
    end
    %show info about number of shuffles
    if opts.if_log | opts.if_ask
        disp(sprintf('list (%s) will have %10.0f shuffles; max allowed for exhaustive list: %10.0f',reduced_string,exhaust,exhaust_max));
    end
    %ask whether exhaustive if size permits, otherwise force random (non-exhaustive)
    if opts.if_ask
        if exhaust>exhaust_max
            disp(exhmsg);
            opts.if_exhaust=0;
        else
            opts.if_exhaust=getinp('1 for exhaustive shuffles, 0 for random','d',[0 1]);
        end
        if opts.if_exhaust==0
            opts.nshuffs=getinp('number of shuffles','d',[0 opts.nshuffs_max],opts.nshuffs);
        end
    else %not interactive, make sure that size is OK and if not, log or warn
        if exhaust>exhaust_max & (opts.if_exhaust==1)
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
        opts.nshuffs=exhaust;
    end
    nshuffs=opts.nshuffs;
    if opts.if_exhaust
        opts_enum=struct;
        opts_enum.if_log=opts.if_log;
        %
        %enumerate all
        %
        if if_tagged %tags present: actor by restriction  and at end apply reduction by symmetry
            opts_enum.if_reduce=0; %reduce at end           
            for itag=1:ntags %treat each tag as a separate subset
                % gps_tag{itag}: the group labels for the subset with this tag
                gp_list_tag=cell(1,ngps);
                for k=1:ngps
                    gp_list_tag{k}=intersect(find(gps==gps_unique(k)),find(opts.tags==tags(itag)));
                end
                items_tag=find(opts.tags==tags(itag));
                disp(sprintf('multi_shuff_enum called for tag %1.0f [applies to %4.0f items]',tags(itag),length(gps_tag{itag})));
                [a,opts_enum_used]=multi_shuff_enum(gp_profile(:,itag),opts_enum);
                %combine this shuffle with previous as a "tensor product"
                 if (itag>1)
                    s0=size(shuffs,1);
                    shuffs=repmat(shuffs,size(a,1),1); %a copy of each of previous shuffles, one for each row of a
                else
                    shuffs=zeros(size(a,1),n);
                    s0=1;
                end
                %convert a into shuffles, using only with the items from each group that have tags(itag) 
                %***this can be optimized*** -- could do a range of values of is_big at once 
                for is_big=1:size(shuffs,1) %fast index
                    for igp=1:ngps
                        is=ceil(is_big/s0); %slow index
                        shuffs(is_big,gp_list_tag{igp})=items_tag(find(a(is,:)==igp));
                    end
                end %is
             end
            %apply opts.if_reduce
            if (opts.if_reduce)
                if (opts.if_log)
                    disp('reducing by symmetry');
                    disp('unique group profiles:')
                    disp(gp_profile_unique);
                    disp('group profiles');
                    disp(gp_profile);
                end
                if (size(gp_profile,1)~=size(gp_profile_unique,1))
                    gp_table=gps(shuffs); %gp_table has the group assignments
                    for iu=1:nu
                        gp_red=find(all(repmat(gp_profile_unique(iu,:),size(gp_profile,1),1)==gp_profile,2));
                        if (opts.if_log)
                            disp(sprintf(' unique group profile %1.0f is %s',iu,sprintf(' %2.0f',gp_profile_unique(iu,:))))
                            disp(sprintf('occurs for group labels %s',sprintf(' %2.0f',gps_unique(gp_red))));
                        end
                        if length(gp_red)>1
                            %do the reduction:
                            %select the rows for which the first appearance of each group in gps_unique(gp_red) occurs in the
                            %order of gps_unique(gps_red)
                            first_igs=zeros(size(gp_table,1),length(gp_red));
                            for ig=1:length(gp_red)
                                first_igs(:,ig)=max(double(gp_table==gps_unique(gp_red(ig))).*repmat(fliplr([1:n]),size(gp_table,1),1),[],2); %first_ig is lage when first occurrence of gp_red(ig) comes early
                            end
                            seq_ok=all(diff(first_igs,1,2)<0,2);
                            shuffs=shuffs(seq_ok,:);
                            gp_table=gp_table(seq_ok,:);
                        end
                    end
                end %are there equivalent groups?
            end
        else
            %untagged:  can apply reduction by symmetry within mult_shuff_enum
            %convert a into shuffles
            % shuffs=zeros(size(a));
            % for is=1:size(shuffs,1)
            %     for igp=1:ngps
            %         shuffs(is,gp_list{igp})=find((a(is,:)==igp));
            %     end
            % end

        end
    else 
    %   
    %random shuffles
    %
    % 
    % If random — factor by restriction subsets, combine by tensor product, apply reduction by symmetry at end
    end
end %just count
opts_used=opts;
return


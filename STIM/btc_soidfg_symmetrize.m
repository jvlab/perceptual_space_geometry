function [dsym,dall,pairings,unpaired,opts_used]=btc_soidfg_symmetrize(d,opts)
% [dsym,dall,pairings,unpaired,opts_used]=btc_soidfg_symmetrize(d,transform_type) symmetrizes a figure-ground database around the origin
%
% This corresponds to the function performed by btc_soid_sym, if symopts.if_symm=1, which uses a linear transformation (arithmetic mean)
% input:
%  d: a database (cell array) of psychophyical data for figure-ground experiments
%    Typical entry (d{k}): 
%    subj_id: 'MC'
%      files: {1×32 cell}
%      ndirs: 8
%       dirs: [1×1 struct]
%     weib_a: [8×3 double]
%     weib_b: [8×3 double]
%   opts: options (all defaults if not supplied)
%     transform_type; how to average ('lin','log','recip','ask'), corresponding to
%       arithmetic mean, geometric mean, or harmonic mean, or ask at keyboard, defaults to 'lin' 
%     tol: tolerance for matching symmetric contrasts, or for accepting weibull b's as equal, defaults to 10^-5
%     weib_b_forceavg: defaults to 0 (combines weib_b parameters by averaging only if they differ, otherwise copy);
%       1->force combining by average
%     verbose: 1->summarize operation to console (default)
%     sym_type: 'sign' (only option at present)
%     pair_select: how to select one from each pair
%       'gpfn': first coord ground positive, if ground is zero, then figure negative -- default -- results in plots in upper half plane
%       'fpgn': first coord figure positive, if figure zero, then if ground is negative -- results in plots in right half half plane
%       'orig': order of original database, might result in different points chosen in each plane
%
% output:
%  dsym:  corresponding structure with thresholds (weib_a) and shape params (weib_b) symmetrized around the origin
%    each segment has a field d{k}.symmetrized='one of each pair'
%    only one of each symmetric pairing are included
%    ndirs and dirs reflects the number of pairs
%    weib_a and weib_b: values combined by averaging, lo and high error bars combined by (1/sqrt(2) RMS deviation about the mean
%  dall:  corresopnding structure with all directions included
%    each segment has a field d{k}.symmetrized='pairs and unpaired'
%    directions with a symmetric pairing are included along with their symmetric pair (so the thresholds are identical)
%    directions without a symmetric pairing are included as is
%  pairings: cell array indicating the paired index of each direction, NaN if unpaired or more than one match
%  unpaired:  array of size length(d) indicating the number unpaired for each segment
%  opts_used: options used
%
% See also:  FIGGND_DBASE_MAKE, BTC_SOIDFG_DEMO, FIGGND_PROC, FILLDEFAULT, BTC_SOID_SYM, MTC_SOID_AVG.
%
transform_type_avail={'lin','log','recip'};
transform_type_avail_string={'arithmetic','geometric','harmonic'};
if (nargin<=1)
    opts=[];
end
opts=filldefault(opts,'transform_type','lin');
opts=filldefault(opts,'tol',10^-5);
opts=filldefault(opts,'weib_b_forceavg',0);
opts=filldefault(opts,'verbose',1);
opts=filldefault(opts,'sym_type','sign');
opts=filldefault(opts,'pair_select','gpfn');
%
if strcmp(opts.transform_type,'ask')
    for itransform=1:length(transform_type_avail)
        disp(sprintf('%1.0f->%5s transform, i.e., %s mean',itransform,transform_type_avail{itransform},transform_type_avail_string{itransform}));
    end
    itransform_choice=getinp('choice','d',[1 length(transform_type_avail)],1);
    opts.transform_type=transform_type_avail{itransform_choice};
end
transform_ptr=strmatch(opts.transform_type,transform_type_avail,'exact');
opts.transform_type_avail=transform_type_avail;
opts.transform_type_avail_string=transform_type_avail_string;
opts_used=opts;
%    
nsegs=length(d);
dsym=cell(1,nsegs);
dall=cell(1,nsegs);
pairings=cell(1,nsegs);
unpaired=zeros(1,nsegs);
ndirs_tot=0;
allpaired_tot=0;
oneofpair_tot=0;
unpaired_tot=0;
%symmetrize each segment
for iseg=1:nsegs
    ndirs=d{iseg}.ndirs;
    ndirs_tot=ndirs_tot+ndirs;
    %collect coord values, across figure and ground
    cvals=zeros(0,ndirs);
    dirs=d{iseg}.dirs;
    dir_fields=fieldnames(dirs); %typically {'fig','gnd'}
    cvals_fg=cell(1,length(dir_fields));
    for ifg=1:length(dir_fields)
        fg=dirs.(dir_fields{ifg});
        cvals_fg{ifg}=zeros(0,ndirs);
        btc_fields=fieldnames(fg);
        for ibtc=1:length(btc_fields)
            cvals_fg{ifg}=[cvals_fg{ifg};fg.(btc_fields{ibtc})];
        end
        cvals=[cvals;cvals_fg{ifg}];
    end
    %find pairings
    p=NaN(1,ndirs);
    for idir=1:ndirs
        matches=find(all(abs(repmat(cvals(:,idir),1,ndirs)+cvals)<=opts.tol),1);
        if length(matches)==1
            p(idir)=matches;
        end
    end
    pairings{iseg}=p;
    switch opts.pair_select
        case 'orig'
            indices_oneofpair=find([1:ndirs]<p);
        case 'gpfn'
            pok=find(~isnan(p));
            gp=pok(find(cvals_fg{2}(1,pok)>0)); %pointers into p where ground is positive
            gz=pok(find(cvals_fg{2}(1,pok)==0)); %pointers into p where ground is zero
            fn=pok(find(cvals_fg{1}(1,pok)<0)); %pointers into p where figure is negative
            indices_oneofpair=[gp intersect(gz,fn)];
        case 'fpgn'
            pok=find(~isnan(p));
            fp=pok(find(cvals_fg{1}(1,pok)>0)); %pointers into p where figure is positive
            fz=pok(find(cvals_fg{1}(1,pok)==0)); %pointers into p where figure is zero
            gn=pok(find(cvals_fg{2}(1,pok)<0)); %pointers into p where ground is negative
            indices_oneofpair=[fp intersect(fz,gn)]
    end
    indices_allpaired=[indices_oneofpair p(indices_oneofpair)]; %first of each pair, then second of each pair
    unpaired(iseg)=sum(isnan(p));
    indices_unpaired=find(isnan(p));
    allpaired_tot=allpaired_tot+length(indices_allpaired);
    oneofpair_tot=oneofpair_tot+length(indices_oneofpair);
    unpaired_tot=unpaired_tot+length(indices_unpaired);
    %verify that all are either paired or unpaired
    indices_check=[indices_allpaired indices_unpaired];
    if length(unique(indices_check))~=ndirs
        p
        indices_oneofpair
        indices_allpaired
        indices_unpaired
        error(sprintf('pairing error, segment %3.0f',iseg));
    end
    %set up metadata
    d_nodata=rmfield(rmfield(d{iseg},'weib_a'),'weib_b');
    dall{iseg}=d_nodata;
    dall{iseg}.symmetrized=cat(2,opts.sym_type,', pairs and unpaired');    
    dsym{iseg}=d_nodata;
    dsym{iseg}.symmetrized=cat(2,opts.sym_type,', one of each pair');
    dsym{iseg}.ndirs=length(indices_oneofpair);
    for ifg=1:length(dir_fields)
        fg=dirs.(dir_fields{ifg});
        btc_fields=fieldnames(fg);
        for ibtc=1:length(btc_fields)
            c=fg.(btc_fields{ibtc});
            dall{iseg}.dirs.(dir_fields{ifg}).(btc_fields{ibtc})=c([indices_allpaired indices_unpaired]);
            dsym{iseg}.dirs.(dir_fields{ifg}).(btc_fields{ibtc})=c(indices_oneofpair);
        end
    end
    %weibull_a
    %for dsym:  just one from each pair
    for iv=1:length(indices_oneofpair)
        weib_a=[d{iseg}.weib_a(indices_oneofpair(iv),:);d{iseg}.weib_a(p(indices_oneofpair(iv)),:)]; %stack Weibull data for member of each pair
        w_avg=mtc_soid_avg(weib_a);
        dsym{iseg}.weib_a(iv,:)=w_avg(:,transform_ptr,1,2)'; %d3=1 for unweighted, d4=2 for NaN if any are NaN
    end
    %for dall: both from each pair and the unpaired
    for iv=1:length(indices_allpaired)
        weib_a=[d{iseg}.weib_a(indices_allpaired(iv),:);d{iseg}.weib_a(p(indices_allpaired(iv)),:)]; %stack Weibull data for member of each pair
        w_avg=mtc_soid_avg(weib_a);
        dall{iseg}.weib_a(iv,:)=w_avg(:,transform_ptr,1,2)'; %d3=1 for unweighted, d4=2 for NaN if any are NaN
    end
    dall{iseg}.weib_a(length(indices_allpaired)+[1:length(indices_unpaired)],:)=d{iseg}.weib_a(indices_unpaired,:);
    %
    %weibull_b
    %
    %to average or copy?
    if opts.weib_b_forceavg==0
        weib_b_range=max(abs(max(d{iseg}.weib_b,[],1)-min(d{iseg}.weib_b,[],1)));
        weib_b_avg=double(weib_b_range>opts.tol);
    else
        weib_b_avg=1;
    end
    if (weib_b_avg)
        %for dsym:  just one from each pair
        for iv=1:length(indices_oneofpair)
            weib_b=[d{iseg}.weib_b(indices_oneofpair(iv),:);d{iseg}.weib_b(p(indices_oneofpair(iv)),:)]; %stack Weibull data for member of each pair
            w_avg=mtc_soid_avg(weib_b);
            dsym{iseg}.weib_b(iv,:)=w_avg(:,transform_ptr,1,2)'; %d3=1 for unweighted, d4=2 for NaN if any are NaN
        end
        %for dall: both from each pair and the unpaired
        for iv=1:length(indices_allpaired)
            weib_b=[d{iseg}.weib_b(indices_allpaired(iv),:);d{iseg}.weib_b(p(indices_allpaired(iv)),:)]; %stack Weibull data for member of each pair
            w_avg=mtc_soid_avg(weib_b);
            dall{iseg}.weib_b(iv,:)=w_avg(:,transform_ptr,1,2)'; %d3=1 for unweighted, d4=2 for NaN if any are NaN
        end
        dall{iseg}.weib_b(length(indices_allpaired)+[1:length(indices_unpaired)],:)=d{iseg}.weib_b(indices_unpaired,:);
    else %copy
        dall{iseg}.weib_b=d{iseg}.weib_b;
        dsym{iseg}.weib_b=d{iseg}.weib_b(indices_oneofpair,:);
    end
end
if (opts.verbose)
    disp(sprintf('database symmetrized by %s with %4s transform. %3.0f segments, %4.0f total directions->%4.0f paired (%4.0f pairs) and %4.0f unpaired',...
        opts.sym_type,opts.transform_type,nsegs,ndirs_tot,allpaired_tot,oneofpair_tot,unpaired_tot));
end
return
end

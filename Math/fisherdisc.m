function [results,opts_used]=fisherdisc(data,tags,opts)
% [results,opts_used]=fisherdisc(data,tags,opts)
% computes the Fisher discriminant and associated statistics, 2-class, with cross-validation
% (no assumption of equal within-class covariances) 
%
%   see FisherDiscrimWikipedia.pdf
%
%  data: size(data)=[nfeats,nsamps], i.e., one sample per column, one row per feature
%  tags: size(tags)=[1 nsamps], containing 1 or 2 to tag each sample
%      Note  20Aug 2016: Can also have a tag of 0.  These items are then classified according to 
%      classifiers built on the nonzero tags.  However, the jackknife or shuffle
%      procedures will include these data entries in the shuffle, so these outputs 
%      should be ignored. Inclusion of data points with tag=0 affects results fields
%      ss_total but not ss_eachclass or ss_class or class_mean or class_cov.
%  opts: options
%    opts.nshuffle: number of shuffles to apply.
%       Defaults to 0.  Can also be Inf -- exhaustive, or -1'
%       if -1, then uses exhausive shuffling, or random selection of opts.nshuffle_max
%       if exhaustive would exceed opts.nshuffle_max.  If number requested is greater
%       than number necessary for exhaustive, then exhaustive is used.
%    opts.nshuffle_max: maximum number of shuffles if opts.nshuffle=-1; defaults to 200
%    opts.xvalid:  1 to do drop-one cross-validation and jackknife stats, defaults to 0
%       If this is 1, then opts.classifiers must be non-empty
%    opts.sub1:   1 to subtract 1 for empirical estimate of within-class covariances. This defaults to 0 for backward compatibility,
%       but a value of 1 is appropriate unless the sampling of the class is complete
%        Note that this only matters if the number of elements in each class is unequal
%    opts.classifiers: cell array of names of classifiers to use, as cell array
%      'halfway' (halfway between the two midpoints)
%      'mapequal' (maximum a posteriori based on equal class occupancies)
%      'mapbayes' (maximum a posteriori based on empiric class occupancies)
%      Cannot be empty if opts.xvalid=1
%    opts.condmax: maximum condition number for covariance matrix, defaults to Inf
%    opts.segs: defaults to [], otherwise a cell array of subsets of trials
%      that are considered to be non-independent because they constitute a
%      contiguous segment
%    opts.delrad:  separation of trials within a segment that are considered dependent.
%      Defaults to Inf, so that everything within a segment is dependent.
%      A value of 1 means that a trial is considered as depedent on its
%      neighbor.  A value of 2 means that a trial is considered as
%      dependent on its neighbor and its next-nearest neighbor.  A value of
%      0 removes the effect of segmentation, since all trials are
%      considered independent even if within the same segment
%    opts.nflips:  number of surrogate datasets to be created by flipping
%      tags within a segment, via findsegflips
%    opts.nflips_maxtries: maximum number of tries to find flips
%    opts.nflips_tol:  tolerance on change in number of trials assigned to each group
%      
%   results: results structure
%      results.discriminant [1 nfeat] a unit-vector discriminant, positive values favor class 2
%      results.projections [1 nsamps] values of the projections of each sample
%      results.ss_eachclass: [1 2] insample sum of squares, within each class
%      results.ss_class: insample sum of squares, total within each class
%      results.ss_total: insample sum of squares, total
%      results.count_eachclass: count of insample values in each class
%      results.var: variance within each class, calculated separately for each class
%      results.var_equated: variance within class, assuming within-class variances are the same in the two classes
%      results.fc_* (not cross-validated)
%      results.fc_[name]_classes: [1 nsamps] 1 or 2 based on classification
%      results.fc_[name]_cmatrix: [2 2] confusion matrix, row=correct class, column=assigned class
%      results.fc_[name]_ncorr: number correct
%      results.fc_[name]_p: fraction of times that a shuffled dataset yields an equal or greater number correct
%      results.fc_[name]_cutpoint, projections greater than cutpoint are assigned to class 2
%            (cutpoint not calculated for mapequal or mapbayes, since the discrimination may be supported by an interval)
%      results.fc_[name]_logprior: log likelihoods of prior (0 for mapequal and halfway, logs of class counts for mapbayes)
%      results.fc_[name]_logposterior: log likelihoods of posterior -- for halfway, uses var_equated, for mapequal or mapbayes, uses in-class estimates of variances
%      results.fc_[name]_logposterior_diff: alternate calculation of log posterior for halfway
%      ***the next fields are only calculated if opts.xvalid=1
%      results.xv_[name]_classes: [1 nsamps] 1 or 2 based on classification
%      results.xv_[name]_cmatrix: [2 2] confusion matrix, row=correct class, column=assigned class
%      results.xv_[name]_ncorr: number correct
%      results.discriminant_jdebiased: [1 nfeat] jackknife-debiased discriminant (not orthonormalized)
%      results.discriminant_jsem: [1 nfeat] std devs from jackknife
%          (_jdebiased is NOT orthonormalized, and wa not used to calculate classes; classes
%          are calculated based on discriminant, or on individual drop-one discriminants)
%    ***the next fields are only calculated if opts.segs is non-empty
%      results.sg_whichseg: [1 nsamps] whichseg(isamp) is the segment that contains a sample; 0 means unassigned
%      results.sg_whichcor: {nsamps} sg_whichcor{isamp} is a list of the samples that are correlated with sample isamp
%    ***the next fields are only calculated if opts.segs is non-empty AND opts.xvalid=1
%      results.sg_[name]_classes: [1 nsamps] 1 or 2 based on classification
%      results.sg_[name]_cmatrix: [2 2] confusion matrix, row=correct class, column=assigned class
%    ***the next fields are only calculated if opts.segs is non-empty AND opts.nflips>0
%      results.segflip.*: ss, fc for segflips
%      results.segflip.segflips: [nflips nsegs] 1 if a segment is to be flipped on a given surrogate
%      results.segflip.aux: auxiliary output from findsegflips
%      results.segflip.msg: message from findsegflips
%      results.segflip.taglist: [nflips ntrials] the tags from the surrogates made by flipping
%    ***the next fields are only calculated if opts.nshuffle is nonzero
%      results.shuffle.*: ss, fc for shuffle
%      results.shuffle.msg: details about shuffle calculation
%    IMPORTANT
%         When the covariances are highly non-circular and signal is small, the discriminant is
%         shaped by the covariances, rather than the signal -- see scenario 2
%         of fisherdisc_test -- and in this case, the error bars may be small,
%         but the out-of-sample testing shows that the discriminator is not valid
%
%   opts_used: options used
%
% 23Dec20: Clarify documentation re sub1, count_eachclass, var_equated
% 30May25: Warn if opts.xvalid=1 and classifier not specified.  Add calculation of number correct and p-value
%
% See also:
%  FISHERDISC_TEST, FISHERDISC_DEF, FINDSEGFLIPS, JACKSA, FISHERDISC_CLASSIFY.
%
if (nargin<=2)
    opts=[];
end
opts=fisherdisc_def(opts);
if opts.xvalid>0 & isempty(opts.classifiers)
    warning('cross-validation requested (opts.xvalid=1) but classifer not specified.  opts.xvalid set to 0.');
    opts.xvalid=0;
end
results=[];
results.errormsg=[];
%
results_all=fisherdisc_do(data,tags,opts);
if ~isempty(results_all.errormsg)
    results.errormsg=strvcat(results.errormsg,results_all.errormsg);
end
nfeat=size(data,1);
nsamps=size(data,2);
for ig=1:2
    ptrs{ig}=find(tags==ig);
    n(ig)=length(ptrs{ig});
end
fn=fieldnames(results_all);
%keep all fields that end in 'in' but strip '_'in'; drop the '_out' fields, keep all others
for ifn=1:length(fn)
    fname=fn{ifn};
    if (strcmp(fname(end-2:end),'_in'))
        results=setfield(results,fname(1:end-3),getfield(results_all,fname));
    elseif (strcmp(fname(end-3:end),'_out'))
    else
        results=setfield(results,fname,getfield(results_all,fname));
    end
    %added 20 Aug 16
    %replace insample and outsample as returned from fisherdisc_do in results_all; 
    %those had insample=[1:length(tags) and outsample=[];
    results.insample=find(tags>0);
    results.outsample=find(tags==0);
end
clist=opts.classifiers;
%
%calculate in-sample fraction correct 
%
for ic=1:length(clist)
    results.(cat(2,'fc_',opts.classifiers{ic},'_ncorr'))=sum(diag(results.(cat(2,'fc_',opts.classifiers{ic},'_cmatrix'))));
end
%
% parse the segments, if supplied
%
if ~isempty(opts.segs)
   % make sure that all sgements are disjoint
   % add singleton segments and figure out pointers from dropped values to each segment
   results.sg_whichseg=zeros(1,nsamps);
   results.sg_whichcor=cell(1,nsamps);
   for iseg=1:length(opts.segs)
        itrials=opts.segs{iseg};
        if (any(results.sg_whichseg(itrials)>0))
            results.errormsg=strvcat(results.errormsg,'all: segment definitions overlap; overlapping segment definition ignored.');
        else
            results.sg_whichseg(itrials)=iseg;
        end
   end
   for itrial=1:nsamps
        if (results.sg_whichseg(itrial)==0) %a singleton
            results.sg_whichseg(itrial)=1+max(results.sg_whichseg); %give it a unique segment number
            results.sg_whichcor{itrial}=itrial;
        else
            trials=opts.segs{results.sg_whichseg(itrial)};
            ptr=find(itrial==trials);
            indices=[max(1,ptr-opts.delrad):min(length(trials),ptr+opts.delrad)]; %use delrad as a range
            results.sg_whichcor{itrial}=trials(indices);
        end
   end
   %create a list of segment flips
   if opts.nflips>0
       [segflips,aux,fsf_opts_used]=findsegflips(tags,results.sg_whichseg,opts);
       segflip_details.segflips=segflips;
       segflip_details.aux=aux;
        if ~isempty(aux.errormsg)
            results.errormsg=strvcat(results.errormsg,sprintf('findsegflips: %s',aux.errormsg));
        end
        nflips=size(segflip_details.segflips,1);
        segflip_details.taglist=repmat(tags,nflips,1);
        %form the tags for each flip surrogate
        for iflip=1:nflips 
            for iseg=find(segflip_details.segflips(iflip,:)==1)
                trials=find(results.sg_whichseg==iseg);
                segflip_details.taglist(iflip,trials)=3-segflip_details.taglist(iflip,trials);
            end
        end
   end
end %opts.segs
%
% if cross-validation is requested: do the drop-one analyses
% for classification, and also for jackknifing
%
if (opts.xvalid) %cross-validation
    jnaive.discriminant=results.discriminant;
    jdrop=[];
    xv=[];
    sg=[];
    for ic=1:length(clist)
        cname=clist{ic};
        sname=cat(2,'fc_',cname);
        name_fc{ic}=cat(2,'fc_',cname,'_classes');
        name_xv_classes{ic}=cat(2,'xv_',cname,'_classes');
        name_xv_cmatrix{ic}=cat(2,'xv_',cname,'_cmatrix');
        name_xv_ncorr{ic}=cat(2,'xv_',cname,'_ncorr');
        xv=setfield(xv,name_xv_classes{ic},zeros(1,nsamps));
        name_sg_classes{ic}=cat(2,'sg_',cname,'_classes');
        name_sg_cmatrix{ic}=cat(2,'sg_',cname,'_cmatrix');
        sg=setfield(sg,name_sg_classes{ic},zeros(1,nsamps));
   end
   %do the Fisher discriminant with each sample flagged as out-of-sample
   for idrop=1:nsamps
        outflag=zeros(1,nsamps);
        outflag(idrop)=1;
        results_drop=fisherdisc_do(data,tags,opts,outflag);
        if ~isempty(results_drop.errormsg)
            results.errormsg=strvcat(results.errormsg,sprintf('drop: %4.0f: %s',idrop,results_drop.errormsg));
        end
        jdrop(idrop).discriminant=results_drop.discriminant;
        for ic=1:length(clist)
            xv=setfield(xv,name_xv_classes{ic},{idrop},getfield(results_drop,name_fc{ic},{idrop}));
        end
   end %idrop
   %confusion matrices for cross-validated classifiers
   for ic=1:length(clist)
       xvc=getfield(xv,name_xv_classes{ic});
       for ig_row=1:2
           for ig_col=1:2
                cmat(ig_row,ig_col)=sum(xvc(ptrs{ig_row})==ig_col);
            end
       end
       xv=setfield(xv,name_xv_cmatrix{ic},cmat);
       xv=setfield(xv,name_xv_ncorr{ic},sum(diag(cmat)));
   end
   xvfn=fieldnames(xv);
   %install the fields
   for ifn=1:length(xvfn)
       results=setfield(results,xvfn{ifn},getfield(xv,xvfn{ifn}));
   end
   %
   % do jackknife
   %
   [jbias,jdebiased,jvar,jsem]=jacksa(jnaive,jdrop);
   results.discriminant_jdebiased=jdebiased.discriminant; %not orthonormal
   results.discriminant_jsem=jsem.discriminant; %not orthonormal
   %
   % do segmented cross-validation
   %
   if ~isempty(opts.segs)
       %
       %When dropping a trial, mark everything within that segment out to opts.delrad as out-of-sample
       %
       %This next loop drops the trials that are correlated with a given trial, and calculates
       %the cross-validated classification.  Note that it is inefficient if opts.delrad=Inf, since it fails
       %to make use of the fact that several classifications can be computed by dropping a complete segment
       for itrial=1:nsamps
            outflag=zeros(1,nsamps);
            outflag(results.sg_whichcor{itrial})=1; %mark a portion of segment as out-of-sample
            results_drop=fisherdisc_do(data,tags,opts,outflag);
            if ~isempty(results_drop.errormsg)
                results.errormsg=strvcat(results.errormsg,sprintf('segd: %4.0f: %s',idrop,results_drop.errormsg));
            end
            for ic=1:length(clist)
                sg=setfield(sg,name_sg_classes{ic},{itrial},getfield(results_drop,name_fc{ic},{itrial}));
            end
        end %itrial
       %confusion matrices for cross-validated classifiers
       for ic=1:length(clist)
           sgc=getfield(sg,name_sg_classes{ic});
           for ig_row=1:2
                for ig_col=1:2
                    cmat(ig_row,ig_col)=sum(sgc(ptrs{ig_row})==ig_col);
                end
           end
           sg=setfield(sg,name_sg_cmatrix{ic},cmat);
       end
       sgfn=fieldnames(sg);
       %install the fields
       for ifn=1:length(sgfn)
            results=setfield(results,sgfn{ifn},getfield(sg,sgfn{ifn}));
       end
    end %isempty(opts.segs)
end %if xvalid
%
% if shuffle is requested, determine how good the classification is
% with real data and with shuffle.  Cross-validation not done on shuffled sets.
%
if ~(opts.nshuffle==0)
%      first, determine whether to do an exhausive shuffle
%      if 'uptomax' then uses exhausive shuffling, or random selection 
%      if 'all' then do all
%
    nexhaust_log=gammaln(nsamps+1)-gammaln(n(1))-gammaln(n(2));
    if (exp(nexhaust_log) > 10^14)
        nexhaust=exp(nexhaust_log);
    else
        nexhaust=nchoosek(nsamps,n(1));
    end
    %
    if opts.nshuffle>0 & opts.nshuffle<Inf
        if opts.nshuffle>nexhaust
            ns=nexhaust;
            ifexhaust=1;
        else
            ns=opts.nshuffle;
            ifexhaust=0;
        end
    elseif opts.nshuffle==Inf
        ns=nexhaust;
        ifexhaust=1;
    else
        if (nexhaust<=opts.nshuffle_max)
            ns=nexhaust;
            ifexhaust=1;
        else
            ns=opts.nshuffle_max;
            ifexhaust=0;
        end
    end
    shuffle_msg=sprintf(' nsamps=%4.0f n=[%4.0f %4.0f] opts.nshuffle=%4.0f nexhaust=%17.10g ->ns=%10.0f ifexhaust=%1.0f',...
        nsamps,n,opts.nshuffle,nexhaust,ns,ifexhaust');
    sorted=[ptrs{1} ptrs{2}];
    if (ifexhaust)
        %establish the list of tags, without assuming that tags were originally in order
        shuffs=nchoosek([1:nsamps],n(1));
    end
    a=[]; %accumulate for shuffle
    ncorr=zeros(ns,length(clist));
    for ishuf=1:ns
        if (ifexhaust)
            tperm=[shuffs(ishuf,:) setdiff([1:nsamps],shuffs(ishuf,:))];
        else
            tperm=randperm(nsamps);
        end
        tags_shuffled=tags(tperm);
        %
        r=fisherdisc_do(data,tags_shuffled,opts);
        if ~isempty(r.errormsg)
            results.errormsg=strvcat(results.errormsg,sprintf('shuf: %4.0f: %s',ishuf,r.errormsg));
        end
        a.ss_total(1,ishuf)=r.ss_total_in;
        a.ss_eachclass(:,ishuf)=r.ss_eachclass_in;
        a.ss_class(:,ishuf)=r.ss_class_in;
        for ic=1:length(clist)
            shuf_cmatrix(:,:,ishuf,ic)=getfield(r,cat(2,'fc_',clist{ic},'_cmatrix_in'));
            ncorr(ishuf,ic)=sum(diag(shuf_cmatrix(:,:,ishuf,ic))); %added 30May25
        end
    end
    for ic=1:length(clist)
        a=setfield(a,cat(2,'fc_',clist{ic},'_cmatrix'),shuf_cmatrix(:,:,:,ic));
        a=setfield(a,cat(2,'fc_',clist{ic},'_ncorr'),ncorr(:,ic)'); %added 30May25
        results.(cat(2,'fc_',clist{ic},'_p'))=sum((ncorr(:,ic)>=results.(cat(2,'fc_',clist{ic},'_ncorr'))))/ns;  %added 30May25
    end
    results.shuffle=a;
    results.shuffle.msg=shuffle_msg;
end
%
% if flips based on segmentation is requested
%
if ~(opts.nflips==0) & ~isempty(opts.segs)
    a=[]; %accumulate for flips
    nflips=size(segflip_details.segflips,1);
    for iflip=1:nflips
        r=fisherdisc_do(data,segflip_details.taglist(iflip,:),opts);
        if ~isempty(r.errormsg)
            results.errormsg=strvcat(results.errormsg,sprintf('segflip: %4.0f: %s',iflip,r.errormsg));
        end
        a.ss_total(1,iflip)=r.ss_total_in;
        a.ss_eachclass(:,iflip)=r.ss_eachclass_in;
        a.ss_class(:,iflip)=r.ss_class_in;
        for ic=1:length(clist)
            flip_cmatrix(:,:,iflip,ic)=getfield(r,cat(2,'fc_',clist{ic},'_cmatrix_in'));
        end
    end
    for ic=1:length(clist)
        a=setfield(a,cat(2,'fc_',clist{ic},'_cmatrix'),flip_cmatrix(:,:,:,ic));
    end
    results.segflip=a;
    fns=fieldnames(segflip_details);
    for ifn=1:length(fns)
        fn=fns{ifn};
        results.segflip=setfield(results.segflip,fn,getfield(segflip_details,fn));
    end
    results.segflip.msg=results.segflip.aux.msg;
end
%
if (size(results.errormsg,1)>1)
    results.errormsg=results.errormsg(2:end,:); %errormsg initiated with a [] but other messages have been added 
end
opts_used=opts;
return

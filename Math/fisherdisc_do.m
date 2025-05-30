function results=fisherdisc_do(data,tags,opts,outflags)
% results=fisherdisc_do(data,tags,opts,data_out) is a helper
% function for fisherdisc, but can be called as a stand-alone
%
%  data: size(data)=[nfeats,nsamps], i.e., one sample per column, one row per feature
%  tags: size(tags)=[1 nsamps], containing 1 or 2 to tag each sample if insample
%  opts: options (see fisherdisc)
%  outflags: size(outflags)=[1 nsamps], containing 1 if out-of-sample, defaults to zero(1,nsamps)
%    note that tags can be 0 if data is out-of-sample
%    if out-of-sample data is tagged with a valid class, then the class is used
%    to accumulate a confusion matrix
%
% See fisherdisc for descriptions of (most of) results structure
% 
% data_out is out-of-sample data, for classification 
%
% 23Dec20: add counts for each class (count_class_in) and var_equated_in, variance within class assuming equal variance per class
% 23Dec20: add prior, posterior, posterior_diff for 'halfway'
%
%    See also:  FISHERDISC, FISHERDISC_TEST, FISHERDISC_CLASSIFY.
%
nfeat=size(data,1);
nsamps=size(data,2);
if (nargin<=2)
    opts=[];
end
opts=fisherdisc_def(opts);
if (nargin<=3)
    outflags=zeros(1,nsamps);
end
results=[];
results_out=[];
results.errormsg=[];
class_mean=zeros(nfeat,2,2); %nfeat, class, in vs out-of-sample
class_cov=zeros(nfeat,nfeat,2,2); %nfeat, class, in vs out-of-sample
for ix=1:2 
    ptrs_inout{ix}=find(outflags==ix-1);
end
for ig=1:2
    ptrs_sample{ig}=find(tags==ig);
end
results.insample=ptrs_inout{1};
results.outsample=ptrs_inout{2};
for ig=1:2
    for ix=1:2
        ptrs{ig,ix}=intersect(ptrs_sample{ig},ptrs_inout{ix});
        n(ig,ix)=length(ptrs{ig,ix});
        if  n(ig,ix)>0
            class_mean(:,ig,ix)=mean(data(:,ptrs{ig,ix}),2);
            class_cov(:,:,ig,ix)=cov((data(:,ptrs{ig,ix}))',1-opts.sub1);
        end
    end
end
results.class_mean=class_mean(:,:,1);  %added 20 Aug 16
results.class_cov=class_cov(:,:,:,1);  %added 20 Aug 16
totcov=sum(class_cov(:,:,:,1),3);
if (cond(totcov)<opts.condmax)
    w=inv(totcov)*(class_mean(:,2,1)-class_mean(:,1,1));
    if (any(w==Inf))
        w=zeros(nfeat,1);
    end
else
    w=zeros(nfeat,1);
end    
if (sum(w.^2)==0)
    results.discriminant=zeros(1,nfeat);
    results.errormsg='Covariance matrix is singular or nearly so';
else
    results.discriminant=w'/sqrt(sum(w.^2));
end
results.projections=results.discriminant*data;
results.ss_total_in=length(results.insample)*var(results.projections(results.insample),0); %variance around mean
for ig=1:2
    results.count_eachclass_in(ig)=n(ig,1);
    results.ss_eachclass_in(ig)=n(ig,1)*var(results.projections(ptrs{ig,1}),0);
    v(ig)=results.discriminant*class_cov(:,:,ig,1)*results.discriminant';
end
results.var_in=v;
ncorr=n(:,1)-opts.sub1; %class count or class count-1
veq=sum(v*ncorr)/sum(ncorr); %combined estimate of within-class variance
results.var_equated_in=veq;
results.ss_class_in=sum(results.ss_eachclass_in);
clist=opts.classifiers;
if (isempty(clist)) | (sum(w.^2)==0)
    return
end
for ic=1:length(clist)
    cname=clist{ic};
    sname=cat(2,'fc_',cname);
    c=[];
    valid=0;
    switch cname
        case 'halfway'
            c.cutpoint=results.discriminant*mean(class_mean(:,:,1),2);
            c.logprior=[0 0];
            for ig=1:2
                lglk(ig,:)=-0.5*log(veq)-(results.projections-results.discriminant*class_mean(:,ig,1)).^2/(2*veq);
            end
            c.logposterior=lglk; %each Gaussian separately
            c.logposterior_diff=(results.projections-results.discriminant*mean(class_mean(:,:,1),2))*(results.discriminant*diff(class_mean(:,:,1),[],2))/veq;
            c.classes=1+double(results.projections>c.cutpoint); %also, c.classes =1 if lglk(1,:)>lglk(2,:), and c.classes=2 if not
            valid=1;
        case {'mapequal','mapbayes'}
            if (strcmp(cname,'mapbayes'))
                logprior=log(n(:,1));
            else
                logprior=[0 0];
            end
            for ig=1:2
                lglk(ig,:)=logprior(ig)-0.5*log(v(ig))-(results.projections-results.discriminant*class_mean(:,ig,1)).^2/(2*v(ig));
            end
            c.classes=1+double(lglk(2,:)>lglk(1,:));
            valid=1;
            c.logprior=logprior;  %added 20 Aug 16
            c.logposterior=lglk; %added 20 Aug 16
    end
    if (valid)
        for ig_row=1:2
            for ig_col=1:2
                c.cmatrix_in(ig_row,ig_col)=sum(c.classes(ptrs{ig_row,1})==ig_col);
                c.cmatrix_out(ig_row,ig_col)=sum(c.classes(ptrs{ig_row,2})==ig_col);
            end
        end
        fn=fieldnames(c);
        for ifn=1:length(fn)
            results=setfield(results,cat(2,sname,'_',fn{ifn}),getfield(c,fn{ifn}));
        end
    end
end %ic
return

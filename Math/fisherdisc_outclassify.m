function r_out=fisherdisc_outclassify(data,results)
% r_out=fisherdisc_outclassify(data,results)
% classifies a data point (a feature vector) using classifier created by fisherdisc,
%  similar to fisherdisc_classify but:
%    * steps through 'halfway' and 'mapbayes'
%    * fc_mapequal_logposterior: conditional log likelihoods  (i.e., without the log likelihood of the prior), class Gaussians calculated with independent variances
%    * fc_mapbayes_logposterior: log likelihoods taking into account the log likelihood of the prior, class Gaussians calculated with independent variances
%    * fc_halfway_logposterior: log likelihood, with Gaussians for each class assumed to have the same variance, and not taking priors into account
%    * returns the output in the form of a results structure
%
%  data: [nfeats,nsamps], i.e., one sample per column, one row per feature
%  results:  results returned by fisherdisc
%    only some fields are used:
%           results.class_mean
%           results.discriminant
%        just for classifier='halfway':
%           results.fc_halfway_cutpoint  
%           results.var_equated (only present after 23Dec30)
%        just for classifier='mapbayes' or 'mapequal'
%           results.var
%           results.fc_{mapbayes|mapequal}.logprior
%
% See fisherdisc for descriptions of (most of) results structure
%
% 23Dec20: multiple typos in documentation corrected, added fc_halfway_logposterior
%
%    See also:  FISHERDISC, FISHERDISC_DO, FISHERDISC_CLASSIFY.
%
nfeat=size(data,1);
nsamps=size(data,2);
classes=zeros(1,nsamps);
projections=zeros(1,nsamps);
lglk=zeros(2,nsamps);
%
r_out=[];
[classes,projections,lglk]=fisherdisc_classify(data,'halfway',results);
r_out.projections=projections;
r_out.fc_halfway_classes=classes;
r_out.fc_halfway_logposterior=lglk;
[classes,projections,lglk]=fisherdisc_classify(data,'mapbayes',results);
r_out.fc_mapbayes_classes=classes;
r_out.fc_mapbayes_logposterior=lglk;
r_out.fc_mapequal_logposterior=lglk-repmat(results.fc_mapbayes_logprior,1,size(lglk,2));
return


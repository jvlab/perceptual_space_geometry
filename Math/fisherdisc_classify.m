function [classes,projections,lglk]=fisherdisc_classify(data,classifier,results)
% [classes,projections,lglk]=fisherdisc_classify(data,classifier,results)
% classifies a data point (a feature vector) using classifier created by fisherdisc.
%
%  data: [nfeats,nsamps], i.e., one sample per column, one row per feature
%  classifier: one of {'halfway','mapbayes','mapequal']
%  results:  results returned by fisherdisc
%    only some fields are used:
%      results.discriminant
%      results.class_mean
%        just for classifier='halfway':
%           results.fc_halfway_cutpoint
%           results.var_equated (only present after 23Dec30)
%        just for classifier='mapbayes' or 'mapequal'
%           results.var
%           results.fc_{mapbayes|mapequal}.logprior
%  
%  classes: classification (row vector of 1's and 2's)
%  projections:  values of the projections onto the discriminant
%  lglk: log likelihoods of the two classes, for mapbayes or mapequal, and,
%  if var_equated is present in results, then also for halfway
%
% See fisherdisc for descriptions of (most of) results structure
%
%  23Dec30:  add lglk for halfway
%
%    See also:  FISHERDISC, FISHERDISC_DO, FISHERDISC_OUTCLASSIFY.
%
nfeat=size(data,1);
nsamps=size(data,2);
classes=zeros(1,nsamps);
projections=zeros(1,nsamps);
lglk=zeros(2,nsamps);
%
projections=results.discriminant*data;
class_mean=results.class_mean;
switch classifier
    case 'halfway'
        cutpoint=results.fc_halfway_cutpoint;
        classes=1+double(projections>cutpoint);
        if ~isfield(results,'var_equated')
            lglk=NaN(2,nsamps);
        else
            veq=results.var_equated;
            for ig=1:2
                lglk(ig,:)=-0.5*log(veq)-(projections-results.discriminant*class_mean(:,ig,1)).^2/(2*veq);
            end
        end
    case {'mapequal','mapbayes'}
        logprior=results.(cat(2,'fc_',classifier,'_logprior'));
        v=results.var;
        for ig=1:2
            lglk(ig,:)=logprior(ig)-0.5*log(v(ig))-(projections-results.discriminant*class_mean(:,ig)).^2/(2*v(ig));
        end       
        classes=1+double(lglk(2,:)>lglk(1,:));
end
return

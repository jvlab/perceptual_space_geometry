function [bbias,bdebiased,bvar,bsem,blower,bupper]=bootsa(bnaive,bresamps,p)
% [bbias,bdebiased,bvar,bsem,blower,bupper]=bootsa(bnaive,bresamps,p) 
% calculates bootstrap statistics from the fields in a structure
%
%   bnaive: a structure with one or more fields; each field can be an array
%   bresamps: a structure array of the parameters in bnaive, each calculated from a bootstrap resample
%     All fields in bnaive must be present in bresamps
%   p: p-value for calculating two-tailed confidence limits, defaults to 0.05.
%      if p='extreme', then lowest and highest values are used
%
%   bbias: the values that must be subtracted from bnaive
%   bdebiased: the debiased values of bnaive (bnaive-bbias)
%   bvar: the bootstrap variance estimates
%   bsem:  the bootstrap standard errors (sqrt(bvar))
%   blower:  lower confidence limit
%   bupper:  upper confidence limit
%
% See entest.doc for partial documentation.
%
% 16Oct24: add option for p='extreme' to give extreme values of bootstrap sample for lower and upper cl
%  See also:  JACKSA.
%
if (nargin<=2) p=0.05; end
bbias=[];bdebiased=[];bvar=[];bsem=[];blower=[];bupper=[];
nresamps=length(bresamps);
names=fieldnames(bnaive);
if isnumeric(p) %added 16Oct24
    indlow=floor(p/2*nresamps); %conservative choices
    indhigh=ceil((1-p/2)*nresamps)+1;
elseif strcmp(p,'extreme')
    indlow=1;
    indhigh=nresamps;
else
    error('unrecognized entry for p');
end
for iname=1:size(names,1);
   u=deblank(names(iname,:));
   fname=u{1}; %convert from cell array to string
   val=getfield(bnaive,fname);
   shape=size(val);
   r=reshape(val,1,prod(shape));
   clear rd
   for k=1:nresamps
      val=getfield(bresamps(k),fname);
      rd(k,:)=reshape(val,1,prod(shape));
   end
   rdmean=mean(rd,1);
   %do bootstrap calcs
   bias=(rdmean-r); %in contrast to jackknife calculation, no multiplier (ndrop-1)
   debiased=r-bias;
   varv=var(rd,1); %in contrast to jackknife calculation
   sem=varv.^(0.5);
   bbias=setfield(bbias,fname,reshape(bias,shape));
   bdebiased=setfield(bdebiased,fname,reshape(debiased,shape));
   bvar=setfield(bvar,fname,reshape(varv,shape));
   bsem=setfield(bsem,fname,reshape(sem,shape));
   clower=zeros(1,size(rd,2));
   cupper=zeros(1,size(rd,2));
   for d2=1:size(rd,2)
       rds=sort(rd(:,d2));
       if (indlow==0) clower(d2)=-Inf; else clower(d2)=rds(indlow); end
       if (indhigh>nresamps) cupper(d2)=Inf; else cupper(d2)=rds(indhigh); end
   end
   blower=setfield(blower,fname,reshape(clower,shape));
   bupper=setfield(bupper,fname,reshape(cupper,shape));
end
return

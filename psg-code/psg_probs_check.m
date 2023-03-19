%psg_probs_check: check alternative calulations for psg_umi_triplike and related
%
%  See also: PSG_UMI_TRIPLIKE, PBETABAYES_COMPARE, PSG_INEQ_APPLY, PSG_PERMUTES_LOGIC, PSG_INEQ_LOGIC.
%
if ~exist('obs') obs=[3 7;6 8;4 5]; end
if ~exist('params') params=setfields([],{'a','h'},{.3,.1}); end
opts_fast=struct;
opts_fast.if_fast=1;
opts_slow=struct;
opts_slow.if_check=0;
[likrat_fast,lik_fast]=psg_umi_triplike(params,obs,opts_fast);
[likrat_slow,lik_slow]=psg_umi_triplike(params,obs,opts_slow);
%
lik_slow_nolabel=lik_slow; %strip label fields
fields=fieldnames(lik_slow);
for ifield=1:length(fields)
    fn=fields{ifield};
    if ~isempty(strfind(fields{ifield},'label'))
        lik_slow_nolabel=rmfield(lik_slow_nolabel,fields{ifield});
    end
end
%
disp('fast calculation:')
disp(likrat_fast);
disp(sprintf('comparing fast and slow calculation: likrat (%3.0f fields)',length(fieldnames(likrat_fast))));
disp(length(compstruct('fast',likrat_fast,'slow',likrat_slow)))
disp(sprintf('comparing fast and slow calculation: lik    (%3.0f fields)',length(fieldnames(lik_fast))));
disp(length(compstruct('fast',lik_fast,'slow',lik_slow_nolabel)))
%
disp('vectorized calculation')
[likrat_vec,lik_vec]=psg_umi_triplike(params,obs,setfield(opts_fast,'if_vec',1));
disp(sprintf('comparing fast and vectorized calculation: likrat (%3.0f fields)',length(fieldnames(likrat_fast))));
disp(length(compstruct('fast',likrat_fast,'vec',likrat_vec)))
disp(sprintf('comparing fast and vectorized calculation: lik    (%3.0f fields)',length(fieldnames(lik_fast))));
disp(length(compstruct('fast',lik_fast,'vec',lik_vec)))
%
disp('comparing with pbetabayes_compare')
opts_check=struct;
opts_check.if_check=-1;
[likrat_check,lik_check]=psg_umi_triplike(params,obs,opts_check);
%
%    orthant_defs=[0 0 0;1 0 0;0 1 0;1 1 0; 0 0 1;1 0 1;0 1 1;1 1 1]; %hard-coded for speed
%    all_match=[1 0 0 0 0 0 0 1]'; %hard-coded for speed
%efficient generation of hard-coded quantities
ncomps=3;
orthant_defs=[0 1]';
nrows=2;
for icomp=2:ncomps
    orthant_defs=[[orthant_defs,zeros(nrows,1)];[orthant_defs,ones(nrows,1)]];
    nrows=2*nrows;
end
all_match=zeros(2^ncomps,1);
all_match([1 end])=1;
%alternate vectorized calculations of data-dependent quantities
lik_alt=struct;
lik_alt.boundaries=1./(2.^obs(:,2))';
disp(sprintf('checking   boundaries: %12.8f',max(abs(lik_alt.boundaries-lik_fast.boundaries))))
%
a=params.a;
betaln_aa=gammaln(a)*2-gammaln(2*a); %10Feb23
a1s=a+obs(:,1);
a2s=a+obs(:,2)-obs(:,1);
bfs=betainc(0.5,a1s,a2s);
mults=exp(gammaln(a1s)+gammaln(a2s)-gammaln(a1s+a2s)-betaln_aa);
lik_alt.intervals=[bfs.*mults,(1-bfs).*mults]';
disp(sprintf('checking    intervals: %12.8f',max(abs(lik_alt.intervals(:)-lik_fast.intervals(:)))))
%
lik_alt.orthants_sum=prod(mults); %normalization for intervals: each interval adds up to mults(icomp)
disp(sprintf('checking orthants_sum: %12.8f',max(abs(lik_alt.orthants_sum-lik_fast.orthants_sum))))
%
%individual orthants: combinatorial approach to multiplication
orthant_prod=lik_alt.intervals(:,1);
for icomp=2:ncomps
    orthant_prod=[orthant_prod.*lik_alt.intervals(1,icomp);orthant_prod.*lik_alt.intervals(2,icomp)];
end
lik_alt.orthants=reshape(orthant_prod,repmat(2,1,ncomps));
disp(sprintf('checking      orthants: %11.8f',max(abs(lik_alt.orthants(:)-lik_fast.orthants(:)))))
%
disp('computations with psg_umi_triplike');
%
%merge orthants and edges to use for psg_ineq_logic
%
margliks=(1-params.h)*[lik_fast.intervals(1,:);zeros(1,3);lik_fast.intervals(2,:)]+params.h*[zeros(1,3);lik_fast.boundaries(1,:);zeros(1,3)];
blockliks=margliks(:,1);
for nd=2:ncomps
    blockliks=[blockliks*margliks(1,nd);blockliks*margliks(2,nd);blockliks*margliks(3,nd)];
end
blockliks=reshape(blockliks,repmat(3,1,ncomps));
partition=struct;
partition.exclude_trans=psg_ineq_logic(ncomps,'exclude_trans');
partition.exclude_umi=psg_ineq_logic(ncomps,'exclude_umi');
partition.exclude_umi_trans=psg_ineq_logic(ncomps,'exclude_umi_trans');
%
t=sum(sum(sum(blockliks.*(1-partition.exclude_trans))));
u=sum(sum(sum(blockliks.*(1-partition.exclude_umi))));
ut=sum(sum(sum(blockliks.*(1-partition.exclude_umi_trans))));
sn=sum(sum(sum(blockliks.*(partition.exclude_trans)))); %inconsistent with symmetry (assuming trans) = inconsistent with transitivity (assuming sym)
disp(sprintf('trans       error: %15.12f',t-lik_slow.trans))
disp(sprintf('umi         error: %15.12f',u-lik_slow.umi))
disp(sprintf('umi trans   error: %15.12f',ut-lik_slow.umi_trans))
disp(sprintf('sym not     error: %15.12f',sn-lik_slow.sym_not));
%
disp('partition calculation:')
opts_partition=struct;
opts_partition.if_partition=1;
opts_partition.partition=partition;
[likrat_partition,lik_partition]=psg_umi_triplike(params,obs,opts_partition);
%
disp(likrat_partition);
disp(sprintf('comparing fast and partition calculation: likrat (%3.0f fields)',length(fieldnames(likrat_partition))));
disp(length(compstruct('fast',likrat_fast,'partition',likrat_partition)))
disp(sprintf('comparing fast and partition calculation: lik    (%3.0f fields)',length(fieldnames(lik_fast))));
disp(length(compstruct('fast',lik_fast,'partition',lik_partition)))
%
%calculations with psg_ineq_apply
disp('calculation with psg_ineq_apply:')
nc=3;
ineqs=cell(0);
ineqs_umi_trans=psg_ineq_logic(nc,'exclude_umi_trans');
ineqs_trans=psg_ineq_logic(nc,'exclude_trans');
liks=psg_ineq_apply(params,obs,{ineqs_umi_trans,ineqs_trans},psg_permutes_logic(nc,'none'));
likrat_ineq.umi_trans=liks(1)/liks(2);
%
ineqs_sym=psg_ineq_logic(nc,'exclude_sym');
liks=psg_ineq_apply(params,obs,{ineqs_sym,1-ineqs_sym},psg_permutes_logic(nc,'none'));
likrat_ineq.sym=liks(1)/(liks(1)+liks(2));
disp(likrat_ineq);
disp(sprintf('comparing fast and ineq calculation: likrat (%3.0f fields)',length(fieldnames(likrat_fast))));
disp(length(compstruct('fast',likrat_fast,'ineq',likrat_ineq)))



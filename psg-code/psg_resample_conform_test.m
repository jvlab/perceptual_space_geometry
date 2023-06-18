% psg_resample_conform_test: test psg_resample_conform
%
%  See also:  PSG_INEQ_LOGIC, PSG_CONFORM, PSG_RESAMPLE_CONFORM.
%
if ~exist('opts')
    opts=struct;
end
if ~exist('params')
    params=struct;
end
rng('default');
%
params=filldefault(params,'a',0.7);
params=filldefault(params,'h',0.1);
%
if ~exist('nprobs') nprobs=7; end
if ~exist('ndraws') ndraws=10; end
if ~exist('hist_bins') hist_bins=100; end
%
logics={'none3','none6','exclude_sym','exclude_umi_trans','exclude_addtree_trans'};
logics_arg=logics;
logics_arg{1}='none';
logics_arg{2}='none';
ncs=[3 6 3 3 6];
%
partitions=struct;
ncloser_resample=struct;
choice_probs=struct;
ou=struct;
figure;
set(gcf,'NumberTitle','off');
set(gcf,'Name','psg resample conform test');
set(gcf,'Position',[50 50 1200 800]);
nlogics=length(logics);
for ics=1:nlogics
    ln=logics{ics};
    ncomps=ncs(ics);
    disp(sprintf(' %s: ncomps=%1.0f',ln,ncomps));
    partitions.(ln)=psg_ineq_logic(ncomps,logics_arg{ics});
    ntrials=[repmat(4,1,ncomps);repmat(10,1,ncomps);repmat(100,1,ncomps);[1:ncomps]];
    [ncloser_resample.(ln),choice_probs.(ln),ou.(ln)]=psg_resample_conform(params,ntrials,partitions.(ln),nprobs,ndraws,opts); 
    disp(ou.(ln));
    nsucc=nprobs*size(ntrials,1);
    natt=sum(ou.(ln).nattempts(:));
    disp(sprintf('success rate: %3.0f / %6.0f (%8.5f)',nsucc,natt,nsucc/natt));
    %
    subplot(2,nlogics,ics);
    hist(choice_probs.(ln)(:),hist_bins);
    set(gca,'XLim',[0 1]);
    title(cat(2,ln,' probs'),'Interpreter','none');
    xlabel(sprintf('%3.0f/%7.0f',nsucc,natt));
    %empiricial probabilities
    subplot(2,nlogics,nlogics+ics);
    p_obs=ncloser_resample.(ln)./repmat(ntrials,[1 1 nprobs ndraws]);
    hist(p_obs(:),hist_bins);
    set(gca,'XLim',[0 1]);
    title('empirical');
    xlabel(sprintf('%3.0f/%7.0f',nsucc,natt));
end

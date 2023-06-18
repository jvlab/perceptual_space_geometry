% psg_resample_conform_demo: demontrate psg_resample_conform
%
% use a range of values for a, and determine the log likelihood ratio per
% triad or tent, assuming an equal number of trials in each component of
% the triad or tent
%
% generates synthetic data and then computes the log likelihood ratio
% the calculated error brs (llrs_std) are to be scaled by 1/sqrt(ntriads or ntents)
% to obtain the standard devs for an experiment.
%
% Note that error bars (for one triad or tent) are plotted after scaling by
% eb_mult, which defaults to 0.1.
%
% it is assumed that the number of triads or tents are so large so that the
% values of a and h used to synthesize the data are also the maximum-likelihood fit
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
if ~exist('a_list') a_list=[0.1:.1:2]; end
if ~exist('h_list') h_list=[0 .1]; end
if ~exist('ntrials_list') ntrials_list=[0 1 2 4 8 16 32]; end
nh=length(h_list);
%
if ~exist('nprobs') nprobs=1024; end
if ~exist('ndraws') ndraws=4; end
if ~exist('eb_mult') eb_mult=0.1; end
%
%
logics={'exclude_sym','exclude_umi_trans','exclude_addtree_trans'};
logics_denom={'none','exclude_trans','exclude_trans_tent'};
logic_brief={'sym','umi','adt'};
req_discrete=[0 1 0];
%
ncs=[3 3 6];
%
nlogics=length(logics);
llrs=zeros(length(ntrials_list),length(a_list),nh,nlogics);
llrs_std=zeros(length(ntrials_list),length(a_list),nh,nlogics);
figure;
set(gcf,'Position',[50 50 1500 900]);
set(gcf,'NumberTitle','off');
tstring=sprintf('nprobs=%6.0f ndraws=%6.0f',nprobs,ndraws);
set(gcf,'name',tstring);
%
for ilogic=1:nlogics
    ln=logics{ilogic};
    ncomps=ncs(ilogic);
    partitions=psg_ineq_logic(ncomps,logics{ilogic});
    ineq_list{1}=partitions;
    ineq_list{2}=psg_ineq_logic(ncomps,logics_denom{ilogic});
    for ih=1:nh
        if (h_list(ih)>0) | req_discrete(ilogic)==0
            for ia=1:length(a_list)
                params=struct;
                params.a=a_list(ia);
                params.h=h_list(ih);
                disp(sprintf(' %s: ncomps=%1.0f, [a,h]=%6.3f %6.3f',ln,ncomps,params.a,params.h));
                for itrial=1:length(ntrials_list)
                    ntrials=repmat(ntrials_list(itrial),1,ncomps);
                    [ncloser_resample,choice_probs,ou]=psg_resample_conform(params,ntrials,partitions,nprobs,ndraws,opts); 
                    nsucc=nprobs*size(ntrials,1);
                    natt=sum(ou.nattempts(:));
                    %disp(sprintf('success rate: %3.0f / %6.0f (%8.5f)',nsucc,natt,nsucc/natt));
                    ncloser=reshape(ncloser_resample,[ncomps 1 nprobs*ndraws]);
                    % obs: two-column array of integers, of length nc.  Each row corresponds to a dimension of a probability
                    %    cube as set up in psg_ineq_logic.
                    %    obs[ic,1) is number of trials in A is seen as closer to B than to C in r(A;B,C), for the triad corresponding to dimension ic
                    %    obs(ic,2) is the total number of trials
                    %  obs can have a third dimension (of size nsets), to carry out the analysis for multiple sets of observations
                    obs=cat(2,ncloser,repmat(ntrials_list(itrial),[ncomps,1,nprobs*ndraws]));
                    liks=psg_ineq_apply(params,obs,ineq_list);
                    llrs(itrial,ia,ih,ilogic)=mean(log(liks(1,:)./liks(2,:)));
                    llrs_std(itrial,ia,ih,ilogic)=std(log(liks(1,:)./liks(2,:)));
                end
            end %ian
            subplot(nh,nlogics,ilogic+(ih-1)*nlogics);
            %plot(a_list,llrs(:,:,ih,ilogic));
            if req_discrete(ilogic)
                llr_sub=log(h_list(ih));
                ylabel_string='llr-log(h)';
                ylims=[-.2 2];
            else
                llr_sub=0;
                ylabel_string='llr';
                ylims=[-0.55 0];
            end
            errorbar(repmat(a_list,length(ntrials_list),1)',llrs(:,:,ih,ilogic)'-llr_sub,eb_mult*llrs_std(:,:,ih,ilogic)');
            title(sprintf('%s h=%6.3f, eb * %6.4f',logic_brief{ilogic},h_list(ih),eb_mult));
            if req_discrete(ilogic)
                legend(num2str(ntrials_list'),'Location','Best','FontSize',7);
            end
            xlabel('a');
            ylabel(ylabel_string);
            set(gca,'XLim',[0 max(a_list)]);
            set(gca,'YLim',ylims);
        end %test h
    end %ih
end %ilogic
axes('Position',[0.01,0.04,0.01,0.01]); %for text
text(0,0,tstring,'Interpreter','none','FontSize',8);
axis off;


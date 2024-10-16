function [lljit,opts_lljit_used]=psg_lljit(ds,typenames,responses,stim_list,opts_lljit)
% [lljit,opts_lljit_used]=psg_lljit(ds,typenames,responses,stim_list,opts_lljit)
%
% computes log likelihood per trial of a set of responses, given coords for
% the stimuli, and added jitter, and uses this to put a confidence limit, as rms jitter per coordinate,
% on the coordinates.
%
%  Note that log likelihoods are computed log 2, for consistency with SAW software
%
% ds: cell array of coordinate data, ds{k} is the k-dimensional model,
%  ds{k}(istim,idim) is the coordinate of stimulus istim on dimension idim
% typenames: cell{nstims,1}, stimulus labels in ds
% responses: choice data, from choice file, 5 or 6 columns
% stim_list: stimulus labels in choice data, character array, size(stim_list,1)=nstims
% opts_lljit: options
%    if_log: 1 to log, defaults to 0
%    ndraws: number of draws for each jitter, defaults to 100 (only one draw if jitter is 0)
%    jit_list: list of jitter values to determine effect on log likelihood, as rms jitter per dimension
%       if jit_list is only 0, then regression of log likelihood vs jitter
%       is not done, and jit_crits, frac_var are not computed
%    if_frozen: 1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed
%        (called for each dataset, each dimension, and each jitter level, defaults to 1)
%    sigma: standard dev in error model, defaults to 1
%    lik_min: minimum nonzero likelihood value (to prevend underflows), defaults to 10^-6;
%    pvals: significance levels to determine confidence circles of jitter
%    jit_type: 'gaussian' for spherical Gaussian (default)
%              'shell' for spherical shell, same rms jitter per coordinate
%
% lljit: results structure
%    lljit.coord_match: which coordinate matches each value of stim_list
%    lljit.jit_list: jitter list used for the p-values
%    lljit.lk2(ijit,k): mean log likelihoods (base 2) per response, for jitter jitlist(ijit) and model ds{k}
%    lljit.jit_crits(ip,k): Gaussian standard dev (per axis) for p-value opts_lljit.pvals(ip)  and model ds{k}
%    lljit.frac_var(k): fraction of variance accounted for by regression of sqrt(delta-log likelihood) vs. jitter
%    lljit.lk2_bootmean(1,k): mean of bootstrapped ll's (done anaytically) for 0 jitter, model ds{k}
%    lljit.lk2_bootstdv(1,k): stdv of bootstrapped ll's (done anaytically) for 0 jitter, model ds{k}
% opts_lljit_used: options used
%
%  See also: PSG_LLJIT_DEMO, PSG_LLJIT_CRIT.
%
if (nargin<=4)
    opts_lljit=struct;
end
%
log2=log(2);
%
opts_lljit=filldefault(opts_lljit,'if_log',0);
opts_lljit=filldefault(opts_lljit,'ndraws',100);
opts_lljit=filldefault(opts_lljit,'jit_list',[0, 0.01*2.^[0:5]]);
opts_lljit=filldefault(opts_lljit,'if_frozen',1);
opts_lljit=filldefault(opts_lljit,'sigma',1);
opts_lljit=filldefault(opts_lljit,'lik_min',10^-6);
opts_lljit=filldefault(opts_lljit,'pvals',[0.5 0.1 0.05 0.01 0.005 0.001]);
opts_lljit=filldefault(opts_lljit,'jit_type','gaussian');
opts_lljit.jit_list=unique([0 opts_lljit.jit_list]); %first element must be zero
%
opts_lljit_used=opts_lljit;
%
jit_list=opts_lljit.jit_list;
%
lljit=struct;
lljit.jit_list=jit_list;
lljit.warnings=[];
%
% coordinates are in the order of sas{iset}.typenames but choices are in the order of stim_list
%
%first check that all tags in stim_list correspond to coordinates
%
nstims=size(stim_list,1);
coord_match=zeros(nstims,1);
for istim=1:length(stim_list)
    idx=strmatch(stim_list(istim,:),typenames,'exact'); 
    if length(idx)==1
        coord_match(istim)=idx;
    else
        wmsg=sprintf('no unique match for stim label %s',stim_list(istim,:));
        if opts_lljit.if_log==1
            disp(wmsg);
        end
        lljit.warnings=strvcat(lljit.warnings,wmsg);
    end
end
lljit.coord_match=coord_match;
lljit.pvals=opts_lljit.pvals;
ndims=length(ds);
njits=length(jit_list);
npvals=length(opts_lljit.pvals);
ncols=size(responses,2);
if ~ismember(ncols,[5 6])
    wmsg='response matrix does not have 5 or 6 columns';
    if opts_lljit.if_log==1
        disp(wmsg);
    end
    lljit.warnings=strvcat(lljit.warnings,wmsg);
end
nresps=sum(responses(:,end));
dim_list=zeros(1,ndims);
%
if isempty(lljit.warnings)
    lk2=zeros(njits,ndims);
    for idim_ptr=1:ndims %models of each dimension
        idim=size(ds{idim_ptr},2);
        dim_list(idim_ptr)=idim;
        for ijit=1:njits %each level of jitter
            if (opts_lljit.if_frozen~=0)
                rng('default');
                if (opts_lljit.if_frozen<0)
                    rand(1,abs(opts_lljit.if_frozen));
                end
            else
                rng('shuffle');
            end
            %
            jit_rms=jit_list(ijit);
            lk2_bootmean=zeros(1,ndims);
            lk2_bootstdv=zeros(1,ndims);
            ndraws=opts_lljit.ndraws;
            if jit_rms==0
                idraw=1;
            end
            for idraw=1:opts_lljit.ndraws
                lk_list=zeros(opts_lljit.ndraws,1);
                switch opts_lljit.jit_type
                    case 'gaussian'
                        jits=jit_rms*randn(nstims,idim); %jitters, as rms per dimension
                    case 'shell'
                        jits=randn(nstims,idim);
                        jits=jits./repmat(sqrt(sum(jits.^2,2)),1,idim); %normalized
                        jits=sqrt(idim)*jit_rms*jits; %needed to keep rms jitter per dimension equal to jit_rms
                    otherwise
                        error(sprintf('unknown jitter type: %s',opts_lljit.jit_type));
                end
                pcoords=ds{idim_ptr}+jits; %perturbed coordscoor
                %
                %compute distances
                switch ncols
                    case 5
                        %  5 columns: triad comparisons
                        % 'ref                               '
                        % 's1                                '
                        % 's2                                '
                        % 'N(D(ref, s1) > D(ref, s2))        '
                        % 'N_Repeats(D(ref, s1) > D(ref, s2))'
                        %
                        dsq(:,1)=sum((pcoords(coord_match(responses(:,1)),:)-pcoords(coord_match(responses(:,2)),:)).^2,2);
                        dsq(:,2)=sum((pcoords(coord_match(responses(:,1)),:)-pcoords(coord_match(responses(:,3)),:)).^2,2);
                        n1_gt_n2=responses(:,4);
                        nr=responses(:,5);
                    case 6
                        %  6 columns: tetrad comparisons
                        % 's1                              '
                        % 's2                              '
                        % 's3                              '
                        % 's4                              '
                        % 'N(D(s1, s2) > D(s3, s4))        '
                        % 'N_Repeats(D(s1, s2) > D(s3, s4))'
                        dsq(:,1)=sum((pcoords(coord_match(responses(:,1)),:)-pcoords(coord_match(responses(:,2)),:)).^2,2);
                        dsq(:,2)=sum((pcoords(coord_match(responses(:,3)),:)-pcoords(coord_match(responses(:,4)),:)).^2,2);
                        n1_gt_n2=responses(:,5);
                        nr=responses(:,6);
                end
                %compute predicted fraction of n1_gt_n2
                dists=sqrt(dsq);
                p=0.5*(1+erf((dists(:,1)-dists(:,2))/(2*opts_lljit.sigma)));
                q=1-p;
                logp=log(max(p,opts_lljit.lik_min));
                logq=log(min(q,1-opts_lljit.lik_min));
                %compute log likelihood_
                lklist(idraw)=sum((n1_gt_n2).*logp)+sum((nr-n1_gt_n2).*logq);
            end %idraw
            lk2(ijit,idim_ptr)=mean(lklist(1:ndraws))/nresps/log2;
            %if jit_rms=0, also compute bootstrapped expected log likelihood and its variance
            if jit_rms==0
                a=n1_gt_n2./nr;
                bootmean=nr.*(a.*logp+(1-a).*logq);
                lk2_bootmean=sum(bootmean)/nresps/log2; %should match lk2(1,idimptr)
                lljit.lk2_bootmean(1,idim_ptr)=lk2_bootmean;
                bootvar=nr.*a.*(1-a).*(logp-logq).^2; %variance of total log likelihood
                lk2_bootstdv=sqrt(sum(bootvar))/nresps/log2;
                lljit.lk2_bootstdv(1,idim_ptr)=lk2_bootstdv;              
            end
        end %ijit
    end %idim_list
    lljit.dim_list=dim_list;
    lljit.lk2=lk2; %log likelihood (base 2) per response, should match SAW raw LL values for a itter of zero
    lljit.desc_lk2={'log likelihood (base 2) for each rms jitter','model dim ptr'};
    %
    %find the rms jitter for each critical p-value if there are at least one nonzero jitter
    %do regressions: change in log likelihood is expected to be a quadratic function of rms jitter
    %
    if sum(opts_lljit.jit_list>0)>1
        jit_crits=zeros(npvals,ndims);
        fracvar=zeros(1,ndims);
        for idim_ptr=1:ndims
            idim=dim_list(idim_ptr);
            dlk2_tots=(lljit.lk2(:,idim_ptr)-lljit.lk2(1,idim_ptr))*nresps;
            if any(dlk2_tots<0)
                sqrt_dlk2=sqrt(-dlk2_tots); %should be a linear function of jitter
                [b,bint,resids]=regress(sqrt_dlk2,jit_list(:));
                jit_crits(:,idim_ptr)=sqrt(-log(opts_lljit.pvals(:))/log(2))/b;
                fracvar(idim_ptr)=1-sum(resids.^2)/sum(sqrt_dlk2.^2);
            else
                jit_crits(:,idim_ptr)=NaN;
                wmsg=sprintf('jitter improves log likelihood for model of dimension %1.0f',idim);
                if opts_lljit.if_log==1
                    disp(wmsg);
                end
                lljit.warnings=strvcat(lljit.warnings,wmsg);
            end
        end
        lljit.jit_crits=jit_crits;
        lljit.desc_jit_crits={'critical rms jit for each p-value','model dim ptr'};
        lljit.fracvar=fracvar;
        lljit.desc_fracvar='frac variance explained for sqrt(delta-log likelihood) as linear function of jitter';
    end %length(nonzero jit list)>1
end %coordmatch is ok
lljit.nresps=nresps;
%
return

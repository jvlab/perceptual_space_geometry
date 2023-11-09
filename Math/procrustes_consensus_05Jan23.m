function [consensus,znew,ts,details,opts_pcon_used]=procrustes_consensus(z,opts_pcon)
% [consensus,znew,ts,details,opts_pcon_used]=procrustes_consensus(z,opts_pcon)
% carries out a Procrustes consensus calculation; see procrustes_notes.docx.
%
% Basic algorithm is, guess a centroid from one of the datasets;
% align all the datasets with that centroid;
% recompute the centroid
% At each stage, settings in opts_pcon are used re scale, reflection, and offset.
%
% z: input data, size is [npts nds nsets] (nds is number of dimensions)
% opts_pcon: options
%   max_niters: maximum number of iterations
%   max_rmstol: tolerance (rms change in all coordinates of consensus) for termination
%   allow_scale: 1 to allow scaling, defaults to 1
%   allow_reflection: 1 to allow reflection, defaults to 1
%   allow_offset: 1 to allow translation, defaults to 1
%   initialize_set: 1 to nsets, choosse the set (z(:,:,intialize_set) for initial guess
%     and repeated alignment.
%     if 0, then opts_pcon.initialize is an array of size [npts nds] for the initial guess
%     and opts_pcon.alignment is an array of size [npts nds] for alignment.
%     If opts_pcon.alignment is not specified, opts_pcon.initialize is used.
%
% consensus: [npts ndims]: the consensus data
% znew: [npts ndims nsets]: each original dataset, after transformation
% ts: the found transformations.  ts is a cell array, size (1,nsets); 
%    znew(:,:,iset)=ts{iset}.scaling*z(:,:,iset)*ts{iset}.orthog+repmat(ts{iset}.translation,npts,1);
%    ts{iset}.scaling is a scalar, ts{iset}.orthog is [nds nds], ts{iset}.translation is [1 nds]
% details: details of convergence
%    details.ts_cum{k}{iset} is the cumulative transformation found at iteration k
%    details.consensus(:,:,k) is the consensus found at iteration k
%    details.z(:,:,iset,k) is the best fit of each dataset at iteration k
%    details.rms_change(k) is the rms change at iteration k
%    details.rms_dev(:,k) is the rms deviation between all fitted points and the consensus at iteration k
% opts_pcon_used: options used
%
% Notes:
%   Variable names in ts do not match those of Matlab's procrustes.m
%   Matlab's procrustes.m does not allow for eliminating translations
%   Matlab's procrustes.m returns the translation as a matrix, rather than just the row in common
%
% See also:  PROCRUSTES_CONSENSUS_TEST, PROCRUSTES, PSG_PROCRUSTES_DEMO, FILLDEFAULT.
%
if (nargin<2)
    opts_pcon=struct;
end
opts_pcon=filldefault(opts_pcon,'max_niters',100);
opts_pcon=filldefault(opts_pcon,'max_rmstol',10^-5);
opts_pcon=filldefault(opts_pcon,'allow_scale',1);
opts_pcon=filldefault(opts_pcon,'allow_reflection',1);
opts_pcon=filldefault(opts_pcon,'allow_offset',1);
opts_pcon=filldefault(opts_pcon,'initialize_set',1);
opts_pcon=filldefault(opts_pcon,'initial_guess',[]);
opts_pcon=filldefault(opts_pcon,'alignment',opts_pcon.initial_guess);
%
opts_pcon_used=opts_pcon;
%
scaling_token=(opts_pcon.allow_scale==1);
if opts_pcon.allow_reflection==1
    reflection_token='best';
else
    reflection_token=false;
end
%
npts=size(z,1);
nds=size(z,2);
nsets=size(z,3);
znew=z;
%
details=struct;
details.ts_cum=cell(0);
details.consensus=zeros(npts,nds,0);
details.z=zeros(npts,nds,nsets,0);
details.rms_change=zeros(0);
details.rms_dev=zeros(nsets,0);
%
niters=0;
rms=Inf;
if opts_pcon.initialize_set>0
    initial_guess=z(:,:,opts_pcon.initialize_set);
    alignment=initial_guess;
else
    initial_guess=opts_pcon.initial_guess;
    alignment=opts_pcon.alignment;
end
consensus=initial_guess; %initial guess of consensus
while (niters<opts_pcon.max_niters & rms>opts_pcon.max_rmstol)
    niters=niters+1;
    zz=zeros(npts,nds,nsets);
    ts_cum=cell(1,nsets);
    for iset=1:nsets
        [d,zz(:,:,iset),t]=procrustes(consensus,z(:,:,iset),'Scaling',scaling_token,'Reflection',reflection_token);
        ts_cum{iset}.scaling=t.b;
        ts_cum{iset}.orthog=t.T;
        if (opts_pcon.allow_offset==0) %remove offset if reqeueted
            zz(:,:,iset)=zz(:,:,iset)-t.c;
            ts_cum{iset}.translation=zeros(1,nds);
        else
            ts_cum{iset}.translation=t.c(1,:);
        end
    end
    consensus_new=mean(zz,3);
    %realign 
    [d,consensus_new]=procrustes(alignment,consensus_new,'Scaling',scaling_token,'Reflection',reflection_token);
    details.ts_cum{niters}=ts_cum;
    details.consensus(:,:,niters)=consensus_new;
    rms=sqrt(mean((consensus_new(:)-consensus(:)).^2));
    details.z(:,:,:,niters)=zz;
    details.rms_change(1,niters)=rms;
    consensus=consensus_new;
    details.rms_dev(:,niters)=reshape(sqrt(mean(mean((zz-consensus_new).^2,1),2)),[nsets 1]);
end   
znew=zz;
ts=ts_cum;
return

function [consensus,znew,ts,details,opts_pcon_used]=procrustes_consensus(z,opts_pcon)
% [consensus,znew,ts,details,opts_pcon_used]=procrustes_consensus(z,opts_pcon)
% carries out a Procrustes consensus calculation; see procrustes_notes.docx.
%
% Basic algorithm is  guess a centroid from one of the datasets;
% * Step 1: align all the datasets with the centroid
% * Step 2: recompute the centroid dataset
% * Step 3: align the centroid dataset as specified by opts_pcon.alignment, to prevent drift
%  "Centroid" means the all point-by-point centroids, considering, for each point the mean
%  of the datasets for which overlap(ipt,:)=1)
%
% At each step, settings in opts_pcon are used for scale, reflection, and offset.
%
% z: input data, size is [npts nds nsets] (nds is number of dimensions)
% opts_pcon: options
%   max_niters: maximum number of iterations
%   max_rmstol: tolerance (rms change in all coordinates of consensus) for termination
%   allow_scale: 1 to allow scaling, defaults to 1
%   allow_reflection: 1 to allow reflection, defaults to 1
%   allow_offset: 1 to allow translation, defaults to 1
%   initialize_set: 1 to nsets, choose the set (z(:,:,intialize_set) for initial guess and repeated alignment.
%     if 0, then opts_pcon.initialize is an array of size [npts nds] for the initial guess
%     and opts_pcon.alignment is an array of size [npts nds] for alignment.
%     If opts_pcon.alignment is not specified, opts_pcon.initialize is used.
%   overlaps: [npts nsets]: binary array, indicating which points should be used in the calculation
%     If omitted, defaults to ones(npts,nsets)
%   if_justcheck: just check if overlap matrix is ok (results in details.warnings)
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
%    details.rms_dev(:,k) is the rms deviation between all fitted points and the
%       consensus at iteration k (points without overlaps are omitted)
%    details.zz_check_diff is a check for the difference between
%       transformed points computed by Procrustes and by hand
%    details.warnings: warnings if not enough overlaps
%    details.overlap_pairs(nsets,nsets: number of pairwise overlaps
%    details.overlap_totals(npts,1): number of overlaps for each data point
%    details.initialize_use(npts,1): which dataset is used to initialize each point
%
% opts_pcon_used: options used
%
% Notes:
%   Variable names in ts do not match those of Matlab's procrustes.m
%   Matlab's procrustes.m does not allow for eliminating translations
%   Matlab's procrustes.m returns the translation as a matrix, rather than just the row in common
%
% 06Nov23: begin to work on version with partial overlaps
%
% See also:  PROCRUSTES_CONSENSUS_TEST, PROCRUSTES, PSG_PROCRUSTES_DEMO, FILLDEFAULT, PROCRUSTES_CONSENSUS_PTL_TEST,
%    CONNCOMP, GRAPH.
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
opts_pcon=filldefault(opts_pcon,'overlaps',[]);
opts_pcon=filldefault(opts_pcon,'if_justcheck',0);
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
if isempty(opts_pcon.overlaps)
    opts_pcon.overlaps=ones(npts,nsets);
end
znew=z;
%
details=struct;
details.ts_cum=cell(0);
details.consensus=zeros(npts,nds,0);
details.z=zeros(npts,nds,nsets,0);
details.rms_change=zeros(0);
details.rms_dev=zeros(nsets,0);
details.overlap_pairs=opts_pcon.overlaps'*opts_pcon.overlaps;
details.overlap_totals=sum(opts_pcon.overlaps,2);
details.warnings=[];
if any(details.overlap_totals==0)
    wmsg='some data points never occur in the overlap matrix';
    warning(wmsg);
    details.warnings=strvcat(details.warnings,wmsg);
end
any_ovlp=double(details.overlap_pairs>0);
graph_ovlp=graph(any_ovlp-diag(diag(any_ovlp)));
conncomps=conncomp(graph_ovlp); %is the overlap graph connected?
if any(conncomps>1)
    wmsg='overlap graph is not connected';
    details.warnings=strvcat(details.warnings,wmsg);
end
if (opts_pcon.if_justcheck)
    consensus=[];
    ts=[];
    return
end
%
niters=0;
rms=Inf;
first_nz=zeros(npts,1);
if opts_pcon.initialize_set>0
    %if opts_pcon.overlaps is not all 1's, find first dataset with a nonzero value of overlaps(iset,:)
    overlaps_permute=opts_pcon.overlaps(:,1+mod(opts_pcon.initialize_set-1+[0:nsets-1],nsets));
    for ipt=1:npts
        first_nz(ipt,1)=mod(opts_pcon.initialize_set-1+min(find(overlaps_permute(ipt,:)>0)-1),nsets)+1;
        initial_guess(ipt,:)=z(ipt,:,first_nz(ipt,1));
    end
    alignment=initial_guess;
else
    initial_guess=opts_pcon.initial_guess;
    alignment=opts_pcon.alignment;
end
details.initial_guess=initial_guess;
details.initialize_use=first_nz;
consensus=initial_guess; %initial guess of consensus
while (niters<opts_pcon.max_niters & rms>opts_pcon.max_rmstol)
    niters=niters+1;
    zz=zeros(npts,nds,nsets);
    zznan=nan(npts,nds,nsets); %will have nan's for all non-selected elements
    ts_cum=cell(1,nsets);
    for iset=1:nsets
        select=find(opts_pcon.overlaps(:,iset)>0);
        %do a Procrustes alignment for all points that have an overlap [Step 1]
        [d,zz_sel,t]=procrustes(consensus(select,:),z(select,:,iset),'Scaling',scaling_token,'Reflection',reflection_token);
        ts_cum{iset}.scaling=t.b;
        ts_cum{iset}.orthog=t.T;       
        % Z = TRANSFORM.b * Y * TRANSFORM.T + TRANSFORM.c.
        c_row=t.c(1,:);
        zz(:,:,iset)=t.b*z(:,:,iset)*t.T+repmat(c_row,npts,1); %but transform all the points
        details.zz_check_diff(niters+1,iset)=max(max(abs(zz(select,:,iset)-zz_sel))); 
        if (opts_pcon.allow_offset==0) %remove offset if requested
            zz(:,:,iset)=zz(:,:,iset)-repmat(c_row,npts,1);
            ts_cum{iset}.translation=zeros(1,nds);
        else
            ts_cum{iset}.translation=c_row;
        end
        zznan(select,:,iset)=zz(select,:,iset); %non-selected elements remain nans
    end
    consensus_new=mean(zznan,3,'omitnan'); %only average non-nans [Step 2]
    %realign [Step 3]
    [d,consensus_new]=procrustes(alignment,consensus_new,'Scaling',scaling_token,'Reflection',reflection_token);
    details.ts_cum{niters}=ts_cum;
    details.consensus(:,:,niters)=consensus_new;
    rms=sqrt(mean((consensus_new(:)-consensus(:)).^2));
    details.z(:,:,:,niters)=zz;
    details.rms_change(1,niters)=rms;
    consensus=consensus_new;
%    details.rms_dev(:,niters)=reshape(sqrt(mean(mean((zz-consensus_new).^2,1),2)),[nsets 1]);
%    details.rms_dev_orig(:,niters)=reshape(sqrt(mean(mean((zz-consensus_new).^2,1),2)),[nsets 1]);
    for iset=1:nsets
        details.rms_dev(iset,niters)=sqrt(mean(mean((zznan(:,:,iset)-consensus_new).^2,1,'omitnan'),2));
    end
end   
znew=zz;
ts=ts_cum;
return

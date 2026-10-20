%procrustes_consensus_test2: test procrustes consensus routine, with
%emphasis on missing data, and inconsistency found 14Oct26 between consensus, and consensus deermined from components, when missing data are present:
%
%When if_pca=0 and offset=1 and scaling=1, ther is a constant difference (different for each coordinate)
%When if_pca=1 and if_scale and offset are 0, there is a constant difference.
%When if_pca=1 and if_scale=1 and offsert=0, there is a constant difference
%When if_pca=1 and if_offset=1, ther difference depends on which datawsets are missing, whether or not scaling is allowed)
%
% in contrast to procrustes_test:
%  * original data need not be centered
%  * noise jitter added at end
%  * smaller default tolerance (10^-6)
%
% need to add missing data mask
%
%   See also: PROCRUSTES_CONSENSUS, PROCRUSTES_CONSENSUS_TEST.
%
if ~exist('if_frozen') if_frozen=1; end
if ~exist('nds') nds=5; end %number of dimensions
if ~exist('nsets') nsets=3; end %number of data sets, don't change this
if ~exist('if_reflection') if_reflection=1; end % allow negation?
if ~exist('sd_logscale') sd_logscale=1; end %standard dev for log of scale multiplier
if ~exist('sd_addnoise') sd_addnoise=0.1; end %standard dev for additive noise on each coordinate
if ~exist('sd_offset') sd_offset=0.5; end %standard dev for offset of pre-noise centroid on each coordinate
%
if_frozen=getinp('1 for frozen randomization','d',[0 1],if_frozen);
sd_logscale=getinp('s.d. for log of scale multiplier','f',[0 3],sd_logscale);
sd_addnoise=getinp('s.d. for additive noise on each coordinate','f',[0 3],sd_addnoise);
sd_offset=getinp('s.d. for offset (of pre-noise centroid) on each coordinate','f',[0 5],sd_offset);
%
if (if_frozen)
    rng('default');
end
%
%create test data
%
if ~exist('npts_common') npts_common=12; end %points in common for all datasets
if ~exist('npts_missing') npts_missing=[5 6 7]; end %points missing from one dataset
if ~exist('npts_only') npts_only=[2 4 3]; end %points just in one dataset
npts=npts_common+sum(npts_missing)+sum(npts_only);
%
masks=ones(npts,nsets);
for iset=1:nsets
    ptr=npts_common;
    for jset=1:nsets
        if jset==iset
            masks(ptr+[1:npts_missing(jset)],iset)=0;
        end
        ptr=ptr+npts_missing(jset);
    end
    for jset=1:nsets
        if (jset~=iset)
            masks(ptr+[1:npts_only(jset)],iset)=0;
        end
        ptr=ptr+npts_only(jset);
    end
end
figure;
set(gcf,'Name','overlap pattern');
spy(masks);
title('masks');
xlabel('set');
ylabel('data point');disp('masks (transposed)');
axis tight;
%
%create transformations in a way that preserves random number generator indepedent of options
%
%compute random rotations
rotations=zeros(nds,nds,nsets);
for iset=1:nsets
     rotations(:,:,iset)=randorthu(nds,1);
end
%compute random reflections, making sure that there is at least one sign flip and one unflipped if reflections requested
reflects=(-1).^[randperm(nsets)-1];
if (~if_reflection)
    reflects=ones(1,nsets);
end
scales=exp(sd_logscale*randn(1,nsets)); %random scales
offset=sd_offset*randn(1,nds,nsets); %random offsets
addnoise=sd_addnoise*randn(npts,nds,nsets); %additive noise
%
%random data
%
z_org=randn(npts,nds); %random points
t_expect=cell(nsets,1); %expected transforms, by inverting the random forward transform created here
for iset=1:nsets
    z_rot(:,:,iset)=z_org*rotations(:,:,iset); %rotate
    z_ref(:,:,iset)=z_rot(:,:,iset)*reflects(iset); %possibly reflect
    z_scl(:,:,iset)=z_ref(:,:,iset)*scales(iset); %scale
    z_off(:,:,iset)=z_scl(:,:,iset)+repmat(offset(1,:,iset),[npts,1,1]); %offset all the points
    %
    A=inv(rotations(:,:,iset));
    t_expect{iset}.translation=-offset(1,:,iset)*A*reflects(iset)/scales(iset);
    t_expect{iset}.orthog=A*reflects(iset);
    t_expect{iset}.scaling=1/scales(iset);
end
z_adn=z_off+addnoise;
disp('test data created.');
%
%do  consensus runs and check consistency
%
if ~exist('max_niters') max_niters=50; end
if ~exist('max_rmstol') max_rmstol=10^-6; end
%if ~exist('allow_scale') allow_scale=double(sd_logscale~=0); end
%if ~exist('allow_reflection') allow_reflection=if_reflection; end
%if ~exist('allow_offset') allow_offset=double(sd_offset~=0); end
%
max_niters=getinp('maximum number of iterations','d',[1 10^4],max_niters);
max_rmstol=getinp('maximum rms change','f',[0 1],max_rmstol);
%allow_scale=getinp('1 to allow scale changes','d',[0 1],allow_scale);
%allow_reflection=getinp('1 to allow reflections','d',[0 1],allow_reflection);
%allow_offset=getinp('1 to allow offset','d',[0 1],allow_offset);
%if_normscale=getinp('1 to normalize scale (only relevant if allow offset=1)','d',[0 1],allow_offset);
%
masks_reshape=repmat(reshape(masks,[npts,1,nsets]),[1 nds 1]);
run=0;
consens_devs=zeros(1,2^5);
transform_devs=zeros(1,2^5);
for if_mask=0:1
    z_data=z_adn;
    if (if_mask)
        z_data(masks_reshape==0)=NaN;
    end
    for allow_scale=0:1
        for allow_offset=0:1
            for allow_reflection=0:1
                for if_normscale=0:1
                    opts_pcon=struct();
                    opts_pcon.max_niters=max_niters;
                    opts_pcon.max_rmstol=max_rmstol;
                    opts_pcon.allow_scale=allow_scale;
                    opts_pcon.allow_reflection=allow_reflection;
                    opts_pcon.allow_offset=allow_offset;
                    opts_pcon.if_normscale=if_normscale;
                    %
                    [consens,z_consensus,transforms,details,opts_pcon_used]=procrustes_consensus(z_data,opts_pcon);
                    %
                    disp('****');
                    disp(sprintf('summary for if_mask=%2.0f allow_scale=%2.0f allow_offset=%2.0f allow_reflection=%2.0f if_normscale=%2.0f',if_mask,allow_scale,allow_offset,allow_reflection,if_normscale));
                    rms_dev_consensus=sqrt(reshape(sum(sum((repmat(consens,[1 1 nsets])-z_consensus).^2,1,'omitnan'),2,'omitnan'),[1 nsets 1])/nds/npts);
                    rms_dev_orig=sqrt(reshape(sum(sum((repmat(consens,[1 1 nsets])-z_data).^2,1,'omitnan'),2,'omitnan'),[1 nsets 1])/nds/npts);
                    disp('niters, original rms devs from consensus by dataset, final dev from consensus by dataset');
                    disp(sprintf('niters: %5.0f',size(details.z,4)));
                    disp(rms_dev_orig)
                    disp(rms_dev_consensus);
                    %disp('by iteration')
                    %disp(details.rms_dev');
                    %recalculate consensus from components
                    consens_recalc=mean(z_consensus,3,'omitnan');
                    consens_dev=max(abs(consens(:)-consens_recalc(:)));
                    disp(sprintf('max dev between consensus and recalculated consensus:                           %12.7f',consens_dev));
                    %recalculate components from transformations
                    z_consensus_recalc=zeros(npts,nds,nsets);
                    for iset=1:nsets
                        z_consensus_recalc(:,:,iset)=transforms{iset}.scaling*z_data(:,:,iset)*transforms{iset}.orthog+repmat(transforms{iset}.translation,npts,1);
                    end
                    transform_dev=max(abs(z_consensus(:)-z_consensus_recalc(:)));
                    disp(sprintf('max dev between transformed components and recalculated transformed components: %12.7f',transform_dev));
                    run=run+1;
                    consens_devs(run)=consens_dev;
                    transform_devs(run)=transform_dev;
                end %if_normscale
            end %allow_reflection
        end %allow_offset
    end %allow_scale
end %if_mask
disp('summary');
disp('consensus recalc devations (rows: allow_scale, groups: if_mask)');
disp('.............allow offset=0............     .............allow_offset=1........');
disp('     reflect=0              reflect=1            reflect=0             reflect=1');
disp('normscale=0    =1      normscale=0    =1    normscale=0    =1   normscale=0   =1');
disp(permute(reshape(consens_devs,[8,2,2]),[2 1 3]));
disp(sprintf('max: %12.7f',max(consens_devs(:))));
disp('transform recalc devations (rows: allow_scale, groups: if_mask)');
disp(permute(reshape(transform_devs,[8,2,2]),[2 1 3]));
disp(sprintf('max: %12.7f',max(transform_devs(:))));

%procrustes_consensus_test: test procrustes consensus routine
%
% z_org(npts,nds): Gaussian random points
% z_ctr(npts,nds):  z_orig after centering
% z_rot(npts,nds,nsets):  z_ctr after random rotations
% z_ref(npts,nds,nsets): z_rot after (optional) reflections
% z_scl(npts,nds,nsets): z_ref after (optional) multiplicative scaling
% z_adn(npts,nds,nsets): z_scl after (optional) additive noise
% z_off(npts,nds,nsets): z_adn after (optional) random offsets
%
% note that random draws are made for flip, whether or not they are used,
% to ensure reproducibility
%
% 06Nov23: show summary of consensus
%
%   See also: PROCRUSTES, RANDORTHU, PROCRUSTES_CONSENSUS.
%
if ~exist('if_frozen') if_frozen=1; end
if ~exist('npts') npts=30; end %number of sample points
if ~exist('nds') nds=5; end %number of dimensions
if ~exist('nsets') nsets=4; end %number of data sets
if ~exist('if_reflection') if_reflection=1; end % allow negation?
if ~exist('sd_logscale') sd_logscale=1; end %standard dev for log of scale multiplier
if ~exist('sd_addnoise') sd_addnoise=0.2; end %standard dev for additive noise on each coordinate
if ~exist('sd_offset') sd_offset=0.5; end %standard dev for offset of pre-noise centroid on each coordinate
%
if_frozen=getinp('1 for frozen randomization','d',[0 1],if_frozen);
npts=getinp('number of data points','d',[5 1000],npts);
nds=getinp('number of dimensions','d',[1 10],nds);
nsets=getinp('number of data sets','d',[2 10],nsets);
if_reflection=getinp('1 to allow reflections','d',[0 1],if_reflection);
sd_logscale=getinp('s.d. for log of scale multiplier','f',[0 3],sd_logscale);
sd_addnoise=getinp('s.d. for additive noise on each coordinate','f',[0 3],sd_addnoise);
sd_offset=getinp('s.d. for offset (of pre-noise centroid) on each coordinate','f',[0 5],sd_offset);
%
%create test data
%
if (if_frozen)
    rng('default');
end
z_org=randn(npts,nds); %random points
%center the data
z_ctr=z_org-repmat(mean(z_org,1),npts,1);
%apply random rotations
rotations=zeros(nds,nds,nsets);
for iset=1:nsets
    rotations(:,:,iset)=randorthu(nds,1);
    z_rot(:,:,iset)=z_ctr*rotations(:,:,iset);
end
%apply random reflections, making sure that there is at least one sign flip and one unflipped
reflects=(-1).^[randperm(nsets)-1];
if (if_reflection)
    for iset=1:nsets
        z_ref(:,:,iset)=z_rot(:,:,iset)*reflects(iset);
    end
else
    z_ref=z_rot;
end
%apply random scale
scales=exp(sd_logscale*randn(1,nsets));
for iset=1:nsets
    z_scl(:,:,iset)=z_ref(:,:,iset)*scales(iset);
end
%apply additive noise
addnoise=sd_addnoise*randn(npts,nds,nsets);
z_adn=z_scl+addnoise;
%apply random offsets
offset=sd_offset*randn(1,nds,nsets);
z_off=z_adn+repmat(offset,[npts,1,1]);
disp('test data created.');
%
%do a consensus run
%
if ~exist('max_niters') max_niters=50; end
if ~exist('max_rmstol') max_rmstol=10^-5; end
if ~exist('allow_scale') allow_scale=double(sd_logscale~=0); end
if ~exist('allow_reflection') allow_reflection=if_reflection; end
if ~exist('allow_offset') allow_offset=double(sd_offset~=0); end
%
max_niters=getinp('maximum number of iterations','d',[1 10^4],max_niters);
max_rmstol=getinp('maximum rms change','f',[0 1],max_rmstol);
allow_scale=getinp('1 to allow scale changes','d',[0 1],allow_scale);
allow_reflection=getinp('1 to allow reflections','d',[0 1],allow_reflection);
allow_offset=getinp('1 to allow offset','d',[0 1],allow_offset);
if_normscale=getinp('1 to normalize scale (only relevant if allow offset=1)','d',[0 1],0);
%
opts_pcon=struct();
opts_pcon.max_niters=max_niters;
opts_pcon.max_rmstol=max_rmstol;
opts_pcon.allow_scale=allow_scale;
opts_pcon.if_normscale=if_normscale;
opts_pcon.allow_reflection=allow_reflection;
opts_pcon.allow_offset=allow_offset;
%
[consens,z_consensus,transforms,details,opts_pcon_used]=procrustes_consensus(z_off,opts_pcon);
%
rms_dev_consensus=sqrt(reshape(sum(sum((repmat(consens,[1 1 nsets])-z_consensus).^2,1),2),[1 nsets 1])/nds/npts);
rms_dev_orig=sqrt(reshape(sum(sum((repmat(consens,[1 1 nsets])-z_off).^2,1),2),[1 nsets 1])/nds/npts);
disp('original rms devs');
disp(rms_dev_orig)
disp('by iteration')
disp(details.rms_dev');
disp('final')
disp(rms_dev_consensus);

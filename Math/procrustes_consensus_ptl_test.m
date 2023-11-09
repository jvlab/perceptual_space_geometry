%procrustes_consensus_ptl_test: test procrustes consensus routine
% for datasets that only partially overlap
%
% derived from procrustes_consensus_test.m
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
% ****need better tests for adequacy of number of params, here and in procrustes_consensus
%
% 06Nov23: show summary of consensus
% 06
%   See also: PROCRUSTES, RANDORTHU, PROCRUSTES_CONSENSUS, PROCRUSTES_CONSENSUS_TEST.
%
if ~exist('if_frozen') if_frozen=1; end
if ~exist('npts') npts=30; end %number of sample points
if ~exist('nds') nds=5; end %number of dimensions
if ~exist('nsets') nsets=4; end %number of data sets
if ~exist('if_reflection') if_reflection=1; end % allow negation?
if ~exist('sd_logscale') sd_logscale=1; end %standard dev for log of scale multiplier
if ~exist('sd_addnoise') sd_addnoise=0.2; end %standard dev for additive noise on each coordinate
if ~exist('sd_offset') sd_offset=0.5; end %standard dev for offset of pre-noise centroid on each coordinate
if ~exist('p_ovlp') p_ovlp=0.5; end % raw probability of being included in an overlap
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
%set up an array of partial overlaps
%
ok_ovlp=0;
ovlp_array=zeros(npts,nsets);
while (ok_ovlp==0)
    for ipt=1:npts
        while (sum(ovlp_array(ipt,:))<2)
            ovlp_array(ipt,:)=double(rand(1,nsets)<p_ovlp);
        end
    end
    ovlp_pairset=ovlp_array'*ovlp_array;
    ok_ovlp=all(ovlp_pairset>nds); %check that every dataset overlaps with each other set in at least nds points        
end
novlps=sum(ovlp_array(:));
ovlp_mask=repmat(reshape(ovlp_array,[npts 1 nsets]),[1 nds 1]);
disp('overlap array')
disp(ovlp_array');
disp('pairwise overlaps')
disp(ovlp_pairset);
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
%
if ~exist('opts_pcon')
    opts_pcon=struct();
end
opts_pcon.max_niters=max_niters;
opts_pcon.max_rmstol=max_rmstol;
opts_pcon.allow_scale=allow_scale;
opts_pcon.allow_reflection=allow_reflection;
opts_pcon.allow_offset=allow_offset;
%%
%standard consensus run
%
disp(' ');
disp('standard consensus run')
[consens,z_consensus,transforms,details,opts_pcon_used]=procrustes_consensus(z_off,opts_pcon);
%
rms_dev_orig=sqrt(reshape(mean(mean((repmat(consens,[1 1 nsets])-z_off).^2,1),2),[1 nsets 1]));
rms_dev_consensus=sqrt(reshape(mean(mean((repmat(consens,[1 1 nsets])-z_consensus).^2,1),2),[1 nsets 1]));
disp('original rms devs');
disp(rms_dev_orig)
disp('by iteration')
disp(details.rms_dev');
disp('final rms devs')
disp(rms_dev_consensus);
%
novmode=3;
consens_ov=cell(1,novmode);
z_consensus_ov=cell(1,novmode);
transforms_ov=cell(1,novmode);
details_ov=cell(1,novmode);
opts_pcon_ov_used=cell(1,novmode);
rms_dev_orig_ov=zeros(novmode,nsets);
rms_dev_consensus_ov=zeros(novmode,nsets);
rms_dev_orig_ov_ovonly=zeros(novmode,nsets);
rms_dev_consensus_ov_ovonly=zeros(novmode,nsets);
%
opts_pcon_ov=setfield(opts_pcon,'overlaps',ovlp_array);
%
z_off_nan=z_off;
z_off_nan(ovlp_mask==0)=NaN;
for im=1:3
    disp(' ');
    switch im
        case 1
            disp('consensus run with non-overlaps present but ignored')
            z_off_use=z_off;
        case 2
            disp('consensus run with non-overlaps set to zero')
            z_off_use=z_off;
            z_off_use(ovlp_mask==0)=0;
        case 3
            disp('consensus run with non-overlaps set to NaN');
            z_off_use=z_off_nan;
    end
    [consens_ov{im},z_consensus_ov{im},transforms_ov{im},details_ov{im},opts_pcon_ov_used{im}]=procrustes_consensus(z_off_use,opts_pcon_ov);
    %
    rms_dev_orig_ov(im,:)=sqrt(reshape(mean(mean((repmat(consens_ov{im},[1 1 nsets])-z_off_use).^2,1),2),[1 nsets 1]));
    rms_dev_consensus_ov(im,:)=sqrt(reshape(mean(mean((repmat(consens_ov{im},[1 1 nsets])-z_consensus_ov{im}).^2,1),2),[1 nsets 1]));
    disp('original rms devs');
    disp(rms_dev_orig_ov(im,:))
    disp('by iteration')
    disp(details_ov{im}.rms_dev');
    disp('final rms devs')
    disp(rms_dev_consensus_ov(im,:));
    rms_dev_orig_ov_ovonly(im,:)=sqrt(reshape(mean(mean(((repmat(consens_ov{im},[1 1 nsets])-z_off_nan).^2),1,'omitnan'),2),[1 nsets 1]));
    z_consensus_ov_nan=z_consensus_ov{im};
    z_consensus_ov_nan(ovlp_mask==0)=NaN;
    rms_dev_consensus_ov_ovonly(im,:)=sqrt(reshape(mean(mean(((repmat(consens_ov{im},[1 1 nsets])-z_consensus_ov_nan).^2),1,'omitnan'),2),[1 nsets 1]));
    disp('original rms devs, overlaps only');
    disp(rms_dev_orig_ov_ovonly(im,:))
    disp('final rms devs, overlaps only')
    disp(rms_dev_consensus_ov_ovonly(im,:));
end

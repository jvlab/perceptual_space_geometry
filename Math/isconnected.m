function [isc,aux]=isconnected(binreg,nbrvec,mindex,cmatrix_full)
% [isc,aux]=isconnected(binreg,nbrvec) determines whether a multidimensional binary region is connected
%
%   This is typically much slower than matlab's bwlabeln, but it does not depend on the image processing toolbox.
%
% binreg: multidimensional array, 0 is background and 1 is region
% nbrvec: array of size [nbrs ndims(binreg)], indicating connectivity
%    OR, if nbrvec is a scalar, then the contact dimensionality required for connectivity
%    e.g., if ndims=3 and nbrvec=2, then it is connectivity across the 6 faces
%          if ndims=3 and nbrvec=1, then it is connectivity across the 6 faces and 12 edges
%    Note that generating nbrvec from a scalar is time-consuming, so it is
%    best if it is passed.  It can be pre-calculated by
%    [isc,aux]=isconnected(zeros(repmat(2,1,ndims)),dcontact)
% mindex: multidimensional indexing array, as computed by mindex_make, or blank or omitted
% cmatrix_full: full connectivity matrix for binreg=1 (see opt_tag_setup), optional: computed if not supplied
% 
% isc: 1 if connected, 0 if not.
% aux: auxiliary outputs
%      aux.mindex: result of mindex_make
%      aux.mults: result of mindex_make
%      aux.nbrvec: nbrvec, as provided
%      aux.nbrvec_used: nbrvec, as expanded
%      aux.coord_ptrs: pointers to coordinates (i.e., the rows of mindex that are part of the region)
%      aux.coords: multi-indices of coordinates (i.e., contents of the rows pointed to by coord_ptrs)
%      aux.cmatrix: connectivity matrix (rows and columns correspond to the rows of coord_ptrs)      
%      aux.niters from isgraphc
%      aux.msg from isgraphc
%
% 10Nov25: fix cocumentation (this determines whether region is connected, not "simply connected")
%   See also:  MINDEX_MAKE, ISGRAPHC, OPT_TAG_SETUP, OPT_TAG_ISCONNECTED, BWLABELN.
%
if nargin<=2
    mindex=[];
end
if isempty(mindex)
    [mindex,mults]=mindex_make(size(binreg));
else
    dprod=[1 cumprod(size(binreg))];
    mults=dprod(1:end-1)';
end
if (nargin<=3)
    cmatrix_full=[];
end
r=ndims(binreg);
aux=[];
aux.mindex=mindex;
aux.mults=mults;
aux.coord_ptrs=find(binreg(:));
aux.coords=mindex(aux.coord_ptrs,:);
npts=length(aux.coord_ptrs);
%
if isscalar(nbrvec)
    nbrvec_used=zeros(0,r);
    %generate neighbor-definition vectors if not supplied
    for d=nbrvec:(r-1)
        dsel=nchoosek([1:r],r-d);
        offsets=[0:size(dsel,1)-1]*r;
        posits=repmat(offsets',1,r-d)+dsel;
        %
        dirs=[];
        dseg=zeros(r,size(dsel,1));
        signs=2*mindex_make(repmat(2,1,r-d))-1; %array of all possible 2^(r-d) sign choices
        for isign=1:size(signs,1)
            for pcol=1:(r-d)
                dseg(posits(:,pcol))=signs(isign,pcol);
            end
            dirs=[dirs;dseg'];
        end
        %
        nbrvec_used=[nbrvec_used;dirs];
    end
else
    nbrvec_used=nbrvec;
end
aux.nbrvec=nbrvec;
aux.nbrvec_used=nbrvec_used;
%
%generate connectivity matrix
%
if isempty(cmatrix_full)
    nnbrs=size(nbrvec_used,1);
    max_mindex=repmat(size(binreg),nnbrs,1);
    cmatrix_expanded=zeros(npts,size(mindex,1));
    for j=1:npts
        candidate_nbrs=repmat(mindex(aux.coord_ptrs(j),:),nnbrs,1)+nbrvec_used; %use the multi=index of point j, and add all offsets
        within=ones(nnbrs,1); %determine if it is within the volume
        within(any(candidate_nbrs<0,2))=0;
        within(any(candidate_nbrs>=max_mindex,2))=0;
        nbrs=candidate_nbrs(within>0,:);
        nbr_coords=nbrs*mults+1;
        cmatrix_expanded(j,nbr_coords)=1;
    %    jnbrs=find(ismember(aux.coord_ptrs,nbr_coords)); %coordinate pointers that are also coordinates of neighbors
    %    cmatrix(j,jnbrs)=1;
    %    cmatrix(jnbrs,j)=1;
    end
    cmatrix=cmatrix_expanded(:,aux.coord_ptrs); %only keep columns corresponding to regions that exist
else
    cmatrix=cmatrix_full(aux.coord_ptrs,aux.coord_ptrs);
end
aux.cmatrix=cmatrix;
%
%determine connectivity
%
[isc,cmsg,niters]=isgraphc(cmatrix,1); %fast option
aux.msg=cmsg;
aux.niters=niters;
return

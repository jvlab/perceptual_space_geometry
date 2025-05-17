function [imgs,optsused,errs,metro]=btc_makemaps(method,opts,dict,auxopts,mtcopts)
% [imgs,optsused,errs,metro]=btc_makemaps(method,opts,dict,auxopts,mtcopts) makes one or more
% textures based on a "method" (from btc_augcoords)
%
% 26Jan15: add extension to multiple gray levels
% 01Oct15: add extension for FLFS
% 14Nov19: documentation fixes, no change in code
%
% method:  a structure returned by btc_augcoords or mtc_augcoords
%   method.name='Pickard' or 'DiagMRF' or 'NoPickTT' or 'NoPickBT' or 'FLFS'
% for Pickard: (Pickard, specify probabilities in 2x2 block, and propagate)
%   method.variant_num=[1 2]
%   method.p2x2: [2x2x2x2] or [ng x ng x ng x ng] array of probabilities,  Pickard conditions and normalization not checked
%   variant_num=1: Pickard condition involving dict.inpickard(1,:) holds; generate top to bottom and *L* to *R* in columns
%   variant_num=2: Pickard condition involving dict.inpickard(2,:) holds; generate top to bottom and *R* to *L* in columns
%     uses genmrfm for Pickard methods 1 and 2
% for DiagMRF: (Pickard on diagonal 2x2 sublattices)
%   method.diagp2x2: [2x2x2x2] array of probabilities to be used in a diagonal box, from btc_augcoords
%   method.p2x2_diag: [ng x ng x ng x ng] array playing same role as method.diagp2x2, from mtc_augcoords
%     Pickard conditions and normalization not checked
% for NoPickTT: (Non-Pickard method for thetas sharing an edge)
%   method.variant_num=[1 2 3 4] for number of rotations w.r.t. standard orientation
%      variant_num=1: standard orientation (theta(ACD) and (BCD).
%      Other variants get rotated by rot90(variant_num-1). 
%   Standard orientation is theta(ACD) and (BCD). Specify probabilities on
%   2x2 blocks, iterate, and then apply Metropolis donut algorithm to
%   reduce resulting next-nearest-neighbor probabilities)
% for NoPickBT: (Non-Pickard method for a theta and a beta coming from its corner)
%   method.variant_num=[1 2 3 4] for number of rotations w.r.t. standard orientation
%      variant_num=1: standard orientation (theta(ABC) and beta(AD).
%      Other variants get rotated by rot90(variant_num-1). 
%   Standard orientation is theta(ABC) and beta(AD). Specify probabilities on
%   a T-shaped region, iterate, and then apply Metropolis donut algorithm to
%   reduce resulting correlatins.  This induces
%   some of the opposite theta (BCD), proportional to theta(ABC)*beta(AD).
%   The T rule has to include a triple-correlation "across the top", which results
%   in a pairwise correlation at a horizontal offset of 3.
% for FLFS: (falling-leaves/falling-sticks method)
%   method.variant_num=1
%   method.flfs:  has information for the construction from mtc_augcoords; subfields are:
%        flfs: [1x1 struct], the flfs setup
%        flfs_configs: [1x14 struct], from flfs_enumerate
%        flfs_factors: [1x1 struct], from flfs_bc_factors
%        flfs_eivec: [14x1 double], the stable set of probabilities of configurations
%        flfs_weight: [1x1 struct], fields b, c, bb, and cc, indicating how much these correlations are diluted
%        flfs_p_hv:  the relative probability of horizontal stick (b); vertical has probability 1-flfs_p_hv
%        flfs_probs: {[ngxng double]  [ngxng double]}: 2-point probabilities for horiz (b) and vertical (c) sticks
% opts:
% opts.nmaps: number of maps to show
% opts.show: list maps to show (a subset of [1:nmaps], defaults to [])
% opts.area: size of map (length, width)
% opts.firstrow: holds the first (top) row, filled randomly if not present
%   if not empty:
%     size(firstrow,1)=area(2)
%     size(firstrow,2)=1 or nmaps
%   CAN ONLY BE SPECIFIED WITH SOME METHODS:  Pickard 1, Pickard 2
% opts.firstcol: holds the first (left) column, filled randomly if not present
%   if not empty:
%     size(firstcol,1)=area(1)
%     size(firstcol,2)=1 or nmaps
%   CAN ONLY BE SPECIFIED WITH SOME METHODS:  Pickard 1
%   firstcol(1,*) must = firstrow(1,*) if both are specified
% opts.lastcol:  holds the last (right) column, filled randomly if not present
%   if not empty:
%     size(lastcol,1)=area(1)
%     size(lasttcol,2)=1 or nmaps
%   CAN ONLY BE SPECIFIED WITH SOME METHODS:  Pickard 2
%   lastcol(1,*) must = firstrow(end,*) if both are specified
% opts.verbose:  verbose output of map making routine, if applicable
% opts.onaxis: 1 if on axis (i.e., at most one param value is nonzero); this bypasses Metropolis
%
% dict: dictionary of corrs; created with btc_define if not passed.
% auxopts: options structure set by BTC_AUXOPTS, auxopts.metro_opts contains Metropolis params
% mtcopts: options structure set by MTC_DEFOPTS, if passed
%
% output:
% imgs: the images requested, 0 to ng-1, size=[opts.area(1) opts.area(2) opts.nmaps]
% optsused: opts, with defaults filled in
% errs: empty if no error, otherwise an error message
% metro: structure with results from Metropolis algorithm
%
% number of gray levels determined from dimension of p2x2 in methods.
% if ng=2, methods may have been generated either by BTC_AUGCOORDS or MTC_AUGCOORDS,
%   so other fields of method are inspected for the non-Pickard methods
%
%   See also BTC_AUGCOORDS, BTC_DEFINE, BTC_TEST, BTC_SURV, GENMRFM,
%   GENMRFM_1D, MAPUBI, BTC_DIMRF, BTC_NOPICK, DONUT_METRO, MTC_FLFS_MAKEMAPS, RFP_BTC_MAKE, BTC_MAKEMAPS_SCALE.
%

method_list={'Pickard','DiagMRF','NoPickTT','NoPickBT','FLFS','DiagFLFS'};
method_nofirstrow={'DiagMRF','NoPickTT','NoPickBT','FLFS','DiagFLFS'};
method_nofirstcol={'DiagMRF','NoPickTT','NoPickBT','FLFS','DiagFLFS'};
method_nolastcol={'DiagMRF','NoPickTT','NoPickBT','FLFS','DiagFLFS'};
%
if (nargin<=4)
    mtcopts=[];
end
mtcopts=mtc_defopts(mtcopts);
if (nargin<=3)
    auxopts=btc_auxopts([],[]);
end
if ~isfield(auxopts,'metro_opts')
    auxopts_std=btc_auxopts;
    auxopts.metro_opts=auxopts_std.metro_opts;
    auxopts.metro_show=auxopts_std.metro_show;
end
%
if (nargin<=2)
    dict=btc_define();
end
if (nargin<=1)
    opts=[];
end
ng=size(method.p2x2,1); %number of gray levels
%
metro=[];
opts=filldefault(opts,'nmaps',1);
opts=filldefault(opts,'show',[]);
opts=filldefault(opts,'area',[256 256]);
opts=filldefault(opts,'firstrow',[]);
opts=filldefault(opts,'firstcol',[]);
opts=filldefault(opts,'lastcol',[]);
opts=filldefault(opts,'verbose',0);
opts=filldefault(opts,'onaxis',0); %1 if on axis (and bypasses need for Metropolis)
optsused=opts;
if length(opts.area)==1
    opts.area=repmat(opts.area,1,2);
end
imgs=zeros([opts.area opts.nmaps]);
errs=[];
if all(strcmp(method.name,method_list)==0)
    errs=sprintf('variant %s not yet implemented.',method.name);
    return
end
if ~isempty(opts.firstrow)
    if ~((size(opts.firstrow,1)==opts.area(2)) & (size(opts.firstrow,2)==1 | size(opts.firstrow,2)==opts.nmaps))
        errs=sprintf('first row is specified, but dimensions (%6.0f %6.0f) are incompatible with area (%6.0f %6.0f) and nmaps (%6.0f)',...
            size(opts.firstrow,1),size(opts.firstrow,2),opts.area(1),opts.area(2),opts.nmaps);
        return
    end
end
if ~isempty(opts.firstcol)
    if ~((size(opts.firstcol,1)==opts.area(1)) & (size(opts.firstcol,2)==1 | size(opts.firstcol,2)==opts.nmaps))
        errs=sprintf('first col is specified, but dimensions (%6.0f %6.0f) are incompatible with area (%6.0f %6.0f) and nmaps (%6.0f)',...
            size(opts.firstcol,1),size(opts.firstcol,2),opts.area(1),opts.area(2),opts.nmaps);
        return
    end
end
if ~isempty(opts.lastcol)
    if ~((size(opts.lastcol,1)==opts.area(1)) & (size(opts.lastcol,2)==1 | size(opts.lastcol,2)==opts.nmaps))
        errs=sprintf('last col is specified, but dimensions (%6.0f %6.0f) are incompatible with area (%6.0f %6.0f) and nmaps (%6.0f)',...
            size(opts.lastcol,1),size(opts.lastcol,2),opts.area(1),opts.area(2),opts.nmaps);
        return
    end
end
if ~isempty(opts.firstrow) & ~isempty(opts.firstcol)
    if ~all(opts.firstrow(1,:)==opts.firstcol(1,:))
        errs=sprintf('specified first row and first column not compatible.');
        return
    end
end
if ~isempty(opts.firstrow) & ~isempty(opts.lastcol)
    if ~all(opts.firstrow(end,:)==opts.lastcol(1,:))
        errs=sprintf('specified first row and last column not compatible.');
        return
    end
end
if ~isempty(opts.firstrow) & any(strcmp(method.name,method_nofirstrow)==1)
    errs=sprintf('cannot specify first row with %s',method.name);
    return
end
if ~isempty(opts.firstcol) & any(strcmp(method.name,method_nofirstcol)==1)
    errs=sprintf('cannot specify first col with %s',method.name);
    return
end
if ~isempty(opts.lastcol) & any(strcmp(method.name,method_nolastcol)==1)
    errs=sprintf('cannot specify last col with %s',method.name);
    return
end
%
gopts=[];
switch method.name
    case 'Pickard'
        if method.variant_num==1 & ~isempty(opts.lastcol)
            errs=sprintf('cannot specify last column with Pickard variant 1');
            return
        end
        if method.variant_num==2 & ~isempty(opts.firstcol)
            errs=sprintf('cannot specify first column with Pickard variant 2');
            return
        end
        for imap=1:opts.nmaps         
            gopts.show=0; %we show it here
            gopts.statsub=0;
            gopts.err_rept=1;
            rspec=[];
            if ~isempty(opts.firstrow)
                rspec=opts.firstrow(:,min(imap,size(opts.firstrow,2)));
            end
            cspec=[];
            if method.variant_num==1
                if ~isempty(opts.firstcol)
                    cspec=opts.firstcol(:,min(imap,size(opts.firstcol,2)));
                end
                p2x2=method.p2x2;
            end
            if method.variant_num==2
                if ~isempty(opts.lastcol)
                    cspec=opts.lastcol(:,min(imap,size(opts.lastcol,2)));
                end
                rspec=flipud(rspec);
                p2x2=permute(method.p2x2,[2 1 4 3]);
            end
            gopts.firstrow=rspec;
            gopts.firstcol=cspec;
            gopts.pblocks=p2x2;
            if (ng==2)
                [map,stats,ou,errmap]=genmrfm(gopts,opts.area);
            else
                gopts.noftps=1; %suppress calclation of Fourier transform coordinates
                [map,stats,ou]=genmrfmg(gopts,opts.area);
                errmap=[];
            end           
            if (method.variant_num==2)
                map=fliplr(map);
            end
            if isempty(errmap)
               imgs(:,:,imap)=map;
               if ismember(imap,opts.show)
                   btc_makemaps_show(map,imap,ng);
               end
            else
                errs=strvcat(errs,cat(2,sprintf(' in map %5.0f',imap),': ',errmap));
            end %errmap
        end %imap
    case 'DiagMRF'
        for imap=1:opts.nmaps
            gopts.show=0; %we show it here
            gopts.statsub=0;
            gopts.err_rept=1;
            gopts.firstrow=[];
            gopts.firstcol=[];
            gopts.lastcol=[];
            if isfield(method,'p2x2_diag') 
                gopts.pblocks=method.p2x2_diag; %this is from mtc_augcoords
            else
                gopts.pblocks=method.diagp2x2; %this is from btc_augcoords
            end
            if strcmp(method.variant_lab,'NWSE') | strcmp(method.variant_lab,'std') %explicitly NWSE or "standard"
                % This works for ng=2 or ng>2.
                [map,stats,ou,errmap]=btc_dimrf(gopts,opts.area);
            else
                %NESW:  run the NWSE in the opposite direction and then flip
                gopts_flip=gopts;
                gopts_flip.pblocks=permute(gopts.pblocks,[2 1 4 3]);
                opts_flip=opts;
                opts_flip.area=fliplr(opts.area);
                % This works for ng=2 or ng>2.
                % For ng=2, stats does not depend on flip.
                % For ng>2, stats are not computed
                [map_flip,stats,ou,errmap]=btc_dimrf(gopts_flip,opts_flip.area); 
                map=map_flip'; %transpose
            end
            if isempty(errmap)
               imgs(:,:,imap)=map;
               if ismember(imap,opts.show)
                   btc_makemaps_show(map,imap,ng);
               end
            else
                errs=strvcat(errs,cat(2,sprintf(' in map %5.0f',imap),': ',errmap));
            end %errmap
        end %imap
    case {'NoPickTT','NoPickBT'}
        gopts.nmaps=opts.nmaps; %all maps made at once
        gopts.show=0; %we show it here
        gopts.statsub=0;
        gopts.err_rept=1;
        gopts.firstrow=[];
        gopts.firstcol=[];
        gopts.lastcol=[];
        gopts.verbose=opts.verbose;
        gopts.rotvar=method.variant_num;
        gopts.metro_opts=auxopts.metro_opts;
        gopts.metro_show=auxopts.metro_show;
        gopts.recurtype=method.name;
        if isfield(auxopts,'burnin') 
            gopts.burnin=auxopts.burnin;
        end
        gopts.onaxis=opts.onaxis;
        %
        if strcmp(method.name,'NoPickTT')
            if (method.variant_num==1) %tu
                indices=[3 4];
            end
            if (method.variant_num==2) %tw
                indices=[4 2];
            end
            if (method.variant_num==3) %vw
                indices=[2 1];
            end
            if (method.variant_num==4) %uv
                indices=[1 3];
            end
            if (ng==2)
                corrs=getcorrs_p2x2(method.p2x2);
                gopts.mACD=corrs.theta(indices(1));
                gopts.mBCD=corrs.theta(indices(2));
                gopts.mABCD=corrs.alpha;
            end
        end
        if strcmp(method.name,'NoPickBT')
            if (method.variant_num==1) %vd
                indices=[1 3];
            end
            if (method.variant_num==4) %ue
                indices=[2 4];
            end
            if (method.variant_num==3) %td
                indices=[4 3];
            end
            if (method.variant_num==2) %we
                indices=[3 4];
            end
            if (ng==2)
                corrs=getcorrs_p2x2(method.p2x2);
                gopts.mABC=corrs.theta(indices(1));
                gopts.mAD=corrs.beta(indices(2));
            end
        end
        gopts.method=method;
        gopts.mtcopts=mtcopts;
        [imgs,stats,ou,errs,metro]=btc_nopick(gopts,opts.area);
        metro.stats=[];
        metro.errs=errs;
        metro.opts_used=ou;
        if ~isempty(errs)
            disp(errs)
        end
        for imap=1:opts.nmaps
            if ismember(imap,opts.show)
                btc_makemaps_show(imgs(:,:,imap),imap,ng);
            end
        end %imap
    case {'FLFS','DiagFLFS'}   
        gopts.show=0; %we show it here
        gopts.err_rept=1;
        gopts.onaxis=opts.onaxis;
        gopts.firstrow=[];
        gopts.firstcol=[];
        gopts.lastcol=[];
        gopts.verbose=opts.verbose;
        gopts.metro_opts=auxopts.metro_opts;
        gopts.metro_show=auxopts.metro_show;
        if strcmp(method.name,'DiagFLFS')
            %adjust number of maps to make
            gopts.nmaps=2*opts.nmaps;
            size_diag=ceil(opts.area(1)+opts.area(2)/2)+1; %adjust area, see btc_dimrf
            flfs_area=[size_diag size_diag];
        else
            gopts.nmaps=opts.nmaps; %all maps made at once
            flfs_area=opts.area;
        end
        [imgs_flfs,ou,errs,metro]=mtc_flfs_makemaps(setfield(method,'name','FLFS'),...
            gopts,flfs_area,dict,mtcopts);
        metro.stats=[];
        metro.errs=errs;
        metro.opts_used=ou;
        if strcmp(method.name,'DiagFLFS')
            %rotate maps back and interleave, using logic of btc_dimrf
            %but treating all maps in parallel
            imgs=zeros(opts.area(1),opts.area(2),opts.nmaps);
            for gx=1:size_diag
                gxc=gx-(size_diag+1)/2; %runs from -(size_diag-1)/2 to (size_diag-1)/2
                for gy=1:size_diag
                    gyc=gy-(size_diag+1)/2;
                    mx=floor(gxc-gyc+opts.area(1)/2);
                    my=floor(gyc+gxc+opts.area(2)/2);
                    if (mx>0) & (my>0)& (mx<=opts.area(1)) & (my<=opts.area(2))
                        imgs(mx,my,:)=imgs_flfs(gx,gy,[1:opts.nmaps]);
                    end
                    if (mx>=0) & (my>0)& (mx<opts.area(1)) & (my<=opts.area(2))
                        imgs(mx+1,my,:)=imgs_flfs(gx,gy,opts.nmaps+[1:opts.nmaps]);
                    end
                end
            end
        else
            imgs=imgs_flfs;
        end
        if ~isempty(errs)
            disp(errs)
        end
        for imap=1:opts.nmaps
            if ismember(imap,opts.show)
                btc_makemaps_show(imgs(:,:,imap),imap,ng);
            end
        end %imap
        if ~isempty(errs)
            disp(errs)
        end
        for imap=1:opts.nmaps
            if ismember(imap,opts.show)
                btc_makemaps_show(imgs(:,:,imap),imap,ng);
            end
        end %imap      
end %switch case
return

function btc_makemaps_show(map,imap,ng)
    tstring=sprintf('map %5.0f',imap);
    figure;
    set(gcf,'Name',tstring);
    imagesc(map,[0 ng-1]);
    axis equal;axis off;colormap('gray');
    title(tstring);
return


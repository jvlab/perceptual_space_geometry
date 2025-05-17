function [ra,ou]=btc_soidfg_anal(ru,opts,dict)
% [ra,ou]=btc_soidfg_anal(ru,opts,dict) analyzes models fit by btc_soidfg_demo
% analyzes:
%   number of parameters
%   fraction of variance explained
%  for the quantities below, quantiles from the surrogates for qfit are also calculated
%   eigenvalues
%   projection into subspace with guaranteed infinite threshold for fig=gnd
%   thresholds for each coordinate
%   lowest threhsold, and lowest threshold for fig=gnd
%
%   See ..\gr18\figgnd_modeling_notes.docx for background
%
%  ru: a single entry in the the cell array r from a workspace created by btc_soidfg_demo,
%   containing model setup, data, quadratic form fit, and basic statistics
%  opts: options
%    if_verbose, 1 (default: 0) show results here
%    fgdirs: a two-column matrix indicating figure and ground multipliers for computing model thresholds
%    quantiles: quantiles to show (defaults to [0.025 0.5 0.975])
%  dict: struture returned by btc_define; computed if omitted
%
%  ra: results
%  ou: options used
%
% BTC_SOIDFG_DEFINE, BTC_SOIDFG_DEMO, BTC_SOIDFG_MODEL, BTC_SOIDFG_MODEL_PROJ, BTC_SOIDFG_ANAL_DEMO,
% BTC_SOIDFG_FIT, BTC_SOIDFG_CALCSTATS, BTC_SOIDFG_ANAL_OF.
%
if (nargin<=2)
    dict=btc_define;
end
btc_n=length(dict.codel);
%
if (nargin<=1)
    opts=[];
end
opts=filldefault(opts,'if_verbose',0);
opts=filldefault(opts,'fgdirs',[1 0;0 1;1 -1;1 1]);
opts=filldefault(opts,'quantiles',[0.025 0.5 0.975]);
ou=opts;
%
ncoords=length(ru.model.btcpos);
nfgdirs=size(opts.fgdirs,1);
%
ra=struct();
nparams=ru.model.nparams;
ra.nparams.tot=nparams;
ra.nparams.self=length(strmatch('self',ru.model.param_type));
ra.nparams.cross=length(strmatch('cross',ru.model.param_type));
qforms=ru.model.qforms;
qforms_diag=zeros(size(qforms));
for iparam=1:nparams
    qforms_diag(:,:,iparam)=diag(diag(qforms(:,:,iparam)));
end
ra.nparams.diag=sum(any(qforms_diag~=0,[1 2]));
ra.nparams.off=ra.nparams.tot-ra.nparams.diag;
%
btcpos_f=ru.model.btcpos;
btcpos_g=ru.model.btcpos+btc_n;
ncoords=length(ru.model.btcpos);
btcpos_all=[btcpos_f btcpos_g];
qforms_reduced=qforms(btcpos_all,btcpos_all,:);
%check that the quadratic form components are zero where there are no coordinates
trivial=setdiff([1:2*btc_n],btcpos_all);
if any(qforms(trivial,trivial,:)~=0,[1 2 3])
    warning('nonzero elements of quadratic form in unexpected places.');
end
if (opts.if_verbose)
    disp(sprintf(' parameter counts: tot: %3.0f  self-terms %3.0f cross-terms %3.0f on-diag %3.0f off-diag %3.0f',...
        ra.nparams.tot,ra.nparams.self,ra.nparams.cross,ra.nparams.diag,ra.nparams.off))
    disp(sprintf(' nontrivial coordinates: %3.0f',ncoords));
    disp(sprintf(' fraction of variance explained by model: %7.3f',ru.results.stats.ss_ratio));
end
%analyses of the quadratic model
qfit_reduced=ru.results.qfit(btcpos_all,btcpos_all);
ra_struct=btc_soidfg_anal_do(qfit_reduced,dict,btcpos_f,opts);
fnames=fieldnames(ra_struct);
for ifn=1:length(fnames)
    ra.(fnames{ifn})=ra_struct.(fnames{ifn});
end
%do surrogates, concatenating along last dimension
if isfield(ru.results,'surrogates')
    ra.surrogates=struct();
    nsurr=size(ru.results.surrogates.qfit,3);
    for isurr=1:nsurr
        qfit_reduced=ru.results.surrogates.qfit(btcpos_all,btcpos_all,isurr);
        ra_struct=btc_soidfg_anal_do(qfit_reduced,dict,btcpos_f,opts);
        for ifn=1:length(fnames)
            val=ra_struct.(fnames{ifn});
            catdim=ndims(val)+1;
            if size(val,ndims(val))==1 & ndims(val)==2 %column vectors are concatenated on dim 2
                catdim=ndims(val);
            end
            if prod(size(val))==1 %scalars are concatenated on dim 1
                catdim=1;
            end
            if (isurr==1)
                ra.surrogates.(fnames{ifn})=val;
            else
                ra.surrogates.(fnames{ifn})=cat(catdim,ra.surrogates.(fnames{ifn}),val);
            end
        end
    end
end
if_quantiles=and(isfield(ru.results,'surrogates'),length(opts.quantiles)>0);
if (opts.if_verbose)
    %eigenvalues and eigenvectors
    disp(sprintf(' eigenvalue range for q: [%9.4f  %9.4f], q_sumbyquad: [%9.4f %9.4f]',ra.q_eivals([1 end]),ra.q_sumbyquad_eivals([1 end])));
    if if_quantiles
        %compute quantiles from surrogates
        q_eivals_sq=quantile(ra.surrogates.q_eivals([1 end],:),opts.quantiles,2);
        q_sumbyquad_eivals_sq=quantile(ra.surrogates.q_sumbyquad_eivals([1 end],:),opts.quantiles,2);
        for iq=1:length(opts.quantiles)
            disp(sprintf('        quantile %6.3f: [%9.4f  %9.4f]               [%9.4f %9.4f]',...
                opts.quantiles(iq),q_eivals_sq([1 end],iq),q_sumbyquad_eivals_sq([1 end],iq)));
        end %quantiles
    end
    %projection
    disp(sprintf(' best projection onto sum-by-quadrant matrix is with h=%7.4f, fraction of unexplained variance by projection is %7.5f',ra.q_hmin,ra.q_vunex));
    if if_quantiles
        %compute quantiles from surrogates
        q_hmin_sq=quantile(ra.surrogates.q_hmin(:),opts.quantiles,1);
        q_vunex_sq=quantile(ra.surrogates.q_vunex(:),opts.quantiles,1);
        for iq=1:length(opts.quantiles)
            disp(sprintf('        quantile %6.3f:                               %7.4f                                                    %7.5f',opts.quantiles(iq),q_hmin_sq(iq),q_vunex_sq(iq)));
        end %quantiles
    end
    %predicted thesholds
    fgstring=' model thresholds             ';
    for ifgdir=1:nfgdirs
        fgstring=cat(2,fgstring,sprintf('[f,g] mults: [%5.2f %5.2f]   ',opts.fgdirs(ifgdir,:)));
    end
    disp(fgstring);
    for icoord=1:ncoords
        fg_coord=sprintf('  coordinate   %s:      ',dict.codel(btcpos_f(icoord)));
        fg_data=sprintf('                    %9.4f',ra.thresh_btc(icoord,:));
        disp(cat(2,fg_coord,fg_data));
        if (if_quantiles)
            %compute quantiles from surrogates
            thresh_sq=quantile(ra.surrogates.thresh_btc(icoord,:,:),opts.quantiles,3);
            for iq=1:length(opts.quantiles)
                fg_data=sprintf('                    %9.4f',thresh_sq(1,:,iq));
                disp(sprintf('     quantile %6.3f:  %s',opts.quantiles(iq),fg_data));
            end %quantiles
        end
    end %icoord
    disp(sprintf(' thresholds from eigenvecs: lowest: %9.4f, lowest with fig=gnd: %9.4f, ratio: %9.4f',...
        ra.thresh_eig,ra.thresh_eig_sumbyquad,ra.thresh_eig_ratio));
    if (if_quantiles)
        %compute quantiles from surrogates
        thresh_eig_sq=quantile(ra.surrogates.thresh_eig,opts.quantiles,1);
        thresh_eig_sumbyquad_sq=quantile(ra.surrogates.thresh_eig_sumbyquad,opts.quantiles,1);
        thresh_eig_ratio_sq=quantile(ra.surrogates.thresh_eig_ratio,opts.quantiles,1);
        for iq=1:length(opts.quantiles)
            fg_data=sprintf('                    %9.4f',thresh_sq(1,:,iq));
            disp(sprintf('     quantile %6.3f:               %9.4f                       %9.4f         %9.4f',...
                opts.quantiles(iq),thresh_eig_sq(iq),thresh_eig_sumbyquad_sq(iq),thresh_eig_ratio_sq(iq)));
        end %quantiles
    end
end
return

function ra_struct=btc_soidfg_anal_do(qfr,dict,btcpos,opts)
%utility routine for btc_soidfg_anal
%qfr: quadratic form with unused coords removed
%dict: results of btc_define
%btcpos: positions of variables in qfr
%opts: options arguments from btc_soidfg_anal
%
ncoords=length(btcpos);
h_init=0.5;
%
ra_struct.q_frobnorm=sqrt(sum(qfr(:).^2));
ra_struct.q_rank=rank(qfr);
[v,d]=eig(qfr);
dvals=diag(d);
ra_struct.q_eivecs=fliplr(v); %eigenvectors in descending order
ra_struct.q_eivals=flipud(dvals);
%
%find best projection -- fminsearch
[h_min,fval,exitflag,output]=fminsearch(@(x) btc_soidfg_anal_of(x,qfr),h_init);
ra_struct.q_hmin=h_min;
[dsq,qh]=btc_soidfg_anal_of(h_min,qfr);
ra_struct.q_hproj=qh; %best projection
ra_struct.q_vunex=dsq/ra_struct.q_frobnorm.^2; %variance unexplained by projection
%
%analyze the summed-by-quadrant 
%
fpos=[1:ncoords];
gpos=ncoords+fpos;
q_sumbyquad=qfr(fpos,fpos)+qfr(fpos,gpos)+qfr(gpos,fpos)+qfr(gpos,gpos); %responsible for nonzero part of projection
ra_struct.q_sumbyquad_frobnorm=sqrt(sum(q_sumbyquad(:).^2));
ra_struct.q_sumbyquad_rank=rank(q_sumbyquad);
[v,d]=eig(q_sumbyquad);
dvals=diag(d);
ra_struct.q_sumbyquad_eivecs=fliplr(v); %eigenvectors in descending order
ra_struct.q_sumbyquad_eivals=flipud(dvals);
%
%thresholds in btc directions
nfgdirs=size(opts.fgdirs,1);
for icoord=1:ncoords
    for ifgdir=1:nfgdirs
        ra_struct.thresh_btc(:,ifgdir)=btc_soidfg_anal_thr(qfr,opts.fgdirs(ifgdir,:));
    end %ifgdir
end
%thresholds in eigen-directions
ra_struct.thresh_eig=1./sqrt(ra_struct.q_eivals(1));
maxeiv_sumbyquad=ra_struct.q_sumbyquad_eivals(1);
if maxeiv_sumbyquad<=0
    ra_struct.thresh_eig_sumbyquad=Inf;
    ra_struct.thresh_eig_ratio=Inf;
else
    ra_struct.thresh_eig_sumbyquad=1/sqrt(maxeiv_sumbyquad);
    ra_struct.thresh_eig_ratio=ra_struct.thresh_eig_sumbyquad/ra_struct.thresh_eig;
end
return

function thrs=btc_soidfg_anal_thr(qfr,fg)
%computes model threshold in a figure-ground direction specified by fg 
% (fg(1)=figure contrast multiplier, fg(2)=ground contrast multiplier),
% for each coordinate direction
ncoords=length(qfr)/2;
thrs=repmat(Inf,ncoords,1); %assume infinite threshold
for icoord=1:ncoords
    v=zeros(1,ncoords);
    v(icoord)=1;
    vfg=[v*fg(1) v*fg(2)];
    signal=vfg*qfr*vfg';
    if signal>0
        thrs(icoord)=1/sqrt(signal);
    end
end
return



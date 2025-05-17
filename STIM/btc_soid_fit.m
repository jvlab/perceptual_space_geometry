function [results,ou_fit,ou_dict]=btc_soid_fit(edirs_setup,edirs_data,opts_fit,dict)
% [results,ou_fit,ou_dict]=btc_soid_fit(edirs_setup,edirs_data,opts_fit,dict)
% fits a set of psychophysical data to an ellipsoid
%
%  modifications Nov 15 2012 to allow for fitting of models even if some thresholds are large
% (lines involving use_large)
%
%  April 7 2013:  noted that the conversion from cl to std omits the factor of 1/2
%   so that the surrogates reflect 98% confidence limits, not 95%, i.e.,
%   the standard deviation used to generate the surrogates is erroneously too large.
%  see eigmink_btc_notes.txt
%
%   structures with fields yx, and subfields that indicate the experimental
%   setup and psychophysical data
% edirs_setup.yx 
%   plane: 'yx', the btc coordinates used
%   opts, a specification of the options used to lay out the lane
%   and at least one of the fields specified in opts_fit.uvec_fields, typically
%     maxvecs: maximum vectors in each direction (size=[ndirs,2])
%     uvecs:   unit vectors in each direction (size=[ndirs,2])
%     thresh_vecs: threshold vectors in each direction (size=[ndirs,2])
% edirs_data.yx must have 
%     thresh_mags: threshold in each direction, as a column vector
%   and, optionally, has
%     thresh_mags_eblo: lower confidence limit on threshold
%     thresh_mags_ebhi: upper confidence limit on threshold
%
%   thresh_vecs, thresh_vecs_eblo, thresho_vecs_ebhi are not used, but 
%     represnt the in-plane vectors associated with the above magnitudes
% Note that for all *vecs* quantities, the following:
%    *vecs(:,1) corresponds to "x", the *second* character of the plane designator
%    *vecs(:,2) corresponds to "y", the *first* character of the plane designator
%
% opts_fit: options used for fitting
%   verbose:  1 for verbose terminal output
%   coords:  the coordinate set for fitting, e.g.,'gbcdetuvwa') (the order determines the order in qfit)
%   symfull: 1 for symmetric setup, 2 for full setup, 3 for symmetric setup without coordinate substitution of te for td
%     (defaults to 1, see BTC_PAIRSNEEDED)
%   consistency_tol: consistency-check tolerance, defaults to 10^-5
%   consistency_chk: fields of edirs.setup to check for consistent directions
%     (defaults to strvcat('uvecs','maxvecs','thresh_vecs'); at least one
%     of these must be present
%   nsurr:  number of surrogates to analyze (surrogates generated based on supplied error bars)
%     defaults to 0
%   ebpval:  error-bar p-value, defaults to 0.05, used to create surrogates
%   minvec: minimum absolute value for a vector to be used as a regressor
%   maxvec: maximum absolute value for a vector to be used as a regressor
%   ifaug: 1 to augment the coords (as is done to create the textures), 0 to not augment
%   
% dict: btc options, used for augmenting coordinates
%
% results: results of the fitting, chiefly
%    results.qfit: the fitted quadratic form, as a symmetric matrix.
%        aug_vecs'*results.qfit*aug_vecs = 1 as nearly as possible
%    results.surrogates: corresponding surrogates, only contains fields that differ from real data
% ou_fit: opts_fit after defaults
% dict: dict after defaults
%
%  Note that edirs_setup and edirs_data can be the same,
%  and that test examples can be made by BTC_SOID_TEST, 
%  with edirs_setup=edirs, edirs_data=edirs_noise_thr
%
%   See also:  BTC_EDIRS, BTC_SOID_PLOT, BTC_AUGCOORDS, BTC_DEFINE, BTC_SOID_TEST, BTC_PAIRSNEEDED, 
%    BTC_PAIRSYMEXTEND, REGRESS, BTC_SOIDFG_FIT.
% 
if (nargin<=2)
    opts_fit=[];
end
opts_fit=filldefault(opts_fit,'verbose',0);
opts_fit=filldefault(opts_fit,'coords','gbcdetuvwa');
opts_fit=filldefault(opts_fit,'symfull',1);
opts_fit=filldefault(opts_fit,'consistency_tol',10^-5);
opts_fit=filldefault(opts_fit,'consistency_chk',strvcat('uvecs','maxvecs','thresh_vecs'));
opts_fit=filldefault(opts_fit,'nsurr',0);
opts_fit=filldefault(opts_fit,'ebpval',0.05);
opts_fit=filldefault(opts_fit,'minvec',0.0001);
opts_fit=filldefault(opts_fit,'maxvec',10);
opts_fit=filldefault(opts_fit,'ifaug',1);
opts_fit=filldefault(opts_fit,'use_large',0); %set to 1 to allow use of out-of-range values
%
if (nargin<=3)
    dict=[];
end
dict=btc_define(dict);
ou_fit=opts_fit;
ou_dict=dict;
qfit=zeros(length(opts_fit.coords));
results=[];
results.errs=[];
%
nsurr_todo=opts_fit.nsurr;
%
% do consistency checks
%
[planes_sym,planes_full,planes_sym_unsub]=btc_pairsneeded(opts_fit.coords,dict);
if (opts_fit.symfull==1) planes=planes_sym; end
if (opts_fit.symfull==2) planes=planes_full; end
if (opts_fit.symfull==3) planes=planes_sym_unsub; end
nplanes=size(planes,1);
results.nplanes=nplanes;
results.planes=planes;
allok=1;
%
% consistency check, and set up unit vectors (soid_fit)
% 
if (opts_fit.verbose>0) disp('checking consistency.');end
soid_fit=[];
for iplane=1:nplanes
    ps=planes(iplane,:);
    if ~isfield(edirs_setup,ps)
        errmsg=sprintf('plane %s: Not present in edirs_setup.',ps);
        results.errs=strvcat(results.errs,errmsg);
        if (opts_fit.verbose>0) disp(results.errs(end,:)); end
    end
    if ~isfield(edirs_data,ps)
        errmsg=sprintf('plane %s: Not present in edirs_data.',ps);
        results.errs=strvcat(results.errs,errmsg);
        if (opts_fit.verbose>0) disp(results.errs(end,:)); end
    end
    if isfield(edirs_setup,ps) & isfield(edirs_data,ps) 
        es=getfield(edirs_setup,ps);
        ed=getfield(edirs_data,ps);
        % now check that  maxvecs, uvecs, or thresh_vecs are present,
        % that they have the same number of elements, and that the directions match, within tol
        reference=[];
        ifok=1;
        chkmsg=[];
        for ifc=1:size(opts_fit.consistency_chk,1)
            fn=deblank(opts_fit.consistency_chk(ifc,:));
            if isfield(es,fn)
                uvecs=getfield(es,fn);
                uvecs=uvecs./repmat(sqrt(sum(uvecs.^2,2)),1,2); % make each row a unit vector
                if (isempty(reference)) %first time that a field with unit-vector data is encountered
                    chkmsg=sprintf('plane %2s: reference field is %s',ps,fn);
                    reference=fn;
                    uvecs=getfield(es,fn);
                    uvecs=uvecs./repmat(sqrt(sum(uvecs.^2,2)),1,2);
                    soid_fit=setfield(soid_fit,ps,setfield([],'uvecs',uvecs));
                else
                    uvecs_ref=getfield(getfield(soid_fit,ps),'uvecs');
                    if all(size(uvecs)==size(uvecs_ref))
                        if (max(max(abs(uvecs-uvecs_ref))))<= opts_fit.consistency_tol
                            chkmsg=cat(2,chkmsg,sprintf(', consistent with %s',fn));
                        else
                            ifok=0;
                            chkmsg=sprintf('plane %s: value mismatch of %s and %s.',ps,reference,fn);
                        end
                    else
                        chkmsg=sprintf('plane %s: size mismatch of %s and %s.',ps,reference,fn);
                        ifok=0;
                    end
                end
            end
        end
        if (isempty(reference))
            chkmsg=sprintf('plane %s: missing unit vector data.',ps);
            ifok=0;
        end
        if (ifok==0)
            results.errs=strvcat(results.errs,chkmsg);
            if (opts_fit.verbose>0) disp(results.errs(end,:)); end
            allok=0;
        end
        if (ifok==1 & opts_fit.verbose>0) disp(chkmsg); end
        if opts_fit.nsurr>0
            if ~isfield(ed,'thresh_mags_eblo') | ~isfield(ed,'thresh_mags_ebhi')
                ebmsg=sprintf('plane %s: surrogates requested but error thresh_mags_[eblo|ebhi] not supplied',ps);
                results.errs=strvcat(results.errs,ebmsg);
                if (opts_fit.verbose>0) disp(results.errs(end,:)); end
                nsurr_todo=0;
            end
        end
    end %edirs setup and data exist
end %iplane
if (allok==0)
    disp(results.errs);
    return
end
if (nsurr_todo>0)
    disp(sprintf('%3.0f surrogates requested, based on supplied error bars for p=%6.3f',opts_fit.nsurr,opts_fit.ebpval));
elseif (opts_fit.nsurr>0)
    disp(sprintf('%3.0f surrogates requested, but no error bars supplied',opts_fit.nsurr));
else
    disp(sprintf('%3.0f surrogates requested',opts_fit.nsurr));
end
%
% set up the axes (pure quadratic terms for regression)
%
qsetup=[];
[axes_sym,axes_full]=btc_axesneeded(opts_fit.coords,dict);
if ismember(opts_fit.symfull,[1 3])
    axes_list=axes_sym;
else
    axes_list=axes_full;
end
results.axes=axes_list;
%
% find how they enter into the quadratic form
%
for iaxis=1:length(axes_list)
    as=axes_list(iaxis,:);
    [vars,qones,qone]=btc_pairsymextend(as,opts_fit.coords,dict);
    setup=[];
    if ismember(opts_fit.symfull,[1 3])
        setup.vars=vars;
        setup.qones=qones;
    else
        setup.vars=as;
        setup.qones=qone;
    end
    qsetup=setfield(qsetup,as,setup);
end
%
% set up how the cross-terms enter into the quadratic form
%
for iplane=1:nplanes
    ps=planes(iplane,:);
    [vars,qones,qone]=btc_pairsymextend(ps,opts_fit.coords,dict);
    setup=[];
    if ismember(opts_fit.symfull,[1 3])
        setup.vars=vars;
        setup.qones=qones;
    else
        setup.vars=ps;
        setup.qones=qone;
    end
    qsetup=setfield(qsetup,ps,setup);
end
results.qsetup=qsetup;
results.soid_fit=soid_fit;
%
%find the subset of btc coordinates that are being examined
%
codel_all=getfield(btc_define([]),'codel');
for ic=1:length(opts_fit.coords)
    results.which_axes(ic)=find(opts_fit.coords(ic)==codel_all);
end
%
surrogates=[];
if nsurr_todo>0
    eb_to_std=2*norminv(1-opts_fit.ebpval); %this converts an error bar range into a std dev
end
for isurr=0:nsurr_todo
    %
    % create regression matrix
    %
    if (isurr==0)
        rstring='';
        edirs_data_use=edirs_data;
    else
        rstring=sprintf(', surrogate %4.0f',isurr);
        edirs_data_use=[]; %for surrogate data
        %create surrogate data, as a Gaussian centered on the mean of the error bar
        %with an stdv appropriate for ebpval
        for iplane=1:nplanes
            ps=planes(iplane,:);
            ebars=[edirs_data.(ps).thresh_mags_ebhi-edirs_data.(ps).thresh_mags_eblo];
            stdvs=abs(ebars)./eb_to_std;
            means=[edirs_data.(ps).thresh_mags_ebhi+edirs_data.(ps).thresh_mags_eblo]/2;
            threshs=edirs_data.(ps).thresh_mags;
            edirs_data_use.(ps).thresh_mags=means+stdvs.*randn(size(stdvs,1),1);
        end
        surrogates.edirs_data_use{isurr}=edirs_data_use;
    end
    if (opts_fit.verbose>0 & isurr==0) disp('creating regression matrix.');end
    [x,regressors,augvecs_all,augvecs_sel]=btc_soid_fit_makex(results,edirs_data_use,opts_fit,dict,isurr);
    if (opts_fit.verbose>0 &isurr==0) disp('regressors'); disp(regressors'); end
    if (isurr==0)
        results.x=x;
        results.regressors=regressors;
        results.augvecs_all=augvecs_all;
        results.augvecs_sel=augvecs_sel;
        %
        % condition numbers (only for real data)
        %
        condn_raw=cond(results.x);
        av=results.x;
        avc=av./repmat(sqrt(diag(av'*av))',size(av,1),1);
        condn_adjcol=cond(avc);
        avr=av./repmat(sqrt(diag(av*av')),1,size(av,2));
        condn_adjrow=cond(avr);
        if (opts_fit.verbose>0)
            disp(sprintf('regression condition numbers: %8.3f (raw) %8.3f (column adjusted) %8.3f (row adjusted)',...
                condn_raw,condn_adjcol,condn_adjrow));
        end
        results.condn_raw=condn_raw;
        results.condn_adjcol=condn_adjcol;
        results.condn_adjrow=condn_adjrow;
    else
        surrogates.x(:,:,isurr)=x;
        % surrogates.regressors(:,:,isurr)=regressors; not needed since these are all the same
        surrogates.augvecs_all(:,:,isurr)=augvecs_all;
        surrogates.augvecs_sel(:,:,isurr)=augvecs_sel;
    end
    %
    % do the regression
    %
    if (isurr==1);hwait=waitbar(0,sprintf('Doing %2.0f surrogate regressions',nsurr_todo));end
    if (isurr>0);waitbar(isurr/nsurr_todo,hwait);end
    if (opts_fit.verbose>0) disp(sprintf('doing regression%s',rstring));end
    %
    b=regress(ones(size(x,1),1),x);
    %
    if (isurr==nsurr_todo) & (isurr>0);close(hwait);end
    %
    % reconstitute the quadratic form from the regression matrix
    %
    qfit=zeros(length(opts_fit.coords));
    for ic=1:length(results.regressors)
        qfit=qfit+b(ic)*getfield(getfield(results.qsetup,results.regressors{ic}),'qones');
    end
    if (isurr==0)
        results.b=b;
        results.qfit=qfit;
    else
        surrogates.b(:,isurr)=b;
        surrogates.qfit(:,:,isurr)=qfit;
    end
%
end %isurr
results.surrogates=surrogates;
%
return

function [x,regressors,augvecs_all,augvecs_sel]=btc_soid_fit_makex(results,edirs_data,opts_fit,dict,isurr)
% now, for each threshold in edirs.data, find the augmented coordinate
% and set up the regression matrix x
%
% each row of x corresponds to a threshold measurement
% each column of x corresponds to a degree of freedom in the quadratic form, i.e., a field of qsetup
%
planes=results.planes;
qsetup=results.qsetup;
which_axes=results.which_axes;
soid_fit=results.soid_fit;
%
nplanes=size(planes,1);
regressors=fieldnames(qsetup);
x=zeros(0,length(regressors));
nx=0;
coords=opts_fit.coords;
for iplane=1:nplanes
    ps=planes(iplane,:);
    ed=getfield(edirs_data,ps);
    thresh_mags=getfield(ed,'thresh_mags');
    uvecs=getfield(getfield(soid_fit,ps),'uvecs');
    thresh_vecs=repmat(thresh_mags,1,2).*uvecs;
    for idir=1:size(thresh_vecs,1)
        spec=[];
        for ix=1:2
            %the "3-ix" is so that the SECOND coordinate of vecs_inplane 
            %is the FIRST btc coordinate (in planes)
            spec=setfield(spec,ps(ix),thresh_vecs(idir,3-ix));
        end
        if (opts_fit.ifaug)
            augcoords=btc_augcoords(spec,dict);
            av=augcoords.method{1}.vec;
        else %just use the coordinates without augmentation
            av=btc_letcode2vec(spec,dict);
            av(find(isnan(av)))=0;
        end
        avmag=sqrt(sum(av.*av));
        if (avmag>opts_fit.minvec) & (avmag<opts_fit.maxvec)
            inrange=1;
        else
            inrange=0;
        end
        if (inrange==1) | (opts_fit.use_large==1)
            nx=nx+1;
            augvecs_all(nx,:)=av; %full augmented vector
            augvecs_sel(nx,:)=augvecs_all(nx,which_axes);
            vec_outer=augvecs_sel(nx,:)'*augvecs_sel(nx,:); %outer product of augvec
            for ir=1:length(regressors)
                qones=getfield(getfield(qsetup,regressors{ir}),'qones');
                x(nx,ir)=sum(sum(vec_outer.*qones));
            end
        end
        if (inrange==0)
            disp(sprintf('warning: surrogate %4.0f: threshold magnitude is out of range',isurr));
            disp(spec)
            disp(augvecs_all(nx,:))
            if (opts_fit.use_large==1)
                disp(' used anywaybecause use_large is set to 1');
            end
        end
    end
end %iplane
return

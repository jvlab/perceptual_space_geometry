function [vecs,augcoords_method,strats,pwrs,aug_opts_used]=btc_augcoords_mults(spec,mults,dict,aug_opts)
% [vecs,augcoords_method,strats,pwrs,aug_opts_used]=btc_augcoords_mults(spec,mults,dict,aug_opts) determines
%    the coordinates of a mixture of textures corresponding to specs
%
%  If an augmentation exists, but btc coords cause probabilities to fall
%  out of the [0 1] range, then this attempts to extrapolate the augmented
%  coordinates by determining the power-law dependence
%
% spec: a structure of up to 10 fields, chosen from the code letters gbcdetuvwa,
%    indicating the pairwise coordinates of the mixtures to be combined
% mults: a column list of multipliers for spec, if empty, assumed to be 1
%    if extrapolation is needed then mults must be positive
% dict:  dictionary of binary correlation names, generated by btc_define
% aug_opts:  a field of options
%   aug_opts.ifstd: 1 if dict is the standard one (saves recomputing), default
%   aug_opts.nocheck: 1 not call getcorrs_p2x2 in btc_augcoords
%      There is no reason to specify nocheck=0, since probability range is checked here.
%   aug_opts.extrap_mults: two multipliers used for extrapolation, default to [0.125 0.25]
%   aug_opts.extrap_dilate_factor:  factor to increase extrapolation multipliers by on each iteration, defaults to sqrt(2)
%   aug_opts.variant_use: method variant to use, one of 'first','last', defaults to 'first'
%   aug_opts.tol: tolerance for checking probabilty >0 and  normalization,
%      and other checks vs. 0.  defaults to 10^-6, 
%
% vecs: array of size [length(mults),10], a set of augmencted coordinages 
% augcoords_method:  results returned by btc_augcoords for spec (with multiplier 1)
%    augcoords.method{:}: a cell array of methods used
% strats: vector of length nmults, strategy used.  Non-negative return indicates that augmented vectors are NaN's
%   -2->extrapolation fails becaaue mults is negative
%   -1->no method available in btc_augcoords, or, if extrapolated, with extrap_mults
%    0->btc_augcoords
%    1->btc_augcoords extrapolated
% pwrs: power laws determined for extrapolation (NaN if not used)
% aug_opts_used: aug_opts with defaults
%
%    See also: BTC_DEFINE, BTC_AUGCOORDS, BTC_AUGCOORDS_MIX, BTC_SOIDFG_DEFINE.
%
if (nargin<3)
    dict=btc_define;
    aug_opts.ifstd=1;
end
if (nargin<4)
    aug_opts=[];
end
if (isempty(dict))
    dict=btc_define;
    aug_opts.ifstd=1;
end
aug_opts=filldefault(aug_opts,'ifstd',0);
aug_opts=filldefault(aug_opts,'nocheck',1);
aug_opts=filldefault(aug_opts,'extrap_mults',[0.125 0.25]);
aug_opts=filldefault(aug_opts,'extrap_iter_factor',sqrt(2));
aug_opts=filldefault(aug_opts,'variant_use','first');
aug_opts=filldefault(aug_opts,'tol',10^-6);
%
if isempty(mults)
    mults=1;
end
btc_n=length(dict.codel);
mults=mults(:);
nmults=length(mults);
vecs=NaN(nmults,btc_n);
pwrs=NaN(1,btc_n);
augcoords_method=[];
aug_opts_used=aug_opts;
augcoords=btc_augcoords(spec,dict,aug_opts);
strats=NaN(nmults,1);
if isempty(augcoords.method)
    strats=repmat(-1,nmults,1);
    return
end
%determine method for spec
switch aug_opts.variant_use
    case 'first'
        m_use=1;
    case 'last'
        m_use=length(augcoords.method);
end
augcoords_method=augcoords.method(m_use);
%do augmentation for spec*mults
need_extrap=[];
spec_fields=fieldnames(spec);
maxabs=0;
for ifield=1:length(spec_fields)
    maxabs=max([maxabs,abs(spec.(spec_fields{ifield}))]);
end
if maxabs<aug_opts.tol
    vecs=zeros(nmults,btc_n);
    strats=zeros(nmults,1);
    return
end
for imult=1:nmults
    if (abs(imult)<aug_opts.tol)
        vecs(imult,:)=0;
        strats(imult)=1;
    else
        spec_mult=mult_spec(spec,mults(imult),spec_fields); %scale spec
        augcoords=btc_augcoords(spec_mult,dict,setfield(aug_opts,'nocheck',1));
        method=augcoords.method{m_use};
        if all(method.p2x2(:)>-aug_opts.tol) & double(abs(sum(method.p2x2(:))-1)<aug_opts.tol)
            vecs(imult,:)=method.vec;
            strats(imult)=0;
        else
            need_extrap=[need_extrap,imult];
            strats(imult)=1;
        end
    end       
end
need_extrap_pos=intersect(need_extrap,find(mults>=aug_opts.tol));
need_extrap_neg=setdiff(need_extrap,need_extrap_pos);
strats(need_extrap_neg)=-2;
if ~isempty(need_extrap_pos)
    %determine vectors for extrapolation, and if extrapolation is
    %successful move twice as far from origin and try again
    ve=NaN(2,btc_n);
    dilfact=1;
    if_trydouble=1;
    while (if_trydouble)
        ve_try=NaN(2,btc_n);
        for ive=1:2
            spec_mult=mult_spec(spec,dilfact*aug_opts.extrap_mults(ive)/maxabs,spec_fields);
            augcoords=btc_augcoords(spec_mult,dict,setfield(aug_opts,'nocheck',1));
            method=augcoords.method{m_use};
            if all(method.p2x2(:)>-aug_opts.tol) & double(abs(sum(method.p2x2(:))-1)<aug_opts.tol)
                ve_try(ive,:)=method.vec;
            end
        end
        if any(isnan(ve_try(:)))
            if_trydouble=0;
        else
            dilfact_use=dilfact; %save this attempt
            ve=ve_try;
            dilfact=dilfact*aug_opts.extrap_iter_factor; %try twice as far from origin
            aug_opts_used.extrap_used=dilfact*aug_opts.extrap_mults;
        end
    end %if_trydouble
    if any(isnan(ve(:)))
        strats(need_extrap)=-1; %cannot extapolate
        return
    end
    %extrapolate by power law
    for ilet=1:btc_n
        if all(abs(ve(:,ilet))>aug_opts.tol)
            pwr=log(ve(2,ilet)/ve(1,ilet))/log(aug_opts.extrap_mults(2)/aug_opts.extrap_mults(1));
            vecs(need_extrap_pos,ilet)=ve(1,ilet)*(mults(need_extrap_pos)/(aug_opts.extrap_mults(1)*dilfact_use/maxabs)).^pwr;
            pwrs(ilet)=pwr;
        else
            vecs(need_extrap_pos,ilet)=0;
        end
    end
end
return


function spec_mult=mult_spec(spec,mult,spec_fields)
%utility to multiply a spec by a vector
spec_mult=[];
for ifield=1:length(spec_fields)
    ilet=spec_fields{ifield};
    spec_mult.(ilet)=spec.(ilet)*mult;
end
return

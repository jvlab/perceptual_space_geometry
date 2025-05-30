function augcoords=btc_augcoords(spec,dict,aug_opts);
% augcoords=btc_augcoords(spec,dict,aug_opts) attempts to augment a subset of 
% 10 coordinates into a full set  of coords so that a 2x2 maxent texture can be generated
%
% spec: a structure of up to 10 fields, chosen from the code letters gbcdetuvwa
% dict:  dictionary of binary correlation names, generated by btc_define;
%    btc_define called if this field is empty or absent
% aug_opts:  a field of options
%   aug_opts.ifstd: 1 if dict is the standard one (saves recomputing)
% aug_opts.nocheck: %don't call getcorrs_p2x2
%   (so normalization is not checked, but assumed valid)
%
% augcoords.aug_opts: options used
% augcoords.dict: dictionary used
% augcoords.method{:}: a cell array of methods used
%    augcoords.method{:}.name={Pickard | DiagMRF | NoPickTT | NoPickBT}
%      Pickard:   Markov propagation in 2x2 block (most rules)
%      DiagMRF:   Markov propagation in checkerboard-interleaved diagonal sublattices (both diagonal betas)
%      NoPickTT: Markov propagation in 2x2 block, followed by Metropolis
%        for thetas overlapping on an edge, since Markov propagation induces next-nearest-neighbor correls
%      NoPickBT: Markov propagation in a T, followed by Metropolis for theta and beta coming from angle,
%        T-propagation induces beta*theta of the opposite theta, and theta*beta^2 across the top of the T
%    Pickard variant 1 must be generated from NW(A) to SE(D); flankers of A and D pixels are conditionally indep of A and D
%    Pickard variant 2 must be generated from NE(B) to SW C); flankers of B and C pixels are conditionally indep of B and C.
% Note that genmrfmg only generates from upper NW (upper left) to SE (lower right)
%
% Policies for creating the coordinates:
%  0.  Consistency:  getp2x2_atg must be consistent with getp2x2_ag, getp2x2_tg, getp2x2_at.
%  1.  The unspecified coordinates should approach zero as o(specified coordinates), and
%      ideally as O(specified coordinates^2)
%  2.  Symmetry should be maintained.  Nonzero even-order correls should not induce nonzero odd=-order correls
%  3.  High-order coordinates should not influence low-order coordinates.
%      For example, specifying nonzero thetas should not induce a nonzero beta
%      or gamma.  EXCEPTION is that thetas sharing an edge may induce a
%      diagonal beta, since this is done via a Metropolis algorithm.
%      Within rank:
%      * h and v betas DO induce the diagonal beta as their product, necessary for Pickard.
%      * h and diag betas do not induce other nonzero betas (but do induce horizontal correls at larger distances)
%      * thetas do not induce other thetas
%      * theta and a spanning beta do not induce other betas
%      * theta and a beta emanating from its corner: depends on details of Metropolis implementation
%  4.  Low-order coordinates should multiply to form starting point for high-order coordinates � 
%      betas should be gamma^2 if not specified; thetas should be gamma^3 if beta and thetas not specified;
%      alphas should be gamma^4 if beta, theta, alpha not specified
%  5.  Low-order coordinates may induce high-order coordinates via maximum entropy.
%      Often this means multiplying them to form starting point for high-order coordinates � betas should be gamma^2
%      if not specified; thetas should be gamma^3 if beta and thetas not specified; alphas should be gamma^4 if
%      beta, theta, alpha not specified.  But for other choices, the alpha
%      is found by a polynomial equation (see btc_meabb, btc_meabt).
%      For a pair of horizontal and vertical betas, alpha (equal to the product) is induced.
%      For a pair of diagonal betas, a 4-element correlation is induced on a larger scale rotated square.
%      For one cardinal and one diagonal beta, a 4-element correlation is induced on a skewed diamond.
%
%  19Mar20:  bug fixed in call to expthame
%
%    See also: BTC_DEFINE, BTC_TEST, GETP2X2_ATG, GETCORRS_P2X2, GETP2X2_CORRS, BTC_MEABB, BTC_MEABT, BTC_PAIRVERIFY.
%
if (nargin<2)
    dict=btc_define;
end
if (nargin<3)
    aug_opts=[];
end
if (isempty(dict))
    dict=btc_define;
    aug_opts.ifstd=1;
end
aug_opts=filldefault(aug_opts,'tol',10^-6);
aug_opts=filldefault(aug_opts,'ifstd',0);
aug_opts=filldefault(aug_opts,'nocheck',0);
if (aug_opts.ifstd==0)
    dict_std=btc_define; %for when we need a standard dictionary
else
    dict_std=dict;
end
if isempty(spec) %added 13Aug18 to prevent errors if spec=[] rather than struct()
    spec=struct();
end
spec=btc_letcode_strip(spec,dict_std);

exptname=btc_exptname(char(fieldnames(spec))',dict); %' added 19Mar20 so that this works with, say, fields g, b, and c of spec
coorkinds=btc_coorkinds(char(fieldnames(spec)),dict);
augcoords=[];
method=[];
%
name1=dict.ordernames{1};
name2=dict.ordernames{2};
name2hv=cat(2,dict.ordernames{2},'_hv');
name2di=cat(2,dict.ordernames{2},'_diag');
name3=dict.ordernames{3};
name4=dict.ordernames{4};
%
n1=length(getfield(coorkinds,name1));
n2hv=length(getfield(coorkinds,name2hv));
n2di=length(getfield(coorkinds,name2di));
n3=length(getfield(coorkinds,name3));
n4=length(getfield(coorkinds,name4));
n=[n1 n2hv+n2di n3 n4];
%
nspec=length(exptname);
% do many cases here
method=[];
corrs_z=btc_vec2corrs(zeros(1,length(dict.codel)),dict);
letcode_z=btc_vec2letcode(btc_corrs2vec(corrs_z,dict),dict);
%
for q=1:length(exptname)
    coordnum(q)=find(dict.codel==exptname(q));
    val(q)=getfield(spec,exptname(q));
    posit(q)=dict.posit(coordnum(q));
end
%
if (nspec==0) %null case
    method=btc_augcoords_fillmeth(corrs_z,[1 2]); %both Pickard branches
end
if (nspec==1) %one-variable cases
    [corrs,vnums]=btc_augcoords_1var(coorkinds,spec,exptname,dict);
    method=btc_augcoords_fillmeth(corrs,vnums);
end
if (nspec==2) %two-variable cases
    if (n4==1)
        %if alpha is present, find the corrs when alpha=0, and then set alpha
        spec_reduced=rmfield(spec,'a');
        exptname_reduced=btc_exptname(char(fieldnames(spec_reduced)),dict);
        coorkinds_reduced=btc_coorkinds(char(fieldnames(spec_reduced)),dict);
        [corrs,vnums]=btc_augcoords_1var(coorkinds_reduced,spec_reduced,exptname_reduced,dict);
        corrs=setfield(corrs,name4,spec.a);
        method=btc_augcoords_fillmeth(corrs,vnums);
    else
        %two variables, both not alpha
        % for example, theta, gamma:  betas are gamma^2, other thetas are gamma^3 and alpha=gamma*theta
        [corrs,vnums,mname]=btc_augcoords_2var(coorkinds,spec,exptname,dict);
        method=btc_augcoords_fillmeth(corrs,vnums,mname);
    end
end
if (nspec==3) %three-variable cases
    if (n4==1)
        %if alpha is present, find the corrs when alpha=0, and then set alpha
        spec_reduced=rmfield(spec,'a');
        exptname_reduced=btc_exptname(char(fieldnames(spec_reduced)),dict);
        coorkinds_reduced=btc_coorkinds(char(fieldnames(spec_reduced)),dict);
        [corrs,vnums,mname]=btc_augcoords_2var(coorkinds_reduced,spec_reduced,exptname_reduced,dict);
        if strcmp(mname,'Pickard') | strcmp(mname,'NoPickTT');
            corrs=setfield(corrs,name4,spec.a);
            method=btc_augcoords_fillmeth(corrs,vnums,mname);
        else
            method=btc_augcoords_fillmeth(corrs,[],[]);
        end
    else
        %three variables, all not alpha
        corrs=corrs_z;
        %no algorithms at present
        method=btc_augcoords_fillmeth(corrs,[],[]);
    end
end
%
% Common completion
%
% If any corrs's exist but not p2x2, make the p2x2.
for im=1:length(method)
    if (isfield(method{im},'corrs') & ~isfield(method{im},'p2x2'))
        if (aug_opts.ifstd==1)
            p2x2=getp2x2_corrs(method{im}.corrs); %standard order-names used
        else
            p2x2=getp2x2_corrs(btc_vec2corrs(btc_corrs2vec(method{im}.corrs,dict),dict_std)); %in case nonstandard order-names used
        end
        method{im}.p2x2=p2x2;
    end
    if (isfield(method{im},'corrs') & ~isfield(method{im},'vec'))
        method{im}.vec=btc_corrs2vec(method{im}.corrs,dict);
    end
    if (isfield(method{im},'corrs') & ~isfield(method{im},'letcode'))
        method{im}.letcode=btc_vec2letcode(method{im}.vec,dict);
    end
end
% also check p-values, normalization, and, if relevant, Pickard conditions
for im=1:length(method)
    if (aug_opts.nocheck==0)
        gc=getcorrs_p2x2(method{im}.p2x2,1);
    else
        gc=[];
        gc.entropy=0;  %entropy as 2x2 MRF
        gc.cig_conds=zeros(2,4);
        gc.norm=1;
        gc.ok_probs=1;
        gc.ok_norm=1;
    end
    method{im}.entropy_mrf22=gc.entropy;  %entropy as 2x2 MRF
    method{im}.cig_conds=gc.cig_conds;
    method{im}.norm=gc.norm;
    method{im}.ok_probs=double(all(method{im}.p2x2(:)>=-aug_opts.tol));
    method{im}.ok_norm=double(abs(gc.norm-1)<aug_opts.tol);
    if (aug_opts.ifstd==1)
        mcorrs=method{im}.corrs;
    else
        mcorrs=btc_vec2corrs(btc_corrs2vec(method{im}.corrs,dict),dict_std); % in case nonstandard order-names used
    end
    switch method{im}.name
        case 'Pickard'
            method{im}.entropy_mrf=method{im}.entropy_mrf22;
        case 'DiagMRF'
            %
            % set up the correlations on the diagonally-interleaved lattices
            %
            b2=mcorrs.beta(3); % see btc_dimrf: horizontal beta corresponds to upper-right to lower-left
            b1=mcorrs.beta(4); % see btc_dimrf: vertical beta corresponds to upper-left to lower-right
            diagcorrs=[];
            diagcorrs.gamma=0; %gamma needs to be 0 so Pickard rules hold on the sublattices
            diagcorrs.theta=[0 0 0 0];
            diagcorrs.beta=[b1 b2 b1*b2 b1*b2]; %for Pickard rules to hold on sublattices
            diagcorrs.alpha=btc_meabb(b1,b2);
            method{im}.diagcorrs=diagcorrs;
            method{im}.diagp2x2=getp2x2_corrs(diagcorrs);
            method{im}.diagvec=btc_corrs2vec(diagcorrs);
            method{im}.diagcorrs=btc_vec2corrs(method{im}.diagvec,dict); %in case nonstandard order-names used
            % this is the entropy per unit area, since it is two independent lattices, each at half the density
            method{im}.entropy_mrf=getfield(getcorrs_p2x2(method{im}.diagp2x2,1),'entropy');
        case 'NoPickTT'
            method{im}.entropy_mrf=method{im}.entropy_mrf22; % an estimate is method{im}.entropy_mrf22 but it could be increased by Metropolis
        case 'NoPickBT'
            method{im}.entropy_mrf=method{im}.entropy_mrf22; % an estimate is method{im}.entropy_mrf22 but it could be changed by t-propagation and Metropolis
            %the 2x2 probs have to be OK, but also, the tee-probs, since
            %this propagates in a tee
            method{im}.ok_probs_2x2=method{im}.ok_probs;
            method{im}.ok_norm_2x2=method{im}.ok_norm;
            b=mcorrs.beta(posit(2));
            t=mcorrs.theta(posit(1));
            teecorrs=[];
            teecorrs.gamma=0;
            teecorrs.beta=[0 0 b 0];
            teecorrs.theta=[t*b^2 0 0 t];
            teecorrs.alpha=0;
            method{im}.teecorrs=teecorrs;
            method{im}.teeprobs=getp2x2_corrs(teecorrs);
            method{im}.ok_probs=double(all(method{im}.teeprobs(:)>=-aug_opts.tol));
            method{im}.ok_norm=double(abs(sum(method{im}.teeprobs(:))-1)<aug_opts.tol);
    end
    %
    %
    % if normalization is within tolerance, try to make it perfect
    if (method{im}.ok_probs==1) & (method{im}.ok_norm==1)
        p2x2=method{im}.p2x2;
        p2x2=max(p2x2,0);
        method{im}.p2x2=p2x2/sum(p2x2(:));
    end
    %
    if strcmp(method{im}.name,'Pickard');
        method{im}.ok_Pickard=0;
        if strcmp(method{im}.variant_lab,'NWSE');
            method{im}.ok_Pickard=(max(max(abs(gc.cig_conds(:,[1 4]))))<aug_opts.tol);
        end
        if strcmp(method{im}.variant_lab,'NESW');
            method{im}.ok_Pickard=(max(max(abs(gc.cig_conds(:,[2 3]))))<aug_opts.tol);
        end
        if (method{im}.ok_Pickard==0)
            disp(sprintf('Pickard condition violated for exptname %s and variant %s',exptname,method{im}.variant_lab));
            disp(spec);
        end
    end
end
%
augcoords.method=method;
augcoords.dict=dict;
augcoords.aug_opts=aug_opts;
augcoords.exptname=exptname;
augcoords.coorkinds=coorkinds;
return

function method=btc_augcoords_fillmeth(corrs,vnums,mname)
%fills in method values from correlations and Pickard branches
if (nargin<=2)
    mname='Pickard';
end
method=[];
for im=1:length(vnums)
    method{im}.name=mname;
    method{im}.corrs=corrs;
    method{im}.variant_num=vnums(im);
    switch mname
        case 'Pickard'
            if (vnums(im)==1)
                method{im}.variant_lab='NWSE';
            end
            if (vnums(im)==2)
                method{im}.variant_lab='NESW';
            end
        case 'DiagMRF'
            method{im}.variant_lab='std';
        case {'NoPickTT','NoPickBT'}
            method{im}.variant_lab=sprintf('rot%1.0f',method{im}.variant_num);
    end
end
return

function [corrs,vnums]=btc_augcoords_1var(coorkinds,spec,exptname,dict);
%get a method for a singleton variable
name1=dict.ordernames{1};
name2=dict.ordernames{2};
name2hv=cat(2,dict.ordernames{2},'_hv');
name2di=cat(2,dict.ordernames{2},'_diag');
name3=dict.ordernames{3};
name4=dict.ordernames{4};
%
n1=length(getfield(coorkinds,name1));
n2hv=length(getfield(coorkinds,name2hv));
n2di=length(getfield(coorkinds,name2di));
n3=length(getfield(coorkinds,name3));
n4=length(getfield(coorkinds,name4));
%
corrs=btc_vec2corrs(zeros(1,length(dict.codel)),dict);
%
%
for q=1:length(exptname)
    coordnum(q)=find(dict.codel==exptname(q));
    val(q)=getfield(spec,exptname(q));
    posit(q)=dict.posit(coordnum(q));
end
%
if (n1==1) %gamma->gamma, gamma^2, gamma^3, gamma^4
    corrs=setfield(corrs,name1,val(1));
    corrs=setfield(corrs,name2,repmat(val(1)^2,1,4));
    corrs=setfield(corrs,name3,repmat(val(1)^3,1,4));
    corrs=setfield(corrs,name4,val(1)^4);
    vnums=[1 2];
end
if (n2hv==1) %alpha gets two copies of a beta_hv
    corrs=setfield(corrs,name2,{posit(1)},val(1));
    corrs=setfield(corrs,name4,val(1)^2);
    vnums=[1 2];
end
if (n2di==1) %alpha is 0 for specified beta_di
    corrs=setfield(corrs,name2,{posit(1)},val(1));
    vnums=find(dict.inpickard(:,coordnum(1))==1);
end
if (n3==1)
    corrs=setfield(corrs,name3,{posit(1)},val(1));
    vnums=find(dict.inpickard(:,coordnum(1))==1);
end
if (n4==1)
    corrs=setfield(corrs,name4,val(1));
    vnums=[1 2];
end
return

function [corrs,vnums,mname]=btc_augcoords_2var(coorkinds,spec,exptname,dict);
%get a method for a pair of variables, neither of which are alpha
name1=dict.ordernames{1};
name2=dict.ordernames{2};
name2hv=cat(2,dict.ordernames{2},'_hv');
name2di=cat(2,dict.ordernames{2},'_diag');
name3=dict.ordernames{3};
name4=dict.ordernames{4};
%
n1=length(getfield(coorkinds,name1));
n2hv=length(getfield(coorkinds,name2hv));
n2di=length(getfield(coorkinds,name2di));
n3=length(getfield(coorkinds,name3));
n4=length(getfield(coorkinds,name4));
%
corrs=btc_vec2corrs(zeros(1,length(dict.codel)),dict);
vnums=[];
%
for q=1:length(exptname)
    coordnum(q)=find(dict.codel==exptname(q));
    val(q)=getfield(spec,exptname(q));
    posit(q)=dict.posit(coordnum(q));
end
if (n1==1) & (n2hv==1) %val(1) is beta_hv, val(2) is gamma -- parameters are for independent 1-d processes
    corrs=setfield(corrs,name1,val(2)); %gamma
    corrs=setfield(corrs,name2,repmat(val(2)^2,1,4)); %non-specified betas=gamma^2
    corrs=setfield(corrs,name2,{posit(1)},val(1)); %the specified beta
    corrs=setfield(corrs,name3,repmat(val(1)*val(2),1,4)); %theta=beta*gamma
    corrs=setfield(corrs,name4,val(1)^2); %alpha=beta^2
    vnums=[1 2];
    mname='Pickard';
end
if (n1==1) & (n2di==1) %val(1) is beta_di, val(2) is gamma -- parameters are for independent 1-d processes
    corrs=setfield(corrs,name1,val(2)); %gamma
    corrs=setfield(corrs,name2,repmat(val(2)^2,1,4)); %non-specified betas=gamma^2
    corrs=setfield(corrs,name2,{posit(1)},val(1)); %the specified beta
    corrs=setfield(corrs,name3,repmat(val(2)^3,1,4)); %start with theta=gamma^3
    if (posit(1)==3)
        thetalist=[2 3];
    end
    if (posit(1)==4)
        thetalist=[1 4];
    end
    corrs=setfield(corrs,name3,{thetalist(1)},val(1)*val(2)); %linked thetas are beta*gamma
    corrs=setfield(corrs,name3,{thetalist(2)},val(1)*val(2)); %linked thetas are beta*gamma
    corrs=setfield(corrs,name4,val(1)*val(2)^2); %alpha=beta*gamma^2
    vnums=find(dict.inpickard(:,coordnum(1))==1);
    mname='Pickard';
end
if (n2hv==2) %val(1) is one beta_hv, val(2) is the other
    corrs=setfield(corrs,name2,{posit(1)},val(1)); %the specified beta_hv
    corrs=setfield(corrs,name2,{posit(2)},val(2)); %the other specified beta_hv
    corrs=setfield(corrs,name2,{3},val(1)*val(2));
    corrs=setfield(corrs,name2,{4},val(1)*val(2));
    corrs=setfield(corrs,name4,btc_meabb(val(1),val(2)));         %alpha=maxent for beta_h, beta_v, beta_diag=beta_h*beta_v
    vnums=[1 2];
    mname='Pickard';
end
if (n2hv==1) & (n2di==1) %val(1) is beta_hv, val(2) is beta_di
    corrs=setfield(corrs,name2,{posit(1)},val(1)); %the specified beta_hv
    corrs=setfield(corrs,name2,{posit(2)},val(2)); %the specified beta_di
    if (posit(2)==3)
        thetalist=[2 3];
    end
    if (posit(2)==4)
        thetalist=[1 4];
    end
    corrs=setfield(corrs,name4,btc_meabb(val(1),0,val(2),0));         %alpha=maxent for beta_h, 0, beta_diag
    vnums=find(dict.inpickard(:,coordnum(2))==1);
    mname='Pickard';
end
%
if (n2di==2) %val(1) is one beta_di, val(2) is the other; do this by Pickard on sublattices
    corrs=setfield(corrs,name2,{posit(1)},val(1)); %one specified beta_di
    corrs=setfield(corrs,name2,{posit(2)},val(2)); %the other specified beta_di
    corrs=setfield(corrs,name4,val(1)*val(2));         %independent sublattices
    vnums=1;
    mname='DiagMRF';
end
%
if (n1==1) & (n3==1) %val(1) is theta, val(2) is gamma
    corrs=setfield(corrs,name1,val(2)); %gamma
    corrs=setfield(corrs,name2,repmat(val(2)^2,1,4)); %betas=gamma^2
    corrs=setfield(corrs,name3,repmat(val(2)^3,1,4)); %non-specified thetas=gamma^3
    corrs=setfield(corrs,name3,{posit(1)},val(1)); %the specified theta
    corrs=setfield(corrs,name4,val(1)*val(2)); %alpha=gamma*theta
    vnums=find(dict.inpickard(:,coordnum(1))==1);
    mname='Pickard';
end
if (n2hv==1) & (n3==1) %val(1) is theta, val(2) beta_hv
    corrs=setfield(corrs,name2,{posit(2)},val(2)); %the specified beta
    corrs=setfield(corrs,name3,{posit(1)},val(1)); %the specified theta
    corrs=setfield(corrs,name4,btc_meabt(val(2),val(1))); %adjust alpha
    vnums=find(dict.inpickard(:,coordnum(1))==1);
    mname='Pickard';
end
if (n2di==1) & (n3==1) %val(1) is theta, val(2) beta_di
    corrs=setfield(corrs,name2,{posit(2)},val(2)); %the specified beta
    corrs=setfield(corrs,name3,{posit(1)},val(1)); %the specified theta
    if all((dict.inpickard(:,coordnum(1))==dict.inpickard(:,coordnum(2))))
    %case in which the checks correlated by the diagonal beta is included in the checks of the theta
        vnums=find(dict.inpickard(:,coordnum(1))==1);
        mname='Pickard';
    else
    %case in which the checks correlated by the diagonal beta is ARE NOT in the checks of the theta
    %do this by Metropolis, seeding with random other than diagonal beta,
    %enforcing gamma=0 and h,v betas=0 by swapping algorithm
        if (posit(1)==1) & (posit(2)==3) %v, d
            vnums=1;
        end
        if (posit(1)==3) & (posit(2)==4) %u, e
            vnums=2;
        end
        if (posit(1)==4) & (posit(2)==3) %t, d
            vnums=3;
        end
        if (posit(1)==2) & (posit(2)==4) %w, e
            vnums=4;
        end
        corrs=setfield(corrs,name3,{5-posit(1)},val(1)*val(2)); %the opposite theta is beta*theta
        mname='NoPickBT';
    end
end
if (n3==2) %val(1) and val(2) are thetas
    corrs=setfield(corrs,name3,{posit(1)},val(1)); %one specified theta
    corrs=setfield(corrs,name3,{posit(2)},val(2)); %the other specified theta
    if all((dict.inpickard(:,coordnum(1))==dict.inpickard(:,coordnum(2))))
    %case in which the thetas share a common corner (Pickard)
        vnums=find(dict.inpickard(:,coordnum(1))==1);
        mname='Pickard';
    else
    %case in which the thetas share a common edge
    %do this by Metropolis, enforcing betas and gammas=0 by swapping algorithm
        if (length(intersect(posit,[3 4])))==2 %t, u
            vnums=1;
        end
        if (length(intersect(posit,[2 4])))==2 %t, w
            vnums=2;
        end
        if (length(intersect(posit,[1 2])))==2 %v, w
            vnums=3;
        end
        if (length(intersect(posit,[1 3])))==2 %u, v
            vnums=4;
        end
        mname='NoPickTT';
    end
end
return

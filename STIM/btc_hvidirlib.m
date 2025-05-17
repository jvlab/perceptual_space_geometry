function lib=btc_hvidirlib(opts);
% lib=btc_hviirlib(opts) generates a list of eigenvectors in the horizontal-vertical invert
% subspace of the btc coords
%
% if opts is present: source of eigenvector data
%     (see btc_eivecs_def for other definitions aqnd defaults)
%     Note that opts.symtype forced to be 'hvi'
%   if file cannot be found or data are not available, generates a warning and uses builtins
%
% if opts is absent, uses builtins
%
% builtin: empirical eigenvectors from model 12 via btc_eivecs_stats
% in SYMMETRIC basis (here transformed back)
%
% builtin values can be recovered for each subject by, e.g.,
% btc_symdirlib(setfields([],{'subjid','raw_res'},{'mc',0}))
%
% also includes Minkowski directions from btc_morphological.doc and btc_euler_table.doc
%
%   See also:  BTC_SYMBASIS, BTC_EIVEC_STATS, BTC_HVIDIRLIB, BTC_EIVECS_READ, BTC_EIVECS_DEF.
%
lib=[];
hvi_eivecs=[];
hvi_eivals=[];
if (nargin>0)
    opts=btc_eivecs_def(opts);
    opts.symtype='hvi';
    [evecs,evals,errs,estruct,ou]=btc_eivecs_read(opts);
    if isempty(errs)
        fieldname=cat(2,estruct.subjid,zpad(opts.variant,2));
        if (opts.raw_res==1)
            fieldname=cat(2,fieldname,'_res');
        end
        hvi_eivecs.(fieldname)=evecs;
        hvi_eivals.(fieldname)=evals;
    else
        disp('Warning: problem accessing eigenvector database, using builtins instead.');
        disp(errs)
    end 
end
%
if (isempty(hvi_eivecs) | isempty(hvi_eivals))
    %use builtins
    hvi_eivecs.mc12=...
      [-0.9769    0.2138;...
       -0.2138   -0.9769];
    hvi_eivals.mc12=[3.8384; -0.2823];
    %
    hvi_eivecs.dt12=...
      [-0.9773    0.2117;...
       -0.2117   -0.9773];
    hvi_eivals.dt12=[3.9621;  0.0752];
    %
    hvi_eivecs.df12=...
      [-0.9816    0.1911;...
       -0.1911   -0.9816];
    hvi_eivals.df12=[4.3510;  0.2664];
    %
    hvi_eivecs.jd12=...
      [-0.9857    0.1688;...
       -0.1688   -0.9857];
hvi_eivals.jd2=[3.1055;  -0.0364];
end %file not available, use builtins
%
%transform the eigenvectors back into standard btc coordinates
names=fieldnames(hvi_eivecs);
for iname=1:length(names)
    eivecs=diag([1/sqrt(2) 1/2])*hvi_eivecs.(names{iname});
    for k=1:size(eivecs,1)
        n=length(lib);
        vals=eivecs(:,k);
        mv=vals(min(find(abs(vals)==max(abs(vals)))));
        if (mv<0)
            vals=-vals;
        end
        lib{n+1}.desc=sprintf('eivec hvi%1.0f (%s)',k,names{iname});
        lib{n+1}.vals=vals';
    end
end
%
% directions corresponding to the Minkowski invariants
%
n=length(lib);
lib{n+1}.desc='Minkowski6L-R';
lib{n+1}.vals=[0 -1]/8;
%
for k=1:length(lib)
    lib{k}.valstring=sprintf('d(-e)=%5.2f t(-u)v(-w)=%5.2f',lib{k}.vals);
end
return

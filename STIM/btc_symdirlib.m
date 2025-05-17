function lib=btc_symdirlib(opts);
% lib=btc_symdirlib(opts) generates a list of eigenvectors in the symmetric
% subspace of the btc coords
%
% if opts is present: source of eigenvector data
%     (see btc_eivecs_def for other definitions aqnd defaults)
%     Note that opts.symtype forced to be 'sym'
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
sym_eivecs=[];
sym_eivals=[];
if (nargin>0)
    opts=btc_eivecs_def(opts);
    opts.symtype='sym';
    [evecs,evals,errs,estruct,ou]=btc_eivecs_read(opts);
    if isempty(errs)
        fieldname=cat(2,estruct.subjid,zpad(opts.variant,2));
        if (opts.raw_res==1)
            fieldname=cat(2,fieldname,'_res');
        end
        sym_eivecs.(fieldname)=evecs;
        sym_eivals.(fieldname)=evals;
    else
        disp('Warning: problem accessing eigenvector database, using builtins used instead.');
        disp(errs)
    end 
end
if (isempty(sym_eivecs) | isempty(sym_eivals))
    %use builtins
    sym_eivecs.mc12=...
      [-0.9708   -0.1192    0.0992   -0.1583   -0.0915;...
        0.1221   -0.9355    0.1428    0.1822   -0.2373;...
       -0.0021    0.1214    0.8948    0.3105    0.2971;...
       -0.2017   -0.0591   -0.4092    0.7569    0.4641;...
        0.0433   -0.3039   -0.0415   -0.5219    0.7947];
    sym_eivals.mc12=[47.5045; 17.9333 ;12.2834; 3.1483; 1.4038];
    %
    sym_eivecs.dt12=...
      [-0.9814   -0.0156   -0.0507   -0.1579   -0.0959;...
        0.0058   -0.9647    0.0825    0.1529   -0.1980;...
       -0.0051   -0.0763   -0.9524    0.2452    0.1640;...
       -0.1920    0.0556    0.2862    0.8113    0.4690;...
       -0.0024   -0.2454    0.0398   -0.4831    0.8395];
    sym_eivals.dt12=[31.5176; 15.2701; 7.3608; 2.0814;  0.5708];
    %
    sym_eivecs.df12=...
      [-0.9778   -0.0458    0.2024    0.0173    0.0221;...
        0.0909   -0.9555    0.2229   -0.1324    0.1072;...
       -0.1865   -0.2183   -0.9109   -0.2338   -0.1822;...
        0.0258    0.1157    0.2767   -0.7873   -0.5381;...
        0.0101   -0.1544    0.0557    0.5547   -0.8157];
    sym_eivals.df12=[33.1123; 18.9530; 14.1614; 4.2306; 0.5748];
    %
    sym_eivecs.jd12=...
      [-0.9961    0.0285   -0.0358    0.0593   -0.0469;...
        0.0213    0.9761    0.1143   -0.1381   -0.1206;...
        0.0164    0.0977   -0.9565   -0.2294    0.1507;...
       -0.0832   -0.0979    0.2611   -0.8581    0.4231;...
       -0.0129    0.1649    0.0517    0.4341    0.8840];
    sym_eivals.jd12=[24.2985; 13.2771; 6.8577; 1.9919; 0.9403];
end %isempty
%
%transform the eigenvectors back into standard btc coordinates
names=fieldnames(sym_eivecs);

for iname=1:length(names)
    eivecs=diag([1 1/sqrt(2) 1/sqrt(2) 1/2 1])*sym_eivecs.(names{iname});
    for k=1:size(eivecs,1)
        n=length(lib);
        vals=eivecs(:,k);
        mv=vals(min(find(abs(vals)==max(abs(vals)))));
        if (mv<0)
            vals=-vals;
        end
        lib{n+1}.desc=sprintf('eivec sym%1.0f (%s)',k,names{iname});
        lib{n+1}.vals=vals';
    end
end
%
% directions corresponding to the Minkowski invariants
%
n=length(lib);
lib{n+1}.desc='Minkowski4';
lib{n+1}.vals=[-4 -2 +1 +1 +1]/16;
n=length(lib);
lib{n+1}.desc='Minkowski8';
lib{n+1}.vals=[-4 +2 -1 +1 -1]/16;
n=length(lib);
lib{n+1}.desc='Minkowski4+8';
lib{n+1}.vals=[-4  0  0 +1  0]/8;
n=length(lib);
lib{n+1}.desc='Minkowski4-8';
lib{n+1}.vals=[ 0 -4 +2  0 +2]/16;
%
for k=1:length(lib)
    lib{k}.valstring=sprintf('g=%5.2f bc=%5.2f de=%5.2f tuvw=%5.2f a=%5.2f',lib{k}.vals);
end
return

% btc_test
% tests the modules of the btc series for texture generation
%
%    See also: BTC_DEFINE, BTC_CORRS2VEC, BTC_VEC2CORRS, BTC_VEC2LETCODE, BTC_LETCODE2VEC,
%      BTC_EXPTNAME, BTC_COORKINDS, BTC_AUGCOORDS.
%
if ~exist('opts') 
    opts=[];
end
opts=filldefault(opts,'ifshow',1);
opts=filldefault(opts,'ifshow_fig',0);
opts=filldefault(opts,'ordernames',{'gamma','beta','theta','alpha'});
%
aug_opts=[];
%
disp(opts)
disp(' exit any time and re-define opts if desired.');
%
dict=btc_define(opts);
%
for k=1:length(dict.checks)
    disp(sprintf('variable %2.0f->%1s (%s)',k,dict.codel(k),dict.name{k}));
end
coordnums=getinp('which coordinate(s) to manipulate, typically one or two, 0 for none','d',[0 length(dict.checks)]);
if (coordnums==0) coordnums=[]; end
coordnums=unique(coordnums);
spec_init=btc_vec2letcode(repmat(NaN,1,length(dict.checks)),dict);
for ic=1:length(coordnums)
    k=coordnums(ic);
    spec_init=setfield(spec_init,dict.codel(k),getinp(sprintf('value for %1s',dict.codel(k)),'f',[-1 1],0));
end
results.spec_vec=btc_letcode2vec(spec_init,dict);
% remove the NaN's
results.spec_stripped=btc_letcode_strip(spec_init,dict);
% create the standard experiment string
results.spec_exptname=btc_exptname(char(fieldnames(results.spec_stripped)),dict);
% determine the kinds of coordinates
results.spec_coorkinds=btc_coorkinds(char(fieldnames(results.spec_stripped)),dict);
% attempt to augment the coordinates
results.augcoords=btc_augcoords(spec_init,dict,aug_opts);
%
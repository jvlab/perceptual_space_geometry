%btc_qform_customize_test: framework for customizing a quadartic form model
%
%  See btc_qform_customize_notes.docx
%
%  See also:  BTC_QFORM_CUSTOMIZE_RUN, BTC_QFORM_CUSTOMIZE, BTC_SOID_FIND.
%
dict=btc_define;
nbtc=length(dict.codel);
%thresholds for debugging -- to be imported from subject-specific threshold data
thr.a=0.6;
thr.t=0.9;
thr.d=0.4;
thr.b=0.25;
thr.g=0.15;
%
qfm_model_file='btc_allraysfixedb_avg_100surrs_madj.mat';
r_orig=getfield(load(qfm_model_file),'r');
imodel=12;
qform_orig=r_orig{imodel}.results.qfit;
%
axes_to_fix=fieldnames(thr);
nfix=length(axes_to_fix);
%
disp(sprintf(' customizing model %2.0f, label %s',imodel,r_orig{imodel}.setup.label));
%
dsq=1;
for k=1:nfix
    ax=axes_to_fix{k};
    spec=struct;
    spec.(ax)=1;
    [specx,avec,res_fzero,ou]=btc_soid_find(spec,dsq,qform_orig,dict);
    disp(sprintf(' axis %s thr desired %8.5f un-customized thr %8.5f',ax,thr.(ax),specx.(ax)));
end
[qform_new,cust_factors]=btc_qform_customize(qform_orig,thr,dict);
%
cust_factors
%verify that qform_new yields desired thresholds
for k=1:nfix
    ax=axes_to_fix{k};
    spec=struct;
    spec.(ax)=1;
    [specx,avec,res_fzero,ou]=btc_soid_find(spec,dsq,qform_new,dict);
    disp(sprintf(' axis %s thr desired %8.5f    customized thr %8.5f',ax,thr.(ax),specx.(ax)));
end 

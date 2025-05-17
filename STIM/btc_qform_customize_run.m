%btc_qform_customize_run: customize a file of quadratic form models
% for individual subjects' thresholds, including surrogates
% 
% Threshold data from btc_thresholds_*subjs_*.xls
%
%  See btc_qform_customize_notes.docx
%
% 11Feb25: add subjects SN and NF from btc_thresholds_29subjs_11Feb25.xls
%
%  See also:  BTC_QFORM_CUSTOMIZE_TEST, BTC_QFORM_CUSTOMIZE.
%
dict=btc_define;
nbtc=length(dict.codel);
%
fields_keep_r={'setup','ou_dict','edirs'};
fields_keep_results={'nplanes','planes','axes','qsetup','soid_fit','which_axes'};
%
%thresholds imported from subject-specific threshold data in
%btc_thresholds_27subjs_31Jul23.xls (thresholds along negative and postive rays)
%
if ~exist('thrs')
    thr_pairs='gbdta';
    raw_thrs.MC =[0.15	0.14	0.24	0.27	0.36	0.35	0.64	0.64	0.46	0.57];
    raw_thrs.DT =[0.19	0.17	0.26	0.26	0.41	0.42	0.81	0.94	0.60	0.72];
    raw_thrs.DF =[0.17	0.17	0.26	0.25	0.36	0.36	0.72	0.88	0.58	0.64];
    raw_thrs.JD	=[0.20	0.20	0.29	0.30	0.42	0.45	0.90	0.98	0.65	0.81];
    raw_thrs.KP	=[0.21	0.19	0.30	0.34	0.46	0.45	0.75	0.87	0.64	0.92];
    raw_thrs.SR	=[0.16	0.17	0.26	0.24	0.48	0.39	0.62	0.75	0.55	0.92];
    raw_thrs.RS	=[0.17	0.15	0.34	0.31	0.47	0.40	0.81	1.05	0.68	0.89];
    raw_thrs.SP	=[0.15	0.14	0.24	0.22	0.36	0.33	0.58	0.65	0.50	0.47];
    raw_thrs.WC	=[0.20	0.19	0.33	0.34	0.45	0.53	0.95	0.98	0.53	0.77];
    raw_thrs.ZA	=[0.20	0.20	0.28	0.39	0.39	0.53	1.22	1.19	0.59	0.91];
    raw_thrs.JWB=[0.18	0.15	0.25	0.25	0.39	0.40	0.82	0.95	0.53	0.79];
    raw_thrs.IL =[0.25	0.27	0.30	0.32	0.44	0.50	1.40	1.64	0.61	0.99];
    raw_thrs.EFV=[0.22	0.18	0.22	0.24	0.36	0.42	1.00	1.00	0.48	0.64];
    raw_thrs.YCL=[0.23	0.22	0.27	0.28	0.40	0.48	0.97	1.09	0.56	0.81];
    raw_thrs.PJ =[0.19	0.20	0.39	0.37	0.50	0.51	0.97	1.00	0.64	1.05];
    raw_thrs.AJT=[0.32	0.34	0.48	0.48	0.64	0.61	0.89	1.12	0.80	1.73];
    raw_thrs.BL	=[0.20	0.18	0.35	0.31	0.50	0.52	0.92	1.07	0.58	0.87];
    raw_thrs.ZK	=[0.24	0.23	0.42	0.41	0.50	0.53	0.94	1.05	0.71	1.10];
    raw_thrs.SN	=[0.24	0.22	0.35	0.43	0.51	0.48	1.43	1.34	0.74	0.94];
    raw_thrs.NF	=[0.19	0.19	0.26	0.27	0.38	0.40	0.79	0.83	0.53	0.65];
    subjs=sort(fieldnames(raw_thrs));
    thrs=struct;
    for isubj=1:length(subjs)
        for k=1:length(thr_pairs)
            thrs.(lower(subjs{isubj})).(thr_pairs(k))=mean(raw_thrs.(subjs{isubj})(2*(k-1)+[1 2]));
        end
    end
end
qfm_template='btc_allraysfixedb_XX_100surrs_madj.mat';
qfm_customized_template='btc_allraysfixedb_XX-YY_100surrs_madj.mat';
if_ok=0;
while (if_ok==0)
    qfm_orig_name=getinp('original quadratic form source (typically ''avg'')','s',[],'avg');
    %
    qfm_model_file=strrep(qfm_template,'XX',qfm_orig_name);
    srcs=fieldnames(thrs);
    for k=1:length(srcs)
        disp(sprintf('%2.0f->subject ID %s',k,srcs{k}));
    end
    subj_ptr=getinp('choice','d',[1 length(srcs)]);
    qfm_customized_name=srcs{subj_ptr};
    thr=thrs.(qfm_customized_name);
    disp('thresholds to fit')
    disp(thr);
    qfm_customized_file=strrep(strrep(qfm_customized_template,'XX',qfm_orig_name),'YY',qfm_customized_name);
    %
    disp(sprintf('         model file to be read: %s',qfm_model_file));
    disp(sprintf(' customized file to be created: %s',qfm_customized_file));
    if_ok=getinp('1 if ok','d',[0 1])
end
r_orig=getfield(load(qfm_model_file),'r');
disp(r_orig);
r=cell(1,length(r_orig));
%
for im=1:length(r_orig)
    r{im}=struct;
    for k=1:length(fields_keep_r)
        r{im}.(fields_keep_r{k})=r_orig{im}.(fields_keep_r{k});
    end
    r{im}.results=struct;
    for k=1:length(fields_keep_results)
        r{im}.results.(fields_keep_results{k})=r_orig{im}.results.(fields_keep_results{k});
    end
    nsurrs=size(r_orig{im}.results.surrogates.qfit,3);
    %
    r{im}.results.qfit_orig=r_orig{im}.results.qfit;
    r{im}.results.surrogates_orig=r_orig{im}.results.surrogates;
    r{im}.results.surrogates.qfit=zeros(nbtc,nbtc,nsurrs);
    %
    disp(sprintf(' customizing model %2.0f and %4.0f surrogates: label %s',im,nsurrs,r_orig{im}.setup.label));
    qform_orig=r_orig{im}.results.qfit;
    [qform_new,cust_factors]=btc_qform_customize(qform_orig,thr,dict);
    r{im}.results.qfit=qform_new;
    r{im}.results.cust_factors=cust_factors;
    disp(' customization factors')
    disp(cust_factors);
    for isurr=1:nsurrs
        r{im}.results.surrogates.qfit(:,:,isurr)=btc_qform_customize(r_orig{im}.results.surrogates.qfit(:,:,isurr),thr,dict);
    end
end
%save with new file name
if_write=getinp(sprintf('1 if ok to write files %s',qfm_customized_file),'d',[0 1]);
if if_write
    save(qfm_customized_file,'r');
end

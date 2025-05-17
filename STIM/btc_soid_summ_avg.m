% btc_soid_summ_avg:  summarize and average results of btc_soid_demo
% and save results in an "average" subject file, usable for
% btc_pred_emink.
%
% Cross-subject pooling is done by simple averaging.
% Each subject's datafile is assumed to have the same number of surrogates, and same setup.
%
%   See also:  BTC_SOID_DEMO, BTC_EIVECS_STATS, COMPSTRUCT.
%
if ~exist('dict') dict=btc_define([]); end
nbtc=length(dict.codel); %10
%
fields_match={'setup','ou_fit','ou_dict','edirs'}; %fields that should match across subjects
fields_scalar={'rmse','mean_data_minus_fit_on_axis','mean_data_minus_fit_off_axis','rmse_med_surrogates'};
fields_results={'results'};
fields_omit={'edirs_data_and_fits'};
%
rfields_match={'nplanes','planes','axes','qsetup','soid_fit','which_axes','regressors'};
rfields_scalar={'condn_raw','condn_adjcol','condn_adjrow'};
rfields_special={'b','qfit','surrogates'}; %the fields that require special treatment
rfields_omit={'errs','x','augvecs_all','augvecs_sel'};
%
etag='_each';
%
if ~exist('results_fn_template') results_fn_template='btc_allraysfixedb_xx_100surrs'; end
disp('append ''_madj'' to end to access the files with diagonals adjusted.');
results_fn_template=getinp('name of results file name template from btc_soid_test','s',[],results_fn_template);
%
if_verbose=getinp('1 for verbose behavior','d',[0 1],0);
subjids={'mc','dt','df','jd'};
for k=1:length(subjids)
    disp(sprintf('%1.0f->%2s',k,subjids{k}));
end
subjid_list=getinp('subject id list','d',[1 length(subjids)],[1:length(subjids)]);
%
r_in=cell(0);
r=cell(0);
for subjid_ptr=1:[length(subjid_list)] %top value is cross-subject average
    subjid=subjids{subjid_list(subjid_ptr)};
    results_fn=strrep(results_fn_template,'xx',subjid);
    r_in{subjid_ptr}=getfield(load(results_fn),'r');
end
nvariants=length(r_in{1});
r=[];
if_mismatch=0;
for subjid_ptr=1:length(subjid_list)
    subjid=subjids{subjid_list(subjid_ptr)};
    for ivariant=1:nvariants
        r{ivariant}.subjid{subjid_ptr}=subjid;
        disp(sprintf('processing subject %1.0f (%s), variant %2.0f',subjid_ptr,subjid,ivariant));
        fnames=fieldnames(r_in{subjid_ptr}{ivariant});
        for ifn=1:length(fnames)
            fname=fnames{ifn};
            if ~isempty(strmatch(fname,fields_match))
                %disp(sprintf(' field that should match across subjects: %s',fname));
                if (subjid_ptr==1)
                    r{ivariant}.(fname)=r_in{subjid_ptr}{ivariant}.(fname);
                    if (if_verbose)
                        disp(sprintf(' value from %s for %s installed (others should match).',subjids{subjid_list(1)},subjid,fname));
                    end
                else
                    ifdif=compstruct(subjids{subjid_list(1)},r{ivariant}.(fname),subjid,r_in{subjid_ptr}{ivariant}.(fname));
                    if isempty(ifdif)
                        if (if_verbose)
                            disp(sprintf(' value from %s and %s for %s match.',subjids{subjid_list(1)},subjid,fname));
                        end
                    else
                        disp(sprintf(' value from %s and %s for %s DO NOT match.',subjids{subjid_list(1)},subjid,fname));
                        if_mismatch=if_mismatch+1;
                    end
                end
            elseif ~isempty(strmatch(fname,fields_scalar))
                %disp(sprintf(' field that should be treated as a scalar: %s',fname));
                r{ivariant}.(cat(2,fname,etag))(1,subjid_ptr)=r_in{subjid_ptr}{ivariant}.(fname);
                if (if_verbose)
                    disp(sprintf(' value from %s for %s installed.',subjid,fname));
                end
            elseif ~isempty(strmatch(fname,fields_results))
                %disp(sprintf(' results field: %s',fname));
                %r{ivariant}.results=r_in{subjid_ptr}{ivariant}.results;
                rfnames=fieldnames(r_in{subjid_ptr}{ivariant}.results);
                for irfn=1:length(rfnames)
                    rfname=rfnames{irfn};
                    if ~isempty(strmatch(rfname,rfields_match))
                        if (subjid_ptr==1)
                            r{ivariant}.results.(rfname)=r_in{subjid_ptr}{ivariant}.results.(rfname);
                            if (if_verbose)
                                disp(sprintf(' results value from %s for %s installed (others should match).',subjids{subjid_list(1)},subjid,rfname));
                            end
                        else
                            ifdif=compstruct(subjids{subjid_list(1)},r{ivariant}.results.(rfname),subjid,r_in{subjid_ptr}{ivariant}.results.(rfname));
                            if isempty(ifdif)
                                if (if_verbose)
                                    disp(sprintf(' results value from %s and %s for %s match.',subjids{subjid_list(1)},subjid,rfname));
                                end
                            else
                                disp(sprintf(' results value from %s and %s for %s DO NOT match.',subjids{subjid_list(1)},subjid,rfname));
                                if_mismatch=if_mismatch+1;
                            end
                        end
                    elseif ~isempty(strmatch(rfname,rfields_scalar))
                    %disp(sprintf(' field that should be treated as a scalar: %s',fname));
                        r{ivariant}.results.(cat(2,rfname,etag))(1,subjid_ptr)=r_in{subjid_ptr}{ivariant}.results.(rfname);
                        if (if_verbose)
                            disp(sprintf(' results value from %s for %s installed.',subjid,fname));
                        end
                    elseif ~isempty(strmatch(rfname,rfields_omit))
                        if (if_verbose)
                            disp(sprintf(' results field to be omitted from summary: %s',rfname));
                        end
                    elseif ~isempty(strmatch(rfname,rfields_special))
                        switch rfname
                            case 'b'
                                r{ivariant}.results.(cat(2,rfname,etag))(:,subjid_ptr)=r_in{subjid_ptr}{ivariant}.results.(rfname);
                            case 'qfit'
                                r{ivariant}.results.(cat(2,rfname,etag))(:,:,subjid_ptr)=r_in{subjid_ptr}{ivariant}.results.(rfname);
                            case 'surrogates'
                                r{ivariant}.results.(rfname).(cat(2,'b',etag))(:,:,subjid_ptr)=r_in{subjid_ptr}{ivariant}.results.(rfname).b;
                                r{ivariant}.results.(rfname).(cat(2,'qfit',etag))(:,:,:,subjid_ptr)=r_in{subjid_ptr}{ivariant}.results.(rfname).qfit;
                        end
                        if (if_verbose==1)
                            disp(sprintf(' results field added to summary for averaging: %s',rfname));
                        end
                    else
                        disp(sprintf(' unknown results field: %s',rfname));
                    end
                end %irfn
            elseif ~isempty(strmatch(fname,fields_omit))
                if (if_verbose)
                    disp(sprintf(' field to be omitted from summary: %s',fname));
                end
            else
                disp(sprintf(' unknown field: %s',fname));
            end
        end %ifn
    end %ivariant
end %subjid_ptr
disp(sprintf('number of mismatches: %4.0f',if_mismatch));
disp('merging results')
for ivariant=1:nvariants
    for ifn=1:length(fields_scalar)
        fname=fields_scalar{ifn};
        r{ivariant}.(fname)=mean(r{ivariant}.(cat(2,fname,etag)));
    end
    for irfn=1:length(rfields_scalar)
        rfname=rfields_scalar{irfn};
        r{ivariant}.results.(rfname)=mean(r{ivariant}.results.(cat(2,rfname,etag)));
    end
    r{ivariant}.results.b=mean(r{ivariant}.results.b_each,2);
    r{ivariant}.results.qfit=mean(r{ivariant}.results.qfit_each,3);
    r{ivariant}.results.surrogates.b=mean(r{ivariant}.results.surrogates.b_each,3);
    r{ivariant}.results.surrogates.qfit=mean(r{ivariant}.results.surrogates.qfit_each,4);
end %ivariant
%
if getinp('1 to write the merged file','d',[0 1])
    results_fn_out=strrep(results_fn_template,'xx','avg');
    results_fn_out=getinp('merged file name','s',[],results_fn_out);
    save(results_fn_out,'r');
end

%psg_like_maketable: create a table from the automated outputs of psg_[umi_trip|tent]_like_demo.
%
% creates table_like, a table of all log likelihoods, surrogate values, a
% priori values for sym/trans, umi, and addtree analyses, and the Dirichlet
% parameters and their likleihoods
%
% some categorical variables are converted to tokens, as defined in the
% 'tokens' structure and kept as UserData in table_like.
%
% llr quantities for umi are corrected, i.e., have log(h) subtracted
%
%   See also:  PSG_UMI_TRIPLIKE_DEMO, PSG_TENTLIKE_DEMO, PST_UMI_TRIP_LIKE_RUN.
%
if ~exist('table_like') 
    disp('Creating table_like from scratch.');
    table_like=table();
else
    disp('Appending to existing table_like');
end
if_done=0;
n_processed=0;
ds_types_allowed={'su','adt'};
paradigm_types_allowed={'animals','btc'};
tokens=struct;
tokens.ipchoice={'fixed h','fitted h'};
tokens.llr_type={'sym','umi','adt'};
tokens.thr_type={'min','max','avg'};
while if_done==0
    fn=getinp('likelihood analysis file to add to table, e.g., [animals|btc]_[umi_trip|tent]like_db.mat','s',[]);
    ds_list=0;
    if exist(fn,'file')
        f=load(fn);
        if isfield(f,'db')
            ds_names=fieldnames(f.db);
            for ids=1:length(ds_names)
                disp(sprintf('ds %3.0f->%s',ids,ds_names{ids}'));
            end
            ds_list=getinp('datasets to use (0 to re-enter a file or end)','d',[0 length(ds_names)],[1:length(ds_names)]);
            ds_list=sort(ds_list);
        else
            disp(sprintf('%s does not contain a variable db',fn));
        end
    else
        disp(sprintf('%s not found',fn));
    end
    paradigm_type=fn([1:min(find(fn=='_'))-1]);
    if isempty(strmatch(paradigm_type,paradigm_types_allowed,'exact'))
        disp('unrecognized file type');
        ds_list=0;
    end
    if max(ds_list)>0
        for ids_ptr=1:length(ds_list)
            ids=ds_list(ids_ptr);
            ds_name=ds_names{ids};
            %get paradigm and subject id
            underscores=find(ds_name=='_');
            switch paradigm_type
                case 'animals'
                   paradigm_name=ds_name(underscores(1)+1:underscores(2)-1);
                   subj_id=ds_name(1:underscores(1)-1);
                case 'btc'
                    paradigm_name=ds_name(1:underscores(1)-1);
                    subj_id=ds_name(underscores(2)+1:underscores(3)-1);
            end
            table_strings=array2table({paradigm_type,paradigm_name,subj_id});
            table_strings.Properties.VariableNames={'paradigm_type','paradigm_name','subj_id'};
            %
            s=f.db.(ds_name).s;
            r=f.db.(ds_name).r;
            h_fixlist=r.h_fixlist;
            dirichlet=r.dirichlet;
            r_fields=fieldnames(rmfield(rmfield(r,'h_fixlist'),'dirichlet'));
            ds_type=r_fields{1}; %adt from psg_tentlike_demo, or su from psg_umi_triplike_demo.
            if isempty(strmatch(ds_type,ds_types_allowed,'exact'))
                disp(sprintf('type %s unrecognized and not processed',ds_type));
            else
                n_processed=n_processed+1;
                switch ds_type
                    case 'su'
                        llr_types_thisdb={'sym','umi'};
                    case 'adt'
                        llr_types_thisdb={'adt'};
                end
                for illr=1:length(llr_types_thisdb)
                    llr_type=llr_types_thisdb{illr};
                    illr_uid=strmatch(llr_type,tokens.llr_type,'exact');
                    for ipchoice=1:length(tokens.ipchoice)
                        %
                        if isfield(s{ipchoice},llr_types_thisdb{illr})
                            s_data=s{ipchoice}.(llr_types_thisdb{illr});
                            for ithr=1:length(tokens.thr_type)
                                table_common=array2table([illr_uid ipchoice,ithr,s_data.params.a,s_data.params.h,s_data.ah_llr,s_data.apriori_vals]);
                                table_common.Properties.VariableNames={'llr_type','ipchoice','thr_type','a','h','ah_llr','apriori_llr'};
                                s_thr=s_data.thr_type{ithr};
                                nfracs=length(s_thr.thr_ptr_use);
                                thr_ptr_use=s_thr.thr_ptr_use;
                                table_data=array2table([s_thr.frac_keep_list(1:nfracs),s_thr.tally_table(thr_ptr_use,:),s_thr.means_per_set_adj(thr_ptr_use,:),s_thr.eb_stds(thr_ptr_use,:)]);
                                table_data.Properties.VariableNames={'frac_keep','thr_val','ntriads','ntrials','llr_data','llr_flip_all','llr_flip_any','llr_data_sd','llr_flip_all_sd','llr_flip_any_sd'};
                                table_row=[table_strings,table_common];
                                table_like=[table_like;[repmat(table_row,nfracs,1) table_data]];
                            end
                        end %does ipchoice have this llr type?
                    end %ipchoice
                end %illr
                disp(sprintf('processed set %3.0f in table: %30s (set %3.0f in %30s), type %s',n_processed,ds_names{ids},ids,fn,ds_type));               
            end
        end
    end
    disp(sprintf('table_like now has %4.0f rows',size(table_like,1)));
    if_done=getinp('1 if done','d',[0 1]);
end
table_like.Properties.UserData.tokens=tokens;
table_like.Properties.UserData.notes='llr quantitites for umi have log(h) subtracted';
disp('remember to save table_like')
    
    
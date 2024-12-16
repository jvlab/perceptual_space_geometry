%psg_raystats_dbplot: plot ray statstics across subjects, paradigms, etc.
%
% Reads one or more tables created by psg_raystats_summ, analyzes, and plots
% After running, can also save raystats_*.mat t_meta_all t_all to save the concatenated data tables
% 
% 12Dec24: start plots by dimension
%
%  See also: PSG_RAYSTATS_SUMM, PSG_LLFITS_SUMM, TABLECOL2CHAR, PSG_RAYSTATS_DBPLOT_TICKS
%   PSG_RAYSTATS_DBPLOT_STYLE.
%
if ~exist('ui_filter') ui_filter='raystats_*.mat'; end
if ~exist('ui_filter_gen') ui_filter_gen='[raystats|llfits|*]_*_ddmmmyy.mat'; end
if_replace_avg=getinp('1 to replace qform[-avg]-XX by XX for qform models','d',[0 1],1);
%
criterion_names={'subj_model_ID','expt_grp','expt_uid'}; %ways to group or plot
ncrits=length(criterion_names);
%
plot_types={...
    'plot a geometric parameter (rows) for a specific model dimension (columns)',...
    'plot a geometric parameter (columns) for a specific model dimension (rows)',...
    'plot a geometric parameter (rows) as function of dimension',...
    'plot a geometric parameter (columns) as function of dimension'};
%
nplot_types=length(plot_types);
if ~exist('plot_major_space') plot_major_space=1; end
if ~exist('plot_minor_space') plot_minor_space=0.3; end
if ~exist('plot_label_height') plot_label_height=0.9; end
if ~exist('plot_ebhw') plot_ebhw=0.3*plot_minor_space; end
%
values_short.brightness='bright';
values_short.constrained_grouping='paired';
values_short.unconstrained_grouping='group';
values_short.similarity='standard';
values_short.dis_similarity='dissim';
values_short.working_memory='workmem';
values_short.threshold='thresh';
%
%read one or more data table and keep originals
%optionally apply relabelling of model id ('avg-bl' -> 'bl','avg-zk'->'zk')
%
%table definitions for t_meta_all
meta_variable_names={'psy_model','subj_model_ID','expt_grp','expt_name','expt_param','expt_uid','coords_source','coords_file','choices_source','choices_file','sess_range'};
%
% fields that do not change within a dataset
% psy_model: {'psy'|'mdl'}
% subj_model_ID: 'mc', ..., for a subject; 'qfm' for a generic quadratic model, qfm_mc for a customized model
% expt_grp: {'threshold,'similarity','dis_similarity','working_memory','constrained_grouping','unconstrained_grouping','brightness'};
% expt_name: {'','','dis','wm','gp','gm',br'}; 
% expt_param: NaN for anything but wm, 1000 for wm1000
% expt_uid: 'bc6pt','bc55qpt','bcpm3pt','bgca3pt',etc.
% coords_source: if psychophysical datas: coordinate file name with path, if qform model, a label that has file name and params
% coords_file: coordinate file name
% choices_source: choice file name with path
% choices_file: choice filename
% sess_range: range of session numbers [lo, hi]
%
%table definitions for t_data and t_all
data_variable_names={'dim','jit_crit','if_bid','neg_pos_bid','ray_num','ray2_num','ray_label','ray2_label','var_name','value_data','value_eblo','value_ebhi','value_sem'};
stats_needed={'data','clo','chi','sem'}; %correspondence between value_[data|eblo|ebhi|sem] and fields of angles_stats, mults_stats
% dim: dimension of model
% jit_crit: rms jitter per coordinate for requested pval (depends on dimension but not on the geometry param)
% if_bid: 0 for a measurement from a (unidirectional) ray, 1 for a bidirectional ray
% neg_pos_bid: m for neg, p for pos, z for bidirectional
neg_pos_bid_text={'neg[-]','pos[+]','bid[-+]'};
%ray_num: ray number (for mult), 
%ray2_no: second ray number (for angle), 0 for mult
%ray_label: ray label (endpoint, eg., g0.40)
%ray2_label: ray lebel for second ray (for angle), '' for mult
%var_name: variable name
%value: value of the variable
%value_eblo, value_ebhi, value_sem: lower and upper confidence limit and
%   standard error of measurement, NaN statistics not calculates
%
underscore='_';
dash='-';
%
%read and concatenate one or more database files, checking for duplicates
%
if_ok=0;
while (if_ok==0)
    ntabs=0;
    t_meta_all_orig=cell(0);
    t_all_orig=cell(0);
    files_orig=cell(0);
    if_more=1;
    while (if_more==1)
        files_req=cell(0);
        n_add_sign=getinp('>0 to add one or more tables by name, <=0 to add one or more from dialog box, 0 if done','d',[-100 100]);
        if (n_add_sign>0) %from console
            n_toadd=abs(n_add_sign);
            files_req=cell(1,n_toadd);
            for iadd=1:abs(n_add_sign)
                files_req{iadd}=getinp(sprintf('table file name, typically %s',ui_filter_gen),'s',[]);
            end
        else %dialog box with logic for specified or unspecified number of files
            [filenames_short,pathname,filter_index]=uigetfile(ui_filter,'select table files','Multiselect','on');            %dialog box here
            if ~iscell(filenames_short) filenames_short={filenames_short}; end
            if (n_add_sign<0) %specified number of files
                if length(filenames_short)~=abs(n_add_sign)
                    n_toadd=0; %will force the dialog box to be ignored
                    disp(sprintf('no files added; %3.0f requested, %3.0f provided',abs(n_add_sign),length(filenames_short)));
                else
                    n_toadd=abs(n_add_sign);
                end
            elseif filter_index==0
                n_toadd=0;
            else
                n_toadd=length(filenames_short);
            end
            for ireq=1:n_toadd
                files_req{ireq}=cat(2,pathname,filesep,filenames_short{ireq});
            end
        end
        for ireq=1:n_toadd
            files_req{ireq}=strrep(files_req{ireq},'\','/');
            if ~exist(files_req{ireq},'file')
                disp(sprintf('%s not found',files_req{ireq}));
            elseif ~isempty(strmatch(files_req{ireq},files_orig,'exact'))
                disp(sprintf('%s already added.',files_req{ireq}));
            else
                s=load(files_req{ireq});
                if ~isfield(s,'t_meta_all') | ~isfield(s,'t_all')
                    disp(sprintf('%s is missing either t_meta_all or t_all',files_req{ireq}));
                else
                    ntabs=ntabs+1;
                    files_orig{ntabs}=files_req{ireq};
                    t_meta_all_orig{ntabs}=s.t_meta_all;
                    t_all_orig{ntabs}=s.t_all;
                    if (ntabs==1)
                        t_meta_all=s.t_meta_all;
                        t_all=s.t_all;
                    else
                        %check that UserData matches (has params for jitter calcs and p-values), but disregard file list if present
                        ifdif_meta=compstruct('current',setfield(t_meta_all.Properties.UserData,'files_orig',[]),'merging',setfield(s.t_meta_all.Properties.UserData,'files_orig',[]));
                        ifdif_all =compstruct('current',setfield(t_all.Properties.UserData,'files_orig',[]),'merging',setfield(s.t_all.Properties.UserData,'files_orig',[]));
                        if ifdif_meta>0
                            disp('Warning: merging a metadata table with different UserData');
                            disp(ifdif_meta);
                            disp('Current:'),
                            disp(t_meta_all.Properties.UserData);
                            disp('Merging:')
                            disp(s.t_meta_all.Properties.UserData);
                        end
                        if ifdif_all>0
                            disp('Warning: merging a data table with different UserData');
                            disp(ifdif_all);
                            disp('Current:'),
                            disp(t_all.Properties.UserData);
                            disp('Merging:')
                            disp(s.t_all.Properties.UserData);
                        end
                        t_meta_all=[t_meta_all;s.t_meta_all];
                        t_all=[t_all;s.t_all];
                    end
                    disp(sprintf('%s added (%4.0f metadata rows, %6.0f data rows), total so far: %4.0f metadata rows, %6.0f data rows.',...
                        files_req{ireq},size(t_meta_all_orig{ntabs},1),size(t_all_orig{ntabs},1),size(t_meta_all,1),size(t_all,1)));
                end %file found?
            end %fields found?
        end %ireq
        disp(sprintf('So far, %3.0f table files concatenated',ntabs));
        if (ntabs>0)
            if (if_replace_avg)
                t_meta_all.subj_model_ID=strrep(strrep(t_meta_all.subj_model_ID,'qform-avg-','qform-'),'qform-','');
                t_all.subj_model_ID=strrep(strrep(t_all.subj_model_ID,'qform-avg-','qform-'),'qform-','');
                disp('replaced qform[-avg]-XX by XX in subj_model_ID');
            end
            %
            %check for duplicates and offer to save the concatenated file
            % note that different file names or different session ranges are not considered,
            %
            underscores=repmat(underscore,size(t_meta_all,1),1);
            concat_meta=cat(2,...
                tablecol2char(t_meta_all,'psy_model'),underscores,...
                tablecol2char(t_meta_all,'subj_model_ID'),underscores,...
                tablecol2char(t_meta_all,'expt_grp'),underscores,...
                tablecol2char(t_meta_all,'expt_name'),underscores,...
                tablecol2char(t_meta_all,'expt_param'),underscores,...
                tablecol2char(t_meta_all,'expt_uid'),underscores);
            nunique_meta=size(unique(concat_meta,'rows'),1);
            ndups_meta=size(t_meta_all,1)-nunique_meta;
            disp(sprintf('%4.0f possibly duplicate metadata rows detected',ndups_meta));
            %
            underscores=repmat(underscore,size(t_all,1),1);
            concat_all=cat(2,...
                tablecol2char(t_all,'psy_model'),underscores,...
                tablecol2char(t_all,'subj_model_ID'),underscores,...
                tablecol2char(t_all,'expt_grp'),underscores,...
                tablecol2char(t_all,'expt_name'),underscores,...
                tablecol2char(t_all,'expt_param'),underscores,...
                tablecol2char(t_all,'expt_uid'),underscores,...
                tablecol2char(t_all,'dim'),underscores,...
                tablecol2char(t_all,'neg_pos_bid'),underscores,...
                tablecol2char(t_all,'ray_label'),underscores,...
                tablecol2char(t_all,'ray2_label'),underscores,...
                tablecol2char(t_all,'var_name'),underscores);
            nunique_all=size(unique(concat_all,'rows'),1);
            ndups_all=size(t_all,1)-nunique_all;
            disp(sprintf('%4.0f possibly duplicate data rows detected',ndups_all));
            %
            if_more=getinp('1 for more files','d',[0 1]);
        end       
    end %if_more
    t_meta_all.Properties.UserData.files_orig=files_orig;
    t_all.Properties.UserData.files_orig=files_orig;
    if_ok=getinp('1 if ok','d',[0 1]);   
end
if getinp('1 to save the cumulative file','d',[0 1],intersect(ntabs-1,0)) %suggest saving if multiple tables combined
    filename_cum=getinp('file name (suggest [raystats|llfits]_cumulative_ddmmmyy.mat)','s',[]);
    save(filename_cum,'t_meta_all','t_all');
end
%
%summarize by criterion_names={'subj_model_ID','expt_grp','expt_uid'}
%
crit_data=struct;
for icrit=1:ncrits
    crit=criterion_names{icrit};
    meta_avail=t_meta_all.(crit);
    data_avail=t_all.(crit);
    crit_data.(crit).values=unique([meta_avail;data_avail]);
    crit_data.(crit).values_short=crit_data.(crit).values;
    for iv=1:length(crit_data.(crit).values)
        if isfield(values_short,crit_data.(crit).values{iv})
            crit_data.(crit).values_short{iv}=values_short.(crit_data.(crit).values{iv});
        end
    end
    disp(sprintf(' for %s',crit));   
    crit_data.(crit).nchoices=length(crit_data.(crit).values);
    for k=1:crit_data.(crit).nchoices %point to all rows in metadata and data that have this value of the criterion      
        crit_data.(crit).meta_pointers{k,1}=strmatch(crit_data.(crit).values{k},t_meta_all.(crit),'exact');
        crit_data.(crit).meta_counts(1,k)=length(crit_data.(crit).meta_pointers{k});
        crit_data.(crit).data_pointers{k,1}=strmatch(crit_data.(crit).values{k},t_all.(crit),'exact');
        crit_data.(crit).data_counts(1,k)=length(crit_data.(crit).data_pointers{k});
        disp(sprintf('     %25s: %5.0f in metadata, %5.0f in data (%s)',...
            crit_data.(crit).values{k},crit_data.(crit).meta_counts(1,k),crit_data.(crit).data_counts(1,k),crit_data.(crit).values_short{k}));
    end   
    disp(sprintf('   total:                      %6.0f             %6.0f',...
        sum(crit_data.(crit).meta_counts),sum(crit_data.(crit).data_counts)))
end
%
if_reselect=1;
plot_type=1;
while (if_reselect==1)
    %
    %select by criterion_names={'subj_model_ID','expt_grp','expt_uid'}
    %
    meta_sel=[1:size(t_meta_all,1)];
    data_sel=[1:size(t_all,1)];
    for icrit=1:ncrits
        crit=criterion_names{icrit};
        crit_data.(crit).meta_pointers_select=[];
        crit_data.(crit).data_pointers_select=[];
        disp(sprintf('available %s',crit));
        for k=1:crit_data.(crit).nchoices
            disp(sprintf(' %2.0f->%s',k,crit_data.(crit).values{k}))
        end
        crit_data.(crit).select=getinp('selection(s)','d',[1 crit_data.(crit).nchoices],[1:crit_data.(crit).nchoices]);
        crit_data.(crit).nselect=length(crit_data.(crit).select);
        crit_data.(crit).data_pointers_select=[];
        for kk=1:crit_data.(crit).nselect
            crit_data.(crit).meta_pointers_select=...
                union(crit_data.(crit).meta_pointers_select,crit_data.(crit).meta_pointers{crit_data.(crit).select(kk)});
            crit_data.(crit).data_pointers_select=...
                union(crit_data.(crit).data_pointers_select,crit_data.(crit).data_pointers{crit_data.(crit).select(kk)});
        end
        meta_sel=intersect(meta_sel,crit_data.(crit).meta_pointers_select);
        data_sel=intersect(data_sel,crit_data.(crit).data_pointers_select);
    end
    %create a table that matches these selections
    if_replot=1;
    while (if_replot)
        t_meta_sel=t_meta_all(meta_sel,:);
        t_data_sel=t_all(data_sel,:);
        nsel=size(t_data_sel,1);
        %add a column that merges var_name ray_label ray2_label neg_pos_bid, with some modifications
        if_ignore_ends=getinp('1 to merge across different values of endpoints of rays','d',[0 1]);
        underscores_dbl=repmat(underscore,nsel,2);
        ray_label=tablecol2char(t_data_sel,'ray_label');
        ray2_label=tablecol2char(t_data_sel,'ray2_label');
        if if_ignore_ends
            ray_label=strvcat(regexprep(cellstr(ray_label),'[0-9].',''));
            ray2_label=strvcat(regexprep(cellstr(ray2_label),'[0-9].',''));
        end
        merged_label=cat(2,tablecol2char(t_data_sel,'var_name'),underscores_dbl,...
           ray_label,underscores_dbl,ray2_label,underscores_dbl,tablecol2char(t_data_sel,'neg_pos_bid'));
        t_merged_label=array2table(cellstr(merged_label),'VariableNames',{'merged_label'});
        t_data_sel=[t_data_sel,t_merged_label];
        %
        for k=1:nplot_types
            disp(sprintf('%1.0f ->%s',k,plot_types{k}));
        end
        %
        plot_type=getinp('plot type','d',[1 nplot_types],plot_type);
        %
        for icrit=1:ncrits
            disp(sprintf('%1.0f->%s',icrit,criterion_names{icrit}));
        end
        switch plot_type
            case {1,2}
            % 1-> 'plot a geometric parameter (rows) for a specific model dimension (columns)',...
            % 2-> 'plot a geometric parameter (columns) for a specific model dimension (rows)',...
                key_strings={'major grouping','minor grouping'};
                nkeys=3;
                nkeys_ask=2;
                dim_range_min=1;
            case {3,4}
            % 3-> 'plot a geometric parameter (rows) as function of dimension'};
            % 4-> 'plot a geometric parameter (columns) as function of dimension'};
                key_strings={'column','within plot'};
                nkeys=3;
                nkeys_ask=2;
                dim_range_min=2;
        end
        %get distinct values for crit_key(1:nkeys_ask)
        crit_key=zeros(1,nkeys);
        for ikey=1:nkeys_ask
            while any(crit_key(ikey)==crit_key(1:ikey-1)) | crit_key(ikey)==0
                crit_key(ikey)=getinp(cat(2,'choice for ',key_strings{ikey}),'d',[1 ncrits]);
            end
        end
        %fill in crit_key(nkeys_ask+1:nkeys)
        for ikey=nkeys_ask+1:nkeys
            crit_key(ikey)=min(setdiff(1:ncrits,crit_key(1:ikey-1)));
        end
        %
        dims_avail=unique(cell2mat(t_data_sel{:,'dim'}));
        dims_sel=[];
        while length(dims_sel)<dim_range_min
            dims_sel=getinp(sprintf('dimension(s) of models to use, at least %2.0f values',dim_range_min),'d',[0 max(dims_avail)]);
        end
        crit_data_structs=cell(nkeys,1);
        for ikey=1:nkeys
            crit_data_structs{ikey}=crit_data.(criterion_names{crit_key(ikey)});
        end
        %
        vars_avail=unique(t_merged_label.merged_label);
        for k=1:length(vars_avail)
            disp(sprintf('%2.0f->%s',k,strrep(vars_avail{k},'__',' ')));
        end
        vars_sel=getinp('choice','d',[1 length(vars_avail)]);
        %
        if_eb=getinp('1 to plot error bars','d',[0 1]);
        if ~isempty(cell2mat(strfind(vars_avail(vars_sel),'cosang')))
            if_angle=getinp('1 to transform from cos to ang','d',[0 1]);
        else
            if_angle=0;
        end
        switch plot_type
            case {1,2}
            % 1-> 'plot a geometric parameter (rows) for a specific model dimension (columns)',...
            % 2-> 'plot a geometric parameter (columns) for a specific model dimension (rows)',...
                figure;
                set(gcf,'Position',[100 100 1200 800]);               
                nmajor=crit_data_structs{1}.nselect;
                nminor=crit_data_structs{2}.nselect;
                tick_posits=zeros(nminor,nmajor);
                tick_labels=cell(nminor,1);
                for iminor=1:nminor
                    tick_labels{iminor}=crit_data_structs{2}.values_short{crit_data_structs{2}.select(iminor)};
                    tick_labels{iminor}=strrep(tick_labels{iminor},underscore,' ');
                    for imajor=1:nmajor
                        tick_posits(iminor,imajor)=(imajor-1)*(plot_major_space+(nminor-1)*plot_minor_space)+(iminor-(nminor+1)/2)*plot_minor_space;
                    end
                end
                xlims=0.5*plot_major_space*[-1 1]+[min(tick_posits(:)) max(tick_posits(:))];
                %create separate tables for each subplot
                t_subplot=cell(length(dims_sel),length(vars_sel));
                for idim_sel=1:length(dims_sel)
                    dim_sel=dims_sel(idim_sel);
                    for ivar_sel=1:length(vars_sel)
                        var_sel=vars_avail{vars_sel(ivar_sel)};
                        t_subplot{idim_sel,ivar_sel}=t_data_sel(strmatch(var_sel,t_data_sel.merged_label,'exact'),:);
                        t_subplot{idim_sel,ivar_sel}=t_subplot{idim_sel,ivar_sel}(find(cell2mat(t_subplot{idim_sel,ivar_sel}.dim)==dim_sel),:);
                        if (plot_type==1)
                            subplot(length(vars_sel),length(dims_sel),idim_sel+(ivar_sel-1)*length(dims_sel));
                        else
                            subplot(length(dims_sel),length(vars_sel),ivar_sel+(idim_sel-1)*length(vars_sel));
                        end
                        %go through the table and plot
                        hl=cell(0);
                        ht=[];
                        t_plot=t_subplot{idim_sel,ivar_sel};
                        for k=1:size(t_plot,1)
                             imajgrp=strmatch(t_plot{k,criterion_names{crit_key(1)}},crit_data_structs{1}.values(crit_data_structs{1}.select),'exact');
                             imingrp=strmatch(t_plot{k,criterion_names{crit_key(2)}},crit_data_structs{2}.values(crit_data_structs{2}.select),'exact');
                             iwithin=strmatch(t_plot{k,criterion_names{crit_key(3)}},crit_data_structs{3}.values(crit_data_structs{3}.select),'exact');
                             % [k imajgrp imingrp iwithin]
                             values_plot=zeros(1,4);
                             if ~isempty(imajgrp) & ~isempty(imingrp) & ~isempty(iwithin)
                                 values_plot(1)=cell2mat(t_plot{k,'value_data'});
                                 values_plot(2)=cell2mat(t_plot{k,'value_eblo'});
                                 values_plot(3)=cell2mat(t_plot{k,'value_ebhi'});
                                 values_plot(4)=cell2mat(t_plot{k,'value_sem'}); 
                                 if ~isempty(strfind(vars_avail{vars_sel(ivar_sel)},'cosang')) &  if_angle==1
                                     sem_ratio=values_plot(4)/(values_plot(3)-values_plot(2)); %to fix sem after cosine transform
                                     values_plot(1:3)=acos(min(max(values_plot(1:3),-1),1))*180/pi;
                                     values_plot(4)=sem_ratio*(values_plot(3)-values_plot(2));
                                 end
                                 hp=plot(tick_posits(imingrp,imajgrp),values_plot(1),'k.');                                     
                                 hold on;
                                 subj_model_ID=cell2mat(t_plot{k,'subj_model_ID'});                                   
                                 expt_grp=cell2mat(t_plot{k,'expt_grp'});                                   
                                 expt_uid=cell2mat(t_plot{k,'expt_uid'});
                                 %
                                 psg_raystats_dbplot_style; %set plot style, uses hp subj_model_ID,expt_grp, expt_uid
                                 if isempty(strmatch(upper(subj_model_ID),ht,'exact'))
                                     ht=strvcat(ht,upper(subj_model_ID));
                                     hl=[hl;hp];
                                 end
                                 if if_eb
                                     heb=plot(tick_posits(imingrp,imajgrp)+plot_ebhw*[-1 1],repmat(values_plot(2),1,2),'k');
                                     set(heb,'Color',style.color);
                                     heb=plot(tick_posits(imingrp,imajgrp)+plot_ebhw*[-1 1],repmat(values_plot(3),1,2),'k');
                                     set(heb,'Color',style.color);
                                     heb=plot(repmat(tick_posits(imingrp,imajgrp),1,2),values_plot(2:3),'k');
                                     set(heb,'Color',style.color);
                                 end
                             end
                        end
                        if ~isempty(ht)
                            legend(hl,ht,'Location','Best');
                        end
                        title_string=strrep(strrep(vars_avail{vars_sel(ivar_sel)},'__',' '),'[]','');
                        if ~isempty(strfind(vars_avail{vars_sel(ivar_sel)},'cosang')) & if_angle==1
                            title_string=strrep(title_string,'cosang','angle');
                        end
                        title(sprintf(' %s from dim %1.0f',title_string,dims_sel(idim_sel)),'Interpreter','none');
                        %
                        psg_raystats_dbplot_ticks; %uses tick_posits,tick_labels,xlims,vars_avail,vars_sel,ivar_sel,if_angle
                        %
                        for imajor=1:nmajor %additional labels
                            maj_label=crit_data_structs{1}.values_short{crit_data_structs{1}.select(imajor)};
                            maj_label=strrep(maj_label,underscore,' ');
                            text(tick_posits(1,imajor),plot_label_height*max(get(gca,'YLim')),maj_label,'FontSize',7);
                        end
                        %
                    end %ivar_sel
                end %idim_sel
            case {3,4}
            % 3-> 'plot a geometric parameter (rows) as function of dimension'};
            % 4-> 'plot a geometric parameter (columns) as function of dimension'};
            %row or column is then used for crit_key(1), and each value of
            %crit_key(2) is shown within plot, crit_key(3) varied by page
                tick_posits=dims_sel;
                tick_labels=tick_posits;
                xlims=[-0.5 0.5]+[min(dims_sel),max(dims_sel)];
                %create separate tables for each subplot
                t_subplot=cell(crit_data_structs{1}.nselect,length(vars_sel),crit_data_structs{3}.nselect);
                for ipage=1:size(t_subplot,3)
                    crit_page=crit_data_structs{3}.values{crit_data_structs{3}.select(ipage)};
                    crit_page_short=crit_data_structs{3}.values_short{crit_data_structs{3}.select(ipage)};
                    figure;
                    set(gcf,'Position',[100 100 1200 800]); 
                    set(gcf,'NumberTitle','off');
                    set(gcf,'Name',crit_page_short);
                    for icrit_sel=1:size(t_subplot,1)
                        crit_val=crit_data_structs{1}.values{crit_data_structs{1}.select(icrit_sel)};
                        crit_val_short=crit_data_structs{1}.values_short{crit_data_structs{1}.select(icrit_sel)};
                        for ivar_sel=1:length(vars_sel)
                            var_sel=vars_avail{vars_sel(ivar_sel)};
                            %
                            t_subplot{icrit_sel,ivar_sel,ipage}=t_data_sel(strmatch(var_sel,t_data_sel.merged_label,'exact'),:); %select for this row and column
                            t_subplot{icrit_sel,ivar_sel,ipage}=t_subplot{icrit_sel,ivar_sel,ipage}(strmatch(crit_page,t_subplot{icrit_sel,ivar_sel,ipage}.(criterion_names{crit_key(3)}),'exact'),:);
                            t_subplot{icrit_sel,ivar_sel,ipage}=t_subplot{icrit_sel,ivar_sel,ipage}(strmatch(crit_val,t_subplot{icrit_sel,ivar_sel,ipage}.(criterion_names{crit_key(1)}),'exact'),:);
                            %
                            t_plot=t_subplot{icrit_sel,ivar_sel,ipage};
                            if (plot_type==3)
                                subplot(length(vars_sel),size(t_subplot,1),icrit_sel+(ivar_sel-1)*size(t_subplot,1));
                            else
                                subplot(size(t_subplot,1),length(vars_sel),ivar_sel+(icrit_sel-1)*length(vars_sel));
                            end
                            %go through the table and plot, based on crit_data_structs{2,3}
                            hl=cell(0);
                            ht=[];
                            %
                            t_plot=t_subplot{icrit_sel,ivar_sel,ipage};
                            for k=crit_data_structs{2}.select %plot each line graph
                                plot_label=crit_data_structs{2}.values_short{k};
                                plot_rows=strmatch(crit_data_structs{2}.values(k),t_plot{:,criterion_names{crit_key(2)}},'exact');
                                dims=cell2mat(t_plot{plot_rows,'dim'});
                                dims_selptrs=find(ismember(dims,dims_sel)); %which dimensions were chosen to plot
                                if ~isempty(dims_selptrs)
                                    dims_plot=dims(dims_selptrs);
                                    values_plot=zeros(length(dims_selptrs),4);
                                    values_plot(:,1)=cell2mat(t_plot{plot_rows(dims_selptrs),'value_data'});
                                    values_plot(:,2)=cell2mat(t_plot{plot_rows(dims_selptrs),'value_eblo'});
                                    values_plot(:,3)=cell2mat(t_plot{plot_rows(dims_selptrs),'value_ebhi'});
                                    values_plot(:,4)=cell2mat(t_plot{plot_rows(dims_selptrs),'value_sem'});
                                    if ~isempty(strfind(vars_avail{vars_sel(ivar_sel)},'cosang')) &  if_angle==1
                                        sem_ratio=values_plot(:,3)-values_plot(:,2);
                                        values_plot(:,1:3)=acos(min(max(values_plot(:,1:3),-1),1))*180/pi;
                                        values_plot(:,4)=sem_ratio.*(values_plot(:,3)-values_plot(:,2));
                                    end
                                    hp=plot(dims_plot,values_plot(:,1),'k');                                     
                                    hold on;
                                        %check that only one subject and one expt_grp are plotted
                                    subj_model_ID=unique(cell2mat(t_plot{plot_rows(dims_selptrs),'subj_model_ID'}),'rows');
                                    expt_grp=unique(cell2mat(t_plot{plot_rows(dims_selptrs),'expt_grp'}),'rows');
                                    expt_uid=unique(cell2mat(t_plot{plot_rows(dims_selptrs),'expt_uid'}),'rows');
                                    %
                                    psg_raystats_dbplot_style; %set plot style, uses hp subj_model_ID,expt_grp, expt_uid
                                    if if_eb
                                        for kd=1:size(values_plot,1);
                                            heb=plot(dims_plot(kd)+plot_ebhw*[-1 1],repmat(values_plot(kd,2),1,2),'k');
                                            set(heb,'Color',style.color);
                                            heb=plot(dims_plot(kd)+plot_ebhw*[-1 1],repmat(values_plot(kd,3),1,2),'k');
                                            set(heb,'Color',style.color);
                                            heb=plot(repmat(dims_plot(kd),1,2),values_plot(kd,[2 3]),'k');
                                            set(heb,'Color',style.color);
                                        end
                                    end
                                    %
                                    if isempty(strmatch(upper(plot_label),ht,'exact'))
                                        ht=strvcat(ht,upper(plot_label));
                                        hl=[hl;hp];
                                    end
                                end
                            end 
                            %
                            psg_raystats_dbplot_ticks; %common code for tickmarks and labels
                            xlabel('dim');
                            if ~isempty(ht)
                                legend(hl,ht,'Location','Best');
                            end
                            title_string=strrep(strrep(vars_avail{vars_sel(ivar_sel)},'__',' '),'[]','');
                            if ~isempty(strfind(vars_avail{vars_sel(ivar_sel)},'cosang')) & if_angle==1
                                title_string=strrep(title_string,'cosang','angle');
                            end
                            title(sprintf(' %s from %s (%s)',title_string,crit_val_short,crit_page_short),'Interpreter','none');
                            %
                        end %icrit_sel
                    end %ivar_sel
                end %ipage
         end  %which plot_type
        if_replot=getinp('1 to replot with these selections','d',[0 1]);
    end
    if_reselect=getinp('1 to choose another selection','d',[0 1]);
end


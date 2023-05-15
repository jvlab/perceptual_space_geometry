%psg_choicedata_merge: a script to read one or more choice files from psg experiments,
% check compatibility, and merge them.
%
% if stim_list match but in wrong order, this is fixed.
% all triads are reorganized so that s1 tag < s2 tag
%
% logic borrowed from figgnd_proc.
%
%   See also: PSG_CHOICEDATA_READ, PSG_UMI_TRIPLIKE_DEMO, PSG_TENTLIKE_DEMO.
%
fields_needed={'stim_list','responses','responses_colnames'};
fields_match={'stim_list','responses_colnames'};
col_trials=5;
col_closer=4;
%
if_dialog=double(~exist('names'));
if ~exist('ui_string') ui_string='*choices*.mat'; end
if_dialog=getinp('1 to specify file names from dialog box (0 if already specified)','d',[0 1],if_dialog);
if (if_dialog)
    if_ok=0;
    while (if_ok==0)
        [fns_raw,pathname,fptrs]=uigetfile(ui_string,'Select files to merge and analyze','Multiselect','on');       
        if iscell(fns_raw)
            fns=fns_raw;
        else
            fns=cell(0);
            fns{1}=fns_raw;
        end
        if fptrs~=0
            disp(sprintf(' specified path:',pathname));
            names=cell(0);
            for ifn=1:length(fns)
                names{ifn}=cat(2,pathname,fns{ifn});
                disp(sprintf(' file %3.0f->%s',ifn,fns{ifn}));
            end
            if_ok=getinp('1 if ok','d',[0 1]);
        end
    end
end
%
nsets=length(names);
s=cell(0);
names_short=cell(0);
%
for iset=1:nsets %first determine unique rays
    s{iset}=load(names{iset});
    name_short=names{iset};
    filesep_loc=union(strfind(name_short,'/'),strfind(name_short,'\'));
    if ~isempty(filesep_loc)
        name_short=name_short(max(filesep_loc)+1:end);
    end
    names_short{iset}=name_short;
    disp(sprintf(' file %3.0f: %40s (short: %40s) read',iset,names{iset},names_short{iset}));
end
merge_ok=ones(1,nsets);
stim_list=s{1}.stim_list;
responses_colnames=s{1}.responses_colnames;
nstims=size(stim_list,1);
set1stim=zeros(nstims,nsets);
resps=cell(nsets,1); %resps always has s1 tag < s2 tag
for iset=1:nsets
    pstring=sprintf('%40s:',names_short{iset}');
    for ishow=1:length(fields_needed)
        fs=fields_needed{ishow};
        if isfield(s{iset},fs)
            pstring=cat(2,pstring,sprintf('%20s: [%6.0f %2.0f]  ',fs,size(s{iset}.(fs),1),size(s{iset}.(fs),2)));
        else
            pstring=cat(2,pstring,sprintf('%20s:   missing    ',fs));
            merge_ok(iset)=0;
        end
    end
    disp(pstring);
    if merge_ok(iset)==1
          %check consistency
        for imatch=1:length(fields_match);
            fm=fields_match{imatch};
            if any(size(s{1}.(fm))~=size(s{iset}.(fm)))
                disp(sprintf('   mismatch in size of %s',fm));
                merge_ok(iset)=0;
            end
        end
        if merge_ok(iset)==1
            if ~strcmp(responses_colnames,s{iset}.responses_colnames)
                disp(sprintf('  mismatch in %s','responses_colnames'));
                merge_ok(iset)=0;
            end
            %create a sort list
            nstims=size(stim_list,1);
            for istim=1:nstims
                ind=strmatch(s{iset}.stim_list(istim,:),stim_list,'exact');
                if length(ind)==1
                    set1stim(istim,iset)=ind;
                end
            end
            if any(set1stim(:,iset)==0)
                disp(sprintf(' mismatch in %s','stim_list'));
                merge_ok(iset)=0;
            elseif any(set1stim(:,iset)~=[1:nstims]')
                disp(' stim_list in a different order but reordered to match first dataset');
                s{iset}.responses(:,[1:col_closer-1])=reshape(set1stim(s{iset}.responses(:,[1:col_closer-1]),iset),[size(s{iset}.responses,1),col_closer-1]);
            end
            %
        end
    end
end
sets_merge=find(merge_ok==1);
sets_merge=intersect(sets_merge,getinp('sets to merge','d',[min(sets_merge) max(sets_merge)],sets_merge));
resps_all=zeros(0,col_trials);
for imerge=1:length(sets_merge)
    iset=sets_merge(imerge);
    %
    %reorganize the triads so that the tag for comparison stim 1 is
    %less than for comparison stim 2, to ensure that equivalent triads will match
    %
    r=s{iset}.responses;
    need_flip=find(r(:,2)>r(:,3));
    r(need_flip,col_closer)=r(need_flip,col_trials)-r(need_flip,col_closer);
    r(need_flip,[2 3])=r(need_flip,[3 2]);
    resps{iset}=r;
    resps_all=[resps_all;resps{iset}];
    disp(sprintf(' merging set %4.0f: %7.0f trials, %6.0f unique triads, yielding %7.0f trials, %6.0f unique triads',iset,...
        sum(resps{iset}(:,col_trials)),size(unique(resps{iset}(:,1:col_closer-1),'rows'),1),...
        sum(resps_all(:,col_trials)),size(unique(resps_all(:,1:col_closer-1),'rows'),1)));
end
resps_sorted=sortrows(resps_all);
%find the rows where at least one of the first 3 columns differs from the next
newrows=[find(any(resps_sorted(1:end-1,1:col_closer-1)~=resps_sorted(2:end,1:col_closer-1),2))];
newrows=[newrows;size(resps_sorted,1)];
triads_merged=resps_sorted(newrows,[1:col_closer-1]); %triads are the unqiue names
%compute sums from cumulative sums before and after
cumsums=cumsum(resps_sorted(:,[col_closer col_trials]));
data_merged=[cumsums(newrows(1),:);cumsums(newrows(2:end),:)-cumsums(newrows(1:end-1),:)];
responses=[triads_merged,data_merged];
%
fn_merged=getinp('merged file name (and path), *choices*.mat','s',[]);
save(fn_merged,'responses','stim_list','responses_colnames');

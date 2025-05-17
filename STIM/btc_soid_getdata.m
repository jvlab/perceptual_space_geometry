function [ed,eb_avail,ds,condno,ou]=btc_soid_getdata(pdata_fn,opts)
% [ed,eb_avail,ds,condno,ou]=btc_soid_getdata(pdata_fn,opts) shows the psychophysical 
% dataset structures in a file, and chooses one
%
% pdata_fn:  file name (requested if not supplied), containing one or more datasets (ds's).
% Each ds contains one or more conditions, each condition is an "edirs" structure
%     low and high error bars cah be returned either as edirs.yx.thresh_mags_[eblo|ebhi]
%     or as second and third columns of edirs.yx.thresh_mags
%     Note that a value of eblo or ebhi of 0, OR, eblo or ebhi=thresh_mags, indicates
%       that this error bar is invalid, and must be filled in
%
% opts.if_choose: 1 to select a file (0: just show the file)
% opts.eb_fill: how to fill in exceptional error bars, defaults to 1
%     0->do not fill
%     n->fill with a fraction of the threshold that is equal to the mean of
%     the n largest fractions (Inf means use all valid error bars in that plane)
%
% ed: an edirs structure, has edirs fields: ed.yx.thresh_mags and optionally thresh_mags_[eblo|ebhi]
% eb_avail: 1 if all error bars are present, 0 otherwise
% ds: the data structure containing ed.  Note that one or more ds's can be in a single file.
%  it has the following fields:
%  ds.subject
%  ds.desc
%  ds.condition{*}, a cell array for each background (or other) condition
%  ds.condition{*}.figgr (0, 1, or 2)
%  ds.condition{*}.desc
%  ds.condition{*}.edirs (has edirs fields: ds.condition{*}.edirs.yx.thresh_mags, returned as ed)
% condno:  ed=ds.condition{condno}.edirs
% ou: options used
%
% 05Mar21:  add option to supply for dsno (dataset number) and condno  (typically  condition  1: bkgd structured,2: bkgd random, 3: combined)
%
%    See also:  BTC_SOID_FIT, BTC_SOID_TEST, BTC_SOID_DEMO, BTC_EDIRS, BTC_SOID_EBFIX, BTC_SOID_ARCH2DS.
%
if (nargin<=1) opts=[]; end
opts=filldefault(opts,'if_choose',1);
opts=filldefault(opts,'eb_fill',1);
opts=filldefault(opts,'dsno',[]);
opts=filldefault(opts,'condno',[]);
%
ou=opts;
if isempty(pdata_fn)
    pdata_fn=getinp('name of file with psychophysical data','s',[],pdata_fn);
end
%
ed=[];
eb_avail=0;
ds=[];
condno=0;
s=load(pdata_fn);
dsets=fieldnames(s);
%
disp(' datasets            name  #conds     subject                 descriptor')
for id=1:length(dsets)
    nconds(id)=0;
    if (isfield(s.(dsets{id}),'subject')) & (isfield(s.(dsets{id}),'desc')) & (isfield(s.(dsets{id}),'condition'))
        if iscell(s.(dsets{id}).condition)
            nconds(id)=length(s.(dsets{id}).condition);
        end
        for icond=1:nconds(id)
            nplanes(id,icond)=0;
            c=s.(dsets{id}).condition{icond};
            if (isfield(c,'figgr')) & (isfield(c,'desc')) & (isfield(c,'edirs'))
                if (isstruct(c.edirs))
                    nplanes(id,icond)=length(fieldnames(c.edirs));
                end
            else
                nconds(id)=0;
            end
        end
    end
    if (nconds(id)>0)
        disp(sprintf(' %2.0f %20s %5.0f %15s %30s',id,dsets{id},nconds(id),s.(dsets{id}).subject,s.(dsets{id}).desc))
        for icond=1:nconds(id)
            c=s.(dsets{id}).condition{icond};
            disp(sprintf('     condition %2.0f: nplanes=%4.0f, figgr=%2.0f, descriptor=%s',...
                icond,nplanes(id,icond),c.figgr,c.desc));
        end
    else
        disp(sprintf('[%2.0f] %20s **missing fields**',id,dsets{id}));
    end
end
if (opts.if_choose>0)
    ifok=0;
    while (ifok==0)
        if (isempty(opts.dsno))
            dsno=getinp('choice of dataset to analyze','d',[1 length(dsets)]);
        else
            dsno=opts.dsno;
        end
        if (nconds(dsno)>0)
            ifok=1;
        else
            disp('That dataset has missing fields.');
        end
    end
    ds=s.(dsets{dsno});
else
    return
end
%
% show the planes in each condition of that dataset
%
fns=[];
for icond=1:nconds(dsno)
    c=s.(dsets{dsno}).condition{icond};
    fns=[fns;fieldnames(c.edirs)];
end
fns=unique(fns);
nplanes_all=length(fns);
ptable=zeros(nplanes_all,nconds(dsno));
cell(nplanes_all,nconds(dsno));
disp('table of planes and conditions, "-" means that no error bars are available')
disp(sprintf(cat(2,'    plane/condition    ',repmat(' %5.0f',1,nconds(dsno))),[1:nconds(dsno)]));
for iplane=1:nplanes_all
    for icond=1:nconds(dsno)
        c=s.(dsets{dsno}).condition{icond};
        if isfield(c.edirs,fns{iplane})
            ep=c.edirs.(fns{iplane});
            if isfield(ep,'thresh_mags')
                ptable(iplane,icond)=size(ep.thresh_mags,1);
            end
            % determine if error bars are present and have the correct length
            ebsign=-1;
            if isfield(ep,'thresh_mags_eblo') & isfield(ep,'thresh_mags_ebhi')
                if (size(ep.thresh_mags_eblo,1)==size(ep.thresh_mags,1)) & ...
                        (size(ep.thresh_mags_ebhi,1)==size(ep.thresh_mags,1))
                    ebsign=1;
                    ep.thresh_mags=ep.thresh_mags(:,1); %make sure it no longer has 3 columns
                end
            elseif size(ep.thresh_mags,2)==3 %if three columns, assume that the cols represent val, eblo, ebhi
                ep.thresh_mags_eblo=ep.thresh_mags(:,2);
                ep.thresh_mags_ebhi=ep.thresh_mags(:,3);
                ep.thresh_mags=ep.thresh_mags(:,1);
                ebsign=1;
            end
            c.edirs.(fns{iplane})=ep; %in case edirs was modified b/o 3-colum format for edirs.thresh_mags
            ptable(iplane,icond)=ptable(iplane,icond)*ebsign;
        end
        s.(dsets{dsno}).condition{icond}=c;
    end
end    
for iplane=1:nplanes_all
    disp(sprintf(cat(2,'     %-18s',repmat(' %5.0f',1,nconds(dsno))),fns{iplane},ptable(iplane,:)));
end
if isempty(opts.condno)
    condno=getinp('choice of condition to analyze','d',[1 nconds(dsno)]);
else
    condno=opts.condno;
end
ed=s.(dsets{dsno}).condition{condno}.edirs;
eb_avail=1;
if any(ptable(:,condno)<0)
    eb_avail=0;
else
    [ed,eb_avail,eb_msg]=btc_soid_ebfix(ed,opts);
    if (eb_avail==0)
        disp('Warning: error bars do  not pass sanity check.');
        disp(eb_msg);
    else
        disp(eb_msg);
    end
end
return

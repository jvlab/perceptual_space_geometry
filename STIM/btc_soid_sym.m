function [esym,ou]=btc_soid_sym(edata,symopts)
% [esym,ou]=btc_soid_sym(edata,symopts) symmetrizes a btc psychophysical
% dataset, and also adds thresh_vecs and thresh_vecs_[eblo|ebhi]
%
%  input
%
%  edata: a btc psychophysical dataset, each field is a plane
%  symopts: all options default to 0
%    symopts.if_symm=1 to symmetrize around origin
%    symopts.if_axes=1 to uniformize data along every occurrence on 
%      an axis (so that the values are identical in all planes)
%      (but thresholds in positive and negative directions may be different)
%      This eliminates differences between the smae condition measured on different
%      days.
%    symopts.if_axes=-1 to keep different directions separate (i.e., not to 
%      uniformize within bc, de, tuvw
%    symopts.tol: tolerance for matching axes
%
%  notes:
% if_symm operates on diagonals *and* axes
% if_axes does not uniformize diagonals, e.g., force tu=tv (would not be appropriate)
%     this is necessary; if_axes combines t,u,v,w (e.g.) by summing, so it couldn't tell
%     the difference between tu and tv; if_axes just finds the antipodes by
%     adding ndirs/2 to each direction index
% thresholds are combined by averaging
% error bars are combined by rms around their respective thresholds
%
%  output
%
%  esym: output structure
%  ou: options used
%    ou.onaxis_list: a list of pointers to all_coord_table that are on-axis
%    ou.all_coord_table:  table of all coordinates
%    ou.onaxis_coord_table: table of on-axis coordinates
%    ou.onaxis_coord_table_uniquerows: the unique rows of the above, after if_axes uniformizing
%    ou.onaxis_coord_table_indices: which of the unique rows correspond to original onaxis_coord_table
%
%   See also:  BTC_SOID_DEMO, BTC_SOIDFG_SYMMETRIZE.
%
if (nargin<=1) symopts=[]; end
symopts=filldefault(symopts,'if_symm',0);
symopts=filldefault(symopts,'if_axes',0);
symopts=filldefault(symopts,'tol',0.0001);
ou=symopts;
%
dict=btc_define;
nbtc=length(dict.codel); %10
%
planes=fieldnames(edata);
nplanes=length(planes);
%
plane_table=[];
dir_table=[];
onaxis_list=[];
onaxis_coord_table=[];
all_coord_table=[];
plane_off=0;
%first accumulate a table of in-plane directions
for iplane=1:nplanes
    pn=planes{iplane};
    uvecs=edata.(pn).uvecs;
    %pull out all unit vectors that are, within tol, along an axis but
    %not at the origin
    onaxis=find(mod(sum(double(abs(uvecs)<symopts.tol),2),2)==1); %the mod(:,2)==1 does an xor
    all_coords=zeros(size(uvecs,1),nbtc);
    all_coords(:,find(dict.codel==pn(1)))=uvecs(:,2); % the second coord of the plane name is in the first column of uvecs
    all_coords(:,find(dict.codel==pn(2)))=uvecs(:,1); % the first coord of the plane name is in the second column of uvecs
    onaxis_coords=all_coords(onaxis,:);
    all_coord_table=[all_coord_table;all_coords];
    onaxis_coord_table=[onaxis_coord_table;onaxis_coords];
    plane_table=[plane_table;[repmat(iplane,length(onaxis),1)]];
    dir_table=[dir_table;onaxis];
    onaxis_list=[onaxis_list;plane_off+onaxis];
    plane_off=plane_off+size(uvecs,1);
end
%take care of uniformizing within [bc] [de] [tuvw]
if (symopts.if_axes==1)
    z=dict.name_order_aug;
    [uz,iz,jz]=unique(z);
    zmat=zeros(length(uz),length(z));
    for k=1:length(uz);
        zm=strmatch(z{iz(k)},z');
        zmat(min(strmatch(z{iz(k)},z')),zm)=1;
    end
    ou.zmat=zmat;
    onaxis_coord_table=onaxis_coord_table*zmat'; %this combines coords of same type (i.e., bc, de, tuvw)
end
%discretize at precision of symopts.tol
onaxis_coord_table=symopts.tol*round(onaxis_coord_table/symopts.tol);
ou.onaxis_list=onaxis_list;
ou.all_coord_table=all_coord_table;
ou.onaxis_coord_table=onaxis_coord_table;
ou.plane_table=plane_table;
ou.dir_table=dir_table;
ou.onaxis_coord_table_uniquerows=ou.onaxis_coord_table;
ou.onaxis_coord_table_indices=[1:size(ou.onaxis_coord_table,1)]';
ou.onaxis_coord_table_ptrs=[1:size(ou.onaxis_coord_table,1)]';
%
if ~(symopts.if_axes==0)
    [ur,ir,jr]=unique(ou.onaxis_coord_table,'rows');
    ou.onaxis_coord_table_uniquerows=ur;
    ou.onaxis_coord_table_indices=jr;
    ou.onaxis_coord_table_ptrs=ir;
    for iu=1:max(jr)
        matches=find(jr==iu);
        sum_mags=0;
        ss_mags_eblo=0;
        ss_mags_ebhi=0;
        eb_count=0;
        %calculate the averages across repeated measures
        for im=1:length(matches)
            pn=planes{plane_table(matches(im))};
            di=dir_table(matches(im));
            sum_mags=sum_mags+edata.(pn).thresh_mags(di);
            if isfield(edata.(pn),'thresh_mags_eblo') & isfield(edata.(pn),'thresh_mags_ebhi')
                eb_count=eb_count+1;
                ss_mags_eblo=ss_mags_eblo+(edata.(pn).thresh_mags_eblo(di)-edata.(pn).thresh_mags(di)).^2; %combine error bars by sum of squares
                ss_mags_ebhi=ss_mags_eblo+(edata.(pn).thresh_mags_ebhi(di)-edata.(pn).thresh_mags(di)).^2; %combine error bars by sum of squares
            end
        end
        av_mags=sum_mags/length(matches);
        %install the averages across repeated measures
        if (eb_count>0)
            av_mags_eblo=av_mags-sqrt(ss_mags_eblo/eb_count);
            av_mags_ebhi=av_mags+sqrt(ss_mags_ebhi/eb_count);
        end
        for im=1:length(matches)
            pn=planes{plane_table(matches(im))};
            di=dir_table(matches(im));
            edata.(pn).thresh_mags(di)=av_mags;
            edata.(pn).thresh_vecs(di,:)=av_mags*edata.(pn).uvecs(di,:);
            if isfield(edata.(pn),'thresh_mags_eblo') & isfield(edata.(pn),'thresh_mags_ebhi')
                edata.(pn).thresh_mags_eblo(di)=av_mags_eblo;
                edata.(pn).thresh_mags_ebhi(di)=av_mags_ebhi;
            end
            if isfield(edata.(pn),'thresh_vecs_eblo') & isfield(edata.(pn),'thresh_vecs_ebhi')
                edata.(pn).thresh_vecs_eblo(di,:)=av_mags_eblo*edata.(pn).uvecs(di,:);
                edata.(pn).thresh_vecs_ebhi(di,:)=av_mags_ebhi*edata.(pn).uvecs(di,:);
            end
        end
    end
end
%
esym=[];
for iplane=1:nplanes
    pn=planes{iplane};
    esym.(pn)=edata.(pn);
    eref=esym.(pn); %keep for reference
    if (symopts.if_symm==1)
        ndirs=length(edata.(pn).thresh_mags);
        for ptr=1:ndirs
            ptr_opp=1+mod(ptr-1+ndirs/2,ndirs);
            esym.(pn).thresh_mags(ptr)=(eref.thresh_mags(ptr)+eref.thresh_mags(ptr_opp))/2;
            if isfield(edata.(pn),'thresh_mags_eblo') & isfield(edata.(pn),'thresh_mags_ebhi')
                ss_mags_eblo=sum((eref.thresh_mags_eblo([ptr ptr_opp])-eref.thresh_mags([ptr ptr_opp])).^2); %combine error bars by sum of squares
                ss_mags_ebhi=sum((eref.thresh_mags_ebhi([ptr ptr_opp])-eref.thresh_mags([ptr ptr_opp])).^2); %combine error bars by sum of squares
                esym.(pn).thresh_mags_eblo(ptr)=esym.(pn).thresh_mags(ptr)-sqrt(ss_mags_eblo/2);
                esym.(pn).thresh_mags_ebhi(ptr)=esym.(pn).thresh_mags(ptr)+sqrt(ss_mags_eblo/2);
            end
        end
        esym.(pn).thresh_vecs=repmat(esym.(pn).thresh_mags,1,2).*esym.(pn).uvecs;
        if isfield(edata.(pn),'thresh_vecs_eblo') & isfield(edata.(pn),'thresh_vecs_ebhi')
            esym.(pn).thresh_vecs_eblo=repmat(esym.(pn).thresh_mags_eblo,1,2).*esym.(pn).uvecs;
            esym.(pn).thresh_vecs_ebhi=repmat(esym.(pn).thresh_mags_ebhi,1,2).*esym.(pn).uvecs;
        end
    end
end
return;

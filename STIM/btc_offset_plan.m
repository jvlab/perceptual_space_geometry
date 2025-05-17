% btc_offset_plan:  planning measurements of thresholds offset from the origin
% based on btc_soid_demo
%
%  modifications Nov 30 2011 to step through multiple fitting options
%  and save a cell array r{k}, with r{k}.results and other fields
%
%  modifications Nov 15 2012 to allow for fitting of models even if some thresholds are large
% (lines here and in btc_soid_fit involving use_large)
%
%   See also:  BTC_SOID_DEMO, BTC_DEFINE, BTC_PAIRSNEEDED, BTC_SOID_PLOT, BTC_AUGCOORDS,
%   BTC_SOID_FIND, BTC_SOID_PLOT3D, BTC_PRED_DEMO.
%

if ~exist('res') res=10^(-5); end
%if ~exist('magfacts_list') magfacts_list=[0.25 0.5 0.75 1 1.25]; end
if ~exist('magfacts_list') magfacts_list=[0.25 0.5 1]; end
%
if ~exist('dict') dict=btc_define([]); end
nbtc=length(dict.codel); %10
if ~exist('csno'); csno=1; end
syms='.ox+*';
%
coords=dict.codel;
%
if ~exist('coordsets'); coordsets={'bg','bc','bd','de','ab'}; end
%
% set up structure of bounds for each pane
btc_bounds=[];
btc_bounds.bg=[];
btc_bounds.bg.coefs=[+2 -1; -2 -1];% b=2*g-1 and b=-2*g-1 are bounds
btc_bounds.bg.signs=[   +1,    +1];% b must be > 2*g-1 and > -2*g-1
btc_bounds.bc=[];
btc_bounds.bd=[];
btc_bounds.bd.coefs=[+1 -1; +1 +1; -1 +1; -1 -1]; %b=d-1, b=d+1,b=-d+1, b=-d-1 are bounds
btc_bounds.bd.signs=[   +1,    -1,    -1,    +1]; %b must be > d-1, < d+1, <-d+1, > -d-1
btc_bounds.de=[];
btc_bounds.ab=btc_bounds.bg;
%
btc_bounds.de=[];
%
for ics=1:length(coordsets)
    disp(sprintf(' coord set %2.0f->%10s',ics,coordsets{ics}));
end
csno=getinp('choice','d',[1 length(coordsets)],csno);
coords=coordsets{csno};
planes=btc_pairsneeded(coords,dict); %always symmetric setup (symmetric vs full doesn't matter here)
ncoords=length(coords); %always 2, but for compatibility with btc_soid_plot
nplanes=size(planes,1); %always 1, but for compatibility with btc_soid_plot
%
magfacts_list=getinp('magnification factors to be applied to spoke design','f',[0 3],magfacts_list);
res=getinp('resolution for determining maximum excursion','f',[10^(-8),10^(-1)],res);
%
% get the data from a named mat file
%
if ~exist('pdata_fn') pdata_fn='btc_allraysfixedb_mc'; end
pdata_fn=getinp('name of file with psychophysical data (e.g., btc_allraysfixedb_xx)','s',[],pdata_fn);
%
[edata,eb_avail,ds,condno]=btc_soid_getdata(pdata_fn,setfield([],'if_choose',1));
%
%set up the experiment directions
%transfer the experimental data
%and check that there are the right number of directions
%
edirs=[];
edirs_data=[];
ifok=1;
for iplane=1:nplanes
    pn=planes(iplane,:);
    edir=btc_edirs(pn);
    edirs=setfield(edirs,pn,edir);
    edirs_data=setfield(edirs_data,pn,edir);
    if isfield(edata,pn)
        if size(edata.(pn).thresh_mags,1)==edir.ndirs
            edirs_data.(pn).thresh_mags=edata.(pn).thresh_mags;
            edirs_data.(pn).thresh_vecs=edirs.(pn).uvecs.*repmat(edata.(pn).thresh_mags,1,2);
            if (eb_avail>0)
                if size(edata.(pn).thresh_mags_eblo,1)==edir.ndirs & size(edata.(pn).thresh_mags_eblo,1)==edir.ndirs
                    edirs_data.(pn).thresh_mags_eblo=edata.(pn).thresh_mags_eblo;
                    edirs_data.(pn).thresh_mags_ebhi=edata.(pn).thresh_mags_ebhi;
                    edirs_data.(pn).thresh_vecs_eblo=edirs.(pn).uvecs.*repmat(edata.(pn).thresh_mags_eblo,1,2);
                    edirs_data.(pn).thresh_vecs_ebhi=edirs.(pn).uvecs.*repmat(edata.(pn).thresh_mags_ebhi,1,2);
                else
                    disp(sprintf('Error bars for plane %2s from dataset in %s have wrong number of directions.',...
                        pn,pdata_fn))
                    ifok=0;
                end
            end
        else
            disp(sprintf('Threshold data for plane %2s from dataset in %s have wrong number of directions.',...
                pn,pdata_fn))
            ifok=0;
        end
    else
        disp(sprintf('Threshold data for plane %2s missing from dataset in %s.',...
            pn,pdata_fn))
        ifok=0;
    end
end
if (ifok==1)
    if (eb_avail>0)
        disp(sprintf('Threshold data and error bars successfully transferred for %2.0f planes.',nplanes))
    else
        disp(sprintf('Threshold data without error bars successfully transferred for %2.0f planes.',nplanes))
    end
else
    disp('Cannot proceed.');
end
%
allvecs=edirs_data.(planes).allvecs;
maxvecs=edirs_data.(planes).maxvecs;
ndirs=length(allvecs);
vlabel=cat(2,' planning from ',pdata_fn);
opts_plot=[];
opts_plot.datafield=strvcat('thresh_vecs','thresh_vecs_eblo','thresh_vecs_ebhi','maxvecs');
opts_plot.marker=strvcat('b','c','c','m');
opts_plot.cyclic=1;
opts_plot.tstring=vlabel;
%
btc_soid_plot(edirs_data,opts_plot);
set(gcf,'Position',[100 100 750 750]);
set(get(gca,'title'),'Interpreter','none')
allvecs=edirs_data.(planes).allvecs;
ndirs=length(allvecs);
for idir=1:ndirs;
    plot(allvecs{idir}(:,1),allvecs{idir}(:,2),'m.');
end
hold on;
bounds=btc_bounds.(planes);
if ~isempty(bounds)
    xvals=[-1:0.01:1];
    for ipoly=1:size(bounds.coefs,1)
        plot(xvals,polyval(bounds.coefs(ipoly,:),xvals),'r');
    end
end
%
%determine maximum offsets
offset_design=[];
hleg=[];
tleg=[];
for im=1:length(magfacts_list)
    magfact=magfacts_list(im);
    disp(sprintf(' computing maximum offset of origin for a magnification %5.2f of radial design in %s plane',magfact,planes));
    offset_design{im}.magfact=magfact;
    %along each direction in original design, determine how far one can
    %move until the boundary is exceeded on any spoke
    for idir=1:ndirs
        center=[0,0];
        offset_dir=maxvecs(idir,:);
        offset_dir_uvec=offset_dir/sqrt(sum(offset_dir.^2));
        dist=0;
        stepsize=1;
        while (stepsize>res)
            %test whether offset_dir*(dist+stepsize)+(any spoke*magfact) exceeds bounds
            all_ok=1;
            center_trial=offset_dir_uvec*(dist+stepsize);
            extremes=repmat(center_trial,ndirs,1)+maxvecs*magfact;
            %check whether all extremes are within [-1 1]
            if max(abs(extremes(:)))>1
                all_ok=0;
            end
            %
            %now check special bounds
            %btc_bounds.bd.coefs=[+1 -1; +1 +1; -1 +1; -1 -1]; %b=d-1, b=d+1,b=-d+1, b=-d-1 are bounds
            %btc_bounds.bd.signs=[   +1,    -1,    -1,    +1]; %b must be > d-1, < d+1, <-d+1, > -d-1     
            %
            %btc_bounds.bg.coefs=[+2 -1; -2 -1];% b=2*g-1 and b=-2*g-1 are bounds
            %btc_bounds.bg.signs=[   +1,    +1];% b must be > 2*g-1 and > -2*g-1
            if ~isempty(btc_bounds.(planes))
                for ibound=1:size(btc_bounds.(planes).coefs,1)
                    test=extremes(:,2)-(extremes(:,1)*btc_bounds.(planes).coefs(ibound,1)+btc_bounds.(planes).coefs(ibound,2));
                    if any(sign(test)~=btc_bounds.(planes).signs(ibound))
                        all_ok=0;
                    end
                end
            end
            if all_ok
                dist=dist+stepsize;
            end
            stepsize=stepsize/2;
        end
        offset_design{im}.dist(idir,:)=dist;
        offset_design{im}.center(idir,:)=offset_dir_uvec*dist;
        %
        disp(sprintf(' In direction %2.0f ([%7.3f %7.3f]), largest offset is [%7.3f %7.3f], dist of %6.4f',...
            idir,offset_dir,offset_dir_uvec*dist,dist));
        %for first mag factor, plot green arms originating at each allowed offset, which should go to the limits of the available gamut
        if (im==1)
            for tdir=1:ndirs
                plot(offset_dir_uvec(1,1)*dist+[0 magfact*maxvecs(tdir,1)],offset_dir_uvec(1,2)*dist+[0 magfact*maxvecs(tdir,2)],'g');
            end
        end
    end %idir
    ndirs_cyc=[[1:ndirs] 1];
    plot(offset_design{im}.center(ndirs_cyc,1),offset_design{im}.center(ndirs_cyc,2),'k');
    symb=syms(1+mod(im-1,length(magfacts_list)));
    hleg=[hleg;plot(offset_design{im}.center(ndirs_cyc,1),offset_design{im}.center(ndirs_cyc,2),cat(2,'k',symb))];
    tleg=strvcat(tleg,sprintf(' %7.4f',magfact));
end %im
hl=legend(hleg,tleg);
set(hl,'FontSize',7);
%
clear im magfact idir tdir dist stepsize offset_dir all_ok center_trial extremes ndirs_cyc symb hleg tleg
clear bounds ipoly xvals
clear idir ics iplane ifok 

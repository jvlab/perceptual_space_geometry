% btc_pairrays
% shows examples of textures generated by two specified parameters, along rays.
%   rays can be defined by boundary of accessible region (iexci=1)
%   or a Minkowski shape (iexci=2)
%
% requires pmethsum, dict, and dict_def from workspace saved after running btc_surv, but it can be run
% at low search resolution
%
%    See also: BTC_VERIFY, BTC_SURV, BTC_TEST, BTC_DEFINE, BTC_AUGCOORDS.
%
if ~exist('btc_surv_filename')
    btc_surv_filename='btc_surv_13Aug09.mat';
end
btc_surv_filename=getinp('file name of workspace saved from btc_surv','s',[],btc_surv_filename);
s=load(btc_surv_filename);
pmethsum=getfield(s,'pmethsum');
dict=getfield(s,'dict');
dict_def=getfield(s,'dict_def');
clear s;
%
if ~exist('mapsize')
    mapsize=32;
end
if ~exist('iexci') % 1 for use extremes, 2 for use specified values
    iexci=1;
end
if ~exist('mink') %Minkowski exponent if specified values used
    mink=2;
end
if ~exist('radfrac')
    radfrac=0.75;
end
if ~exist('nsteps')
    nsteps=5;
end
if ~exist('rayori')
    rayori=1;
end
if ~exist('ifpointplot')
    ifpointplot=1;
end
mapsize=getinp('sample map size','d',[4 1024],mapsize);
nsteps=getinp('nsteps','d',[1 20],nsteps);
rayori=getinp('1 for rays going down, 2 for rays going across','d',[1 2],rayori);
iexci=getinp('1 for rays to extremes, 2 for rays to specified values','d',[1 2],iexci);
ifpointplot=getinp('1 to make a plot of the parameter values used','d',[0 1],ifpointplot);
if (iexci==1)
    radfrac=getinp('fraction of distance to extreme paramater value','f',[0 1],radfrac);
    nrays=8;
end
if (iexci==2)
    if ~exist('nrays')
        nrays=8;
    end
    nrays=getinp('number of rays','d',[4 16],nrays);
    mink=getinp('Minkowski exponent (2=circle or ellipse, 1=diamond)','f',[0.01 100],mink);
end
%
opts_augcoords=[];
%
opts_makemaps=[];
opts_makemaps.show=0;
opts_makemaps.area=[mapsize mapsize];
opts_makemaps.nmaps=1;
auxopts=btc_auxopts;
auxopts.metro_opts.numiters=getinp('number of Metropolis iterations','d',[0 10000],auxopts.metro_opts.numiters);
%
npar=2;
coordnums_list=nchoosek([1:length(dict.codel)],npar);
allowed_pairs=[];
disp(' id expt   method(s).....');
for icomb=1:size(coordnums_list,1)
    exptname=btc_exptname(dict.codel(coordnums_list(icomb,:)),dict);
    outstring=sprintf('%2.0f->%5s',icomb,exptname);
    if (isfield(pmethsum,exptname))
        allowed_pairs=[allowed_pairs,icomb];
        methods=getfield(pmethsum,exptname);
        for im=1:length(methods)
            outstring=cat(2,outstring,'   ',sprintf('%s %s',methods{im}.name,methods{im}.variant_lab));
        end
        disp(outstring);
    end
end
which_pairs=getinp('pairs to examine','d',[1 size(coordnums_list,1)]);
if ~isempty(setdiff(which_pairs,allowed_pairs));
    disp('the following diallowed pairs are removed.');
    disp(setdiff(which_pairs, allowed_pairs));
    which_pairs=intersect(which_pairs,allowed_pairs);
end
pvsum=[];
for icomb=which_pairs
    exptname=btc_exptname(dict.codel(coordnums_list(icomb,:)),dict);
    methods=getfield(pmethsum,exptname);
    for iv=1:npar
        icn(iv)=dict.codel(coordnums_list(icomb,iv));
    end
    ic=char(icn);
    for im=1:length(methods)
        tstring=sprintf('%2.0f->%s %s %s',icomb,exptname,methods{im}.name,methods{im}.variant_lab);
        xymax=[];
        disp(cat(2,'processing ',tstring));
        %
        if (iexci==2)
            for iray=1:8
               xy=methods{im}.extremes(iray,:);
               disp(sprintf('direction %1.0f: %2s=[%5.2f %5.2f]',iray,ic,xy));
            end
            for iv=1:npar
               ax(iv)=getinp(sprintf('hemiaxis length (maximum abs value) for %s',ic(iv)),'f',[0 1]);
            end
        end
        %
        figure;
        set(gcf,'Position',[100 100 1200 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',tstring);
        for iray=1:nrays
           if (iexci==1)
               xymax(iray,:)=radfrac*methods{im}.extremes(iray,:);
           else
               theta=2*pi*(iray-1)/nrays;
               ct=cos(theta);
               st=sin(theta);
               nfac=(abs(ct)^mink+abs(st)^mink)^(1/mink);
               stfix=st/nfac;
               ctfix=ct/nfac;
               xymax(iray,:)=[ax(1)*ctfix ax(2)*stfix];
           end
           for istep=1:nsteps
                xy=(istep/nsteps)*xymax(iray,:);
                maps=repmat(0.5,mapsize,mapsize);
                spec=[];
                for iv=1:npar
                    spec=setfield(spec,ic(iv),xy(iv));
                end
                rv=btc_augcoords(spec,dict,opts_augcoords);
                if (rv.method{im}.ok_norm==1) & (rv.method{im}.ok_probs==1)
                    opts_makemaps.onaxis=0;
                    if (abs(xy(1))<=0.01) | (abs(xy(2))<=0.01)
                        opts_makemaps.onaxis=1;
                    end
                    [maps,ou,errs]=btc_makemaps(rv.method{im},opts_makemaps,dict,auxopts);
                    if ~isempty(errs)
                        disp(sprintf('btc_makemaps finds an error for ray %1.0f step %1.0f',iray,istep));
                        spec
                        disp(errs);
                        disp(ou);
                        maps=zeros(mapsize,mapsize);
                    end
                else
                    disp(sprintf('btc_augcoords cannot find parameters for ray %1.0f step %1.0f',iray,istep));
                    spec
                    disp(rv.method{im})
                end
                if (rayori==1)
                    subplot(nsteps,nrays,iray+(istep-1)*nrays,'align');
                end
                if (rayori==2)
                    subplot(nrays,nsteps,istep+(iray-1)*nsteps,'align');
                end
                imagesc(maps,[0 1]);
                axis equal;axis tight;colormap('gray');
                xlabel(sprintf('%2s=[%5.2f %5.2f]',ic,xy),'FontSize',7);
                set(gca,'XTick',[]);
                set(gca,'YTick',[]);
            end %istep
        end %iray
        if (ifpointplot)
            figure;
            set(gcf,'Position',[300 200 600 600]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',cat(2,tstring,' point plot'));
            for iray=1:nrays
                plot(xymax(iray,1)*[1:nsteps]/nsteps,xymax(iray,2)*[1:nsteps]/nsteps,'k.','MarkerSize',12); hold on;
                plot(xymax(iray,1)*[0 1],xymax(iray,2)*[0 1],'k','LineWidth',2); hold on;
            end
            plot([-1 1],[0 0],'k');
            plot([0 0],[-1 1],'k');
            for irad=1:4;
                plot(irad/4*cos(pi*[0:48]/24),irad/4*sin(pi*[0:48]/24),'k');
            end
            set(gca,'XTick',[-1:.5:1]);
            set(gca,'YTick',[-1:.5:1]);
            xlabel(ic(1));
            ylabel(ic(2));
            set(gca,'XLim',[-1 1]);
            set(gca,'YLim',[-1 1]);
            axis square;
        end
    end %im
end %icomb
clear maps ic icn icomb iray im iv tstring
clear methods ou outstring istep
clear rv spec which_pairs xy

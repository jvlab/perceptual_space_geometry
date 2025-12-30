%btc_mix2tex_demo.m: script to demonstrate donut Metropolis algorith as applied to mixtures for
% btc stimulus generation, and accumulate correlation and block count statistics
%
% This is derived from btc_metro_mix_demo, but uses btc_vec2tex_pickard to
% create the textures to be mixed, in addition to btc_augcoords.
%
% WARNING:  2-d autocorrelograms have H and V axes interchanged
%
% Notes
%  Computing the entropies and autocorrelations takes a noticeable amount of time; this can be eliminated
%     or reduced by limiting erange*d.
%  If entered parameters do not lead to non-negative probabilities, a new
%     set of parameters is requested.
%  2D Autocorrelogram scaling differs from that of btc_metro_mix_demo.
%
% Entropies are Treves-Panzeri debiased
% Error bars represent standard deviations, not standard errors
%
%   See also:  DONUT_METRO, GENMRFM, GETP2X2_ATG, BTC_NOPICK, ENTHIST,
%   BTC_METRO_DEMO, BTC_METRO_MIX_DEMO, BTC_AUGCOORDS, BTC_MIX2TEX_TUVW_SHOWSAMP.
%
disp('Warning:  H and V are interchanged on 2-d autocorrelation plot.')
%
switch getinp('0 for fast (debug), 1 for intermediate, 2 for slow','d',[0 2],1)
    case 0
        erange2d=[2:4];
        erange1d=[3:5];
        if ~exist('nruns') nruns=2; end
        if ~exist('numiters_def') numiters_def=100;end
        if ~exist('size_recur') size_recur=128; end %size of map to generate via recursion (only middle is used for Metropolis)
        if ~exist('size_metro') size_metro=96; end %size of map to run Metropolis algorithm on, and for stats
    case 1
        if ~exist('numiters_def') numiters_def=1000;end
    case 2
        if ~exist('numiters_def') numiters_def=10000;end
end
if ~exist('erange2d') erange2d=[2:4]; end
if ~exist('erange1d') erange1d=[3:10]; end
if ~exist('nruns') nruns=4; end
if ~exist('size_recur') size_recur=256; end %size of map to generate via recursion (only middle is used for Metropolis)
if ~exist('size_metro') size_metro=size_recur-24; end %size of map to run Metropolis algorithm on, and for stats
if ~exist('size_show') size_show=64; end %size of map sample to show by this routine

nruns=getinp('nruns','d',[1 100],nruns);
runs_show=getinp('list of runs to show','d',[1 nruns],1);
size_recur=getinp('size of map to generate via recursion','d',[64 1024],size_recur);
size_metro=getinp('size of map to run Metropolis donut algorithm on','d',[40 size_recur],min(size_metro,size_recur));
size_show=getinp('size of map to show','d',[16 size_metro],min(size_show,size_metro));
%
if ~exist('nmix')
    nmix=3;
end
nmix=getinp('number of textures to mix','d',[1 10],nmix);
%
% list of block sizes for entropies
%
if ~exist('entlist')
    elist1d=[erange1d;ones(1,length(erange1d))]';
    entlist=[[erange2d;erange2d]';elist1d;fliplr(elist1d)];
end
ncorrs=16;
disp('entlist')
disp(entlist)
%
select_metro=[1:size_metro]+round((size_recur-size_metro)/2);
select_show=[1:size_show]+round((size_metro-size_show)/2);
%
dict=btc_define;
tstringc=[];
%
subset_desc{1}='ga';
subset_desc{2}='bduwa';
subset_desc{3}='cduwa';
subset_desc{4}='betva';
subset_desc{5}='cetva';
subset_desc{6}='any pair';
use_aug=[0 0 0 0 0 1];
%
nsubset_types=length(subset_desc);
%
itypes=zeros(1,nmix);
tstring=cell(1,nmix);
tstringc=[];
tparams=cell(1,nmix);
tvec=zeros(nmix,length(dict.codel));
p2x2=cell(1,nmix);
opts_v2t_used=cell(1,nmix);
opts_v2t=struct;
opts_v2t.if_log=1;
opts_v2t.if_show=1;
map_components=zeros(size_recur,size_recur,nruns,nmix);
corrs=cell(1,nmix);
btc_types=cell(1,nmix);
nparams=zeros(1,nmix);
for imix=1:nmix
    if_ok=0;
    while (if_ok==0)
        for itype=1:nsubset_types
            disp(sprintf(' subset %1.0f: %s',itype,subset_desc{itype}));
        end
        itypes(imix)=getinp(sprintf('choice for component %1.0f',imix),'d',[1 nsubset_types]);
        if ~use_aug(itypes(imix))
            nparams(imix)=length(subset_desc{itypes(imix)});
            tparams{imix}=zeros(1,nparams(imix));
            %use btc_vec2tex_pickard: handles more than two params, but must be pickard
            tvec(imix,:)=0;
            tstring{imix}=[];
            for id=1:nparams(imix)
                let=subset_desc{itypes(imix)}(id);
                tparams{imix}(id)=getinp(sprintf('texture parameter value for %s',let),'f',[-1 1],0);
                tvec(imix,find(dict.codel==let))=tparams{imix}(id);
                tstring{imix}=cat(2,tstring{imix},sprintf('%s=%5.3f ',let,tparams{imix}(id)));
            end
            if (itypes(imix)==1)
                tvec(imix,find(dict.order==2))=tvec(imix,1).^2; %second-order params
                tvec(imix,find(dict.order==3))=tvec(imix,1).^3; %third-order params
            end
            %create all of the maps
            [map_components(:,:,:,imix),p2x2{imix},opts_v2t_used{imix}]=btc_vec2tex_pickard(tvec(imix,:),[size_recur size_recur nruns],opts_v2t);
            %
            if ~any(any(any(isnan(map_components(:,:,:,imix)))))
                 if_ok=1;
            else
                disp('***map generation fails in btc_vec2tex_pickard***');
                disp(opts_v2t_used{imix});
                disp(opts_v2t_used{imix}.msgs);
                disp('try again...')
            end
        else
            %use btc_augcoords: any two coordinates
            nparams(imix)=2;
            tparams{imix}=zeros(1,nparams(imix));
            spec=struct;
            btc_types{imix}=getinp(sprintf('texture parameter pair %2.0f (''tu'' or ''dt'', any two from %s)',imix,dict.codel),'s',[]);
            for id=1:nparams(imix)
                tparams{imix}(id)=getinp(sprintf('texture parameter value %2.0f for %s (dim %1.0f)',imix,btc_types{imix}(id),id),'f',[-1 1],0);
                spec.(btc_types{imix}(id))=tparams{imix}(id);
            end
            tstring{imix}=sprintf('%s=[%7.3f, %7.3f]',btc_types{imix},tparams{imix}(:));
            rv=btc_augcoords(spec);
            p2x2{imix}=rv.method{1}.p2x2;
            tvec(imix,:)=rv.method{1}.vec;
            errmsg=[];
            if all(p2x2{imix}(:)>=0)
                %generate the maps but skip Metropolis mixing, since this will happen when we combine the components
                btc_nometro_opts=btc_auxopts();
                btc_nometro_opts.metro_opts.numiters=0;
                irun=0;
                if_ok=1;
                while (irun<nruns & if_ok==1)
                    irun=irun+1;
                    [map_components(:,:,irun,imix),errmsg,metro_array]=btc_genstim(size_recur,size_recur,reshape(tparams{imix}(:),[1 1 2]),...
                        setfield([],'type',btc_types{imix}),btc_nometro_opts);
                    if_ok=and(if_ok,double(isempty(errmsg)));
                end
            else
                if_ok=0;
            end
            if (if_ok==1) %btc_genstim does not plot
                figure;
                set(gcf,'Position',[100 100 1200 800]);
                imagesc(map_components(:,:,1,imix),[0 1]);
                colormap gray;
                axis equal;
                axis tight;
                title(tstring{imix});
            else
                disp('***map generation fails in btc_augcoords***');
                disp(rv);
                disp(errmsg);
                disp('try again...')
            end
        end
    end %if_ok
    set(gcf,'Name',sprintf('component %1.0f',imix));
    set(gcf,'NumberTitle','off');
    %
    corrs{imix}=getcorrs_p2x2(p2x2{imix});
    disp(sprintf('%4.0f map(s) made for mixture component %1.0f (%s)',nruns,imix,tstring{imix}));
    disp(sprintf('entropy per unit area (if Pickard): %7.4f',corrs{imix}.entropy));
    disp(sprintf('minimum probability of any 2x2 block: %7.4f',min(p2x2{imix}(:))));
    tstringc=cat(2,tstringc,'[',tstring{imix},']+');
end %imix
tstringc=tstringc(1:end-1);
%
%use same defaults as btc modules
if exist('auxopts')
    auxopts=btc_auxopts(auxopts);
else
    auxopts=btc_auxopts;
end
%
auxopts.metro_opts.numiters=getinp('number of Metropolis iterations','d',[0 10^6],numiters_def);
auxopts.metro_opts.sampfreq_map=round(auxopts.metro_opts.numiters/10); %frequency to save sampled maps during Metropolis, also used to calculate statistics
auxopts.metro_opts.sampfreq_stat=0; %never do stats during Metropolis; here we calculate statistics are based on sampled maps
auxopts.metro_opts.showfreq_map=0;
auxopts.metro_opts.sampfreq_map=getinp('frequency to sample stats and map for final display','d',[0 10^6],auxopts.metro_opts.sampfreq_map);
auxopts.metro_opts.nf_frac=(1/2^8)*getinp('target flip fraction*256','f',[0 2^8],auxopts.metro_opts.nf_frac*256);
auxopts.metro_opts.nf_dist=getinp('nf_dist (0->const, 1->unif, 2->Bern, 3->Pois, 4->exp)','d',[0 4],auxopts.metro_opts.nf_dist);
%
%setup for donut algorithm
%
donut=[];
donut.name='donut';
donut.matrix=[1 1 1;1 0 1; 1 1 1];
donut=glider_addcoords(donut);
metro_opts=auxopts.metro_opts;
%
nmap_samps=1+floor(auxopts.metro_opts.numiters/auxopts.metro_opts.sampfreq_map);
entropies=zeros(size(entlist,1),nmap_samps,nmix,nruns);
xcs=zeros(1+2*ncorrs,1+2*ncorrs,nmap_samps,nmix,nruns);
autocors=zeros(1+ncorrs,2,nmap_samps,nmix,nruns); %second dimension is H or V, third is sample 
%
vecs=zeros(length(dict.codel),nmix,nmap_samps,nruns); %empirical coordinates
for irun=1:nruns
    btc_nometro_opts=btc_auxopts();
    btc_nometro_opts.metro_opts.numiters=0;
    %
    map_recur=squeeze(map_components(:,:,irun,:));
    map_start=map_recur(select_metro,select_metro,:);
    [map_donut,metro_samp,metro_optsused]=donut_metro(2,donut,map_start,setfield(metro_opts,'map_show_name',tstringc));
    for ishow=1:nmap_samps
        if (ishow==1)
            mapstack=map_start;
            tlabel='start';
        else
            mapconcat=metro_samp.map.map_synth(:,:,ishow-1);
            for imix=1:nmix
                mapstack(:,:,imix)=mapconcat(:,[1:size(mapconcat,1)]+(imix-1)*size(mapconcat,1)); %extract from concatenated map
            end
            tlabel=sprintf('%5.0f gens',metro_samp.map.iter(ishow-1));
        end
        %do statistical calculations
        %image statistics
        for imix=1:nmix
            counts=btc_map2counts(mapstack(:,:,imix));
            vecs(:,imix,ishow,irun)=btc_corrs2vec(getcorrs_p2x2(counts/sum(counts(:))))';
        end
        %entropies
        for ient=1:size(entlist,1)
            for imix=1:nmix
                blockcounts=mapubi(mapstack(:,:,imix),entlist(ient,:));%creates a list of unique blocks,
                bps=blockcounts(:,end);
                entropies(ient,ishow,imix,irun)=enthist(bps'); %Treves-Panzeri debiased value
            end
        end
        for imix=1:nmix
            map=mapstack(:,:,imix);
            xc=xcorr2(map-mean(map(:)));
            xc_needed=xc(length(map)+[-ncorrs:ncorrs],length(map)+[-ncorrs:ncorrs]);
            xc_needed=xc_needed./((length(map)-abs([-ncorrs:ncorrs]))'*(length(map)-abs(-ncorrs:ncorrs)));
            xcs(:,:,ishow,imix,irun)=xc_needed/xc_needed(ncorrs+1,ncorrs+1);
            autocors(:,1,ishow,imix,irun)=xcs(ncorrs+1+[0:ncorrs],ncorrs+1,ishow,imix,irun);
            autocors(:,2,ishow,imix,irun)=xcs(ncorrs+1,ncorrs+1+[0:ncorrs],ishow,imix,irun)';
        end
    end %ishow
    if ismember(irun,runs_show)
        %
        %maps
        %
        for imix=1:nmix
            figure;
            set(gcf,'Position',[50 200 1000 700]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',sprintf('map (run %1.0f) %s',irun,tstring{imix}));
            [nr,nc]=nicesubp(nmap_samps,0.7);
            for ishow=1:nmap_samps
                if (ishow==1)
                    map=map_start(:,:,imix);
                    tlabel='start';
                else
                    mapconcat=metro_samp.map.map_synth(:,:,ishow-1);
                    map=mapconcat(:,[1:size(mapconcat,1)]+(imix-1)*size(mapconcat,1)); %extract from concatenated map
                    tlabel=sprintf('%5.0f gens',metro_samp.map.iter(ishow-1));
                end
                subplot(nr,nc,ishow);
                imagesc(map(select_show,select_show),[0 1]);
                hold on;
                colormap gray;axis equal;axis tight;
                set(gca,'XTick',[]);
                set(gca,'YTick',[]);
                xlabel(tstring{imix});
                ylabel(tlabel);
            end
        end %imix
        %
        %entropies
        %
        for imix=1:nmix
            figure;
            set(gcf,'Position',[50 200 1000 700]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',sprintf('entropy per pixel (run %1.0f) %s',irun,tstring{imix}));
            [nr,nc]=nicesubp(size(entlist,1),0.7);
            for ient=1:size(entlist,1)
                subplot(nr,nc,ient);
                hold on;
                plot([0 metro_samp.map.iter],squeeze(entropies(ient,:,imix,irun))/prod(entlist(ient,:)),'k.-');
                set(gca,'YLim',[0 1.1]);
                set(gca,'YTick',[0 0.5 1]);
                set(gca,'XTick',[0 metro_samp.map.iter]);
                set(gca,'XLim',[0 metro_samp.map.iter(end)]);
                xlabel(tstring{imix});
                ylabel(sprintf('%2.0f x %2.0f',entlist(ient,:)));
            end   
            %
            %autocors
            %
            figure;
            set(gcf,'Position',[50 200 1000 700]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',sprintf('autocors (run %1.0f) %s',irun,tstring{imix}));
            [nr,nc]=nicesubp(nmap_samps,0.7);
            for ishow=1:nmap_samps
                subplot(nr,nc,ishow);
                hold on;
                plot([0:ncorrs],autocors(:,1,ishow,imix,irun),'k.-');
                plot([0:ncorrs],autocors(:,2,ishow,imix,irun),'ko:');
                legend('dim 1','dim 2');
                set(gca,'YLim',[-0.25 1]);
                set(gca,'YTick',[0 0.5 1]);
                set(gca,'XLim',[-0.5 ncorrs+0.5]);
                xlabel(tstring{imix});
                if (ishow==1) ylabel('gen 0'); end
                if (ishow>1) ylabel(sprintf('gen %1.0f',metro_samp.map.iter(ishow-1))); end
            end
        end %imix
    end  %ismember(iruns,runs_show) 
end %end of run
%
%show averages across runs
%
for imix=1:nmix
    %entropies
    figure;
    set(gcf,'Position',[50 200 1000 700]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',cat(2,'entropy per pixel (avgs) ',tstring{imix}));
    [nr,nc]=nicesubp(size(entlist,1),0.7);
    for ient=1:size(entlist,1)
        subplot(nr,nc,ient);
        hold on;
        ents=squeeze(mean(entropies(ient,:,imix,:),4))/prod(entlist(ient,:));
        if (nruns>1)
            errorbar([0 metro_samp.map.iter],ents,squeeze(std(entropies(ient,:,imix,:),[],4)/prod(entlist(ient,:))),'k.-');
        else
            plot([0 metro_samp.map.iter],ents,'k.-');
        end
        set(gca,'YLim',[0 1.1]);
        set(gca,'YTick',[0 0.5 1]);
        set(gca,'XTick',[0 metro_samp.map.iter]);
        set(gca,'XLim',[0 metro_samp.map.iter(end)]);
        xlabel(tstring{imix});
        ylabel(sprintf('%2.0f x %2.0f',entlist(ient,:)));
    end   
    %autocors
    figure;
    set(gcf,'Position',[50 200 1000 700]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',cat(2,'autocors  (avgs) ',tstring{imix}));
    [nr,nc]=nicesubp(nmap_samps,0.7);
    for ishow=1:nmap_samps
        subplot(nr,nc,ishow);
        hold on;
        if (nruns>1)
            errorbar([0:ncorrs],mean(autocors(:,1,ishow,imix,:),5),std(autocors(:,1,ishow,imix,:),[],5),'k.-'); 
            errorbar([0:ncorrs],mean(autocors(:,2,ishow,imix,:),5),std(autocors(:,2,ishow,imix,:),[],5),'ko:'); 
        else
            plot([0:ncorrs],mean(autocors(:,1,ishow,imix,:),5),'k.-');
            plot([0:ncorrs],mean(autocors(:,2,ishow,imix,:),5),'ko:');
        end
        legend('dim 1','dim 2');
        set(gca,'YLim',[-0.25 1]);
        set(gca,'YTick',[0 0.5 1]);
        set(gca,'XLim',[-0.5 ncorrs+0.5]);
        xlabel(tstring{imix});
        if (ishow==1) ylabel('gen 0'); end
        if (ishow>1) ylabel(sprintf('gen %1.0f',metro_samp.map.iter(ishow-1))); end
    end   
    %2-d autocors
    figure;
    set(gcf,'Position',[50 200 1200 750]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',cat(2,'2d autocors  (avgs) ',tstring{imix}));
    [nr,nc]=nicesubp(nmap_samps,0.7);
    u=mean(xcs(:,:,:,imix,:),5);
    u(ncorrs+1,ncorrs+1,:)=0;
    umax=max(abs(u(:)));
    for ishow=1:nmap_samps
        subplot(nr,nc,ishow);
 %       imagesc([-ncorrs:ncorrs],[-ncorrs:ncorrs],(mean(xcs(:,:,ishow,imix,:),5))',max(abs(tvec(:)))*[-1 1]); hold on;
        imagesc([-ncorrs:ncorrs],[-ncorrs:ncorrs],(mean(xcs(:,:,ishow,imix,:),5))',umax*1.1*[-1 1]); hold on;
        colormap jet;
        colorbar;
        axis square;
        xlabel('dim 1');
        ylabel('dim 2');
        xlabel(tstring{imix});
        if (ishow==1) ylabel('gen 0'); end
        if (ishow>1) ylabel(sprintf('gen %1.0f',metro_samp.map.iter(ishow-1))); end
    end
end %imix
%
%compare empirical coordinates with desired (tvec)
%
tvec_goal=mean(tvec,1);
vecs_avg=mean(vecs,4); %mean across runs; d1=coord vec, d2=mix component, d3=sample number
figure;
set(gcf,'Position',[50 200 1000 700]);
set(gcf,'NumberTitle','off');
set(gcf,'Name','convergence of each coordinate');
[nr,nc]=nicesubp(length(dict.codel),0.7);
for ilet=1:length(dict.codel)
    subplot(nr,nc,ilet);
    plot([0:nmap_samps-1],squeeze(vecs_avg(ilet,:,:))');
    hold on;
    plot([0 nmap_samps-1],repmat(tvec_goal(ilet),[1 2]),'k--');
    set(gca,'XTick',[0:nmap_samps-1]);
    set(gca,'XTickLabel',[0 metro_samp.map.iter])
    set(gca,'XLim',[0 nmap_samps-1]);
    xlabel('iterations');
    ylabel(dict.codel(ilet));
    legend(strvcat(num2str([1:nmix]'),'target'),'FontSize',7,'Location','best');
    title(sprintf('target: %6.3f',mean(tvec(:,ilet))));
end %ilet
%
%project coords at each sample on line between sample and target
%separately for each run
%
figure;
set(gcf,'Position',[150 200 1000 600]);
set(gcf,'NumberTitle','off');
set(gcf,'Name','convergence by run');
for imix=1:nmix
    tvec_diff=(tvec(imix,:)-tvec_goal)';
    if max(abs(tvec_diff))>0.01       
        subplot(1,nmix,imix);
        for irun=1:nruns
            diffs=reshape(vecs(:,imix,:,irun),[length(dict.codel),nmap_samps])-repmat(tvec_goal(:),[1 nmap_samps]);
            projs=diffs'*tvec_diff./(tvec_diff'*tvec_diff);
            plot([0:nmap_samps-1],projs);
        hold on;
        end
        plot([0 nmap_samps-1],zeros(1,2),'k--');
        set(gca,'XTick',[0:nmap_samps-1]);
        set(gca,'XTickLabel',[0 metro_samp.map.iter])
        set(gca,'XLim',[0 nmap_samps-1]);
        set(gca,'YLim',[-.1 1]);
        xlabel('iterations');
        ylabel('fractional distance');
        title(sprintf('component %1.0f',imix));
    end
end %ilet
%
%global summary
%
disp('target statistics for each component prior to mixing')
disp(tvec);
disp('actual statistics for each component prior to mixing')
disp(mean(vecs(:,:,1,:),4)');
disp('difference (component actual-component target)')
disp(mean(vecs(:,:,1,:),4)'-tvec);
disp(' ');
disp('target statistics for final mixture')
disp(mean(tvec,1));
disp('actual statistics for final mixture')
disp(mean(mean(vecs(:,:,end,:),4),2)');
disp('difference (actual-target)')
disp(mean(mean(vecs(:,:,end,:),4),2)'-mean(tvec,1));
disp(' ');
disp('actual statistics for each component after mixing')
disp(mean(vecs(:,:,end,:),4)');
disp('difference (component actual-mixture target)')
disp(mean(vecs(:,:,end,:),4)'-repmat(mean(tvec,1),nmix,1));
disp(' ')
disp('maximum difference between components prior to mixing')
disp((max(mean(vecs(:,:,1,:),4),[],2)-min(mean(vecs(:,:,1,:),4),[],2))');
disp('maximum difference between components after mixing')
disp((max(mean(vecs(:,:,end,:),4),[],2)-min(mean(vecs(:,:,end,:),4),[],2))');



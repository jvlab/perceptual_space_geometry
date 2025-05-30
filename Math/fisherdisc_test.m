%fisherdisc_test: test routine for Fisher discriminant
%
%  See also:
%  FISHERDISC, GNORMCOR, FISHERDISC_DEMO.
%
if ~exist('scenarios') scenarios=cell(0); end
if ~exist('opts')
    opts=[];
    opts.nshuffle=-1; %could try Inf, 100
    opts.xvalid=1; %could try 0
    opts.nshuffle_max=200; %could try 10000;
    opts.sub1=1;
    opts.classifiers={'halfway','mapequal','mapbayes'};
end
%
if getinp('1 to reset the random number generator','d',[0 1],1);
    rand('state',0);
    randn('state',0);
end
%
% the covariances are chosen once and for as a random symmetric matrix with all-positive eigenvalues
% via z=rand(nfeat,nfeat)-0.5; covars=z*z';
%
scenarios{1}.name='equal spherical covariances';
scenarios{1}.nfeat=5;
scenarios{1}.signal=[-.2 .3 .4 -.1 .5]; 
scenarios{1}.sigsize=[0 0.1 0.3 1 3 10];
scenarios{1}.nsamps=[10 20];
scenarios{1}.covars=repmat(eye(scenarios{1}.nfeat),[1 1 2]);
%
scenarios{2}=scenarios{1};
scenarios{2}.name='equal ellipsoidal covariances';
scenarios{2}.covars=3*repmat([...
    0.4443   -0.2201    0.0548   -0.1587    0.2033;...
   -0.2201    0.5462   -0.1381    0.0892   -0.0824;...
    0.0548   -0.1381    0.0523   -0.0628    0.0572;...
   -0.1587    0.0892   -0.0628    0.3795    0.0234;...
    0.2033   -0.0824    0.0572    0.0234    0.3279],[1 1 2]);
%
scenarios{3}=scenarios{2};
scenarios{3}.name='unequal ellipsoidal covariances';
scenarios{3}.covars(:,:,1)=0.5*scenarios{2}.covars(:,:,1);
scenarios{3}.covars(:,:,2)=5*[...
    0.2104   -0.1499    0.1738    0.0467    0.1005;...
   -0.1499    0.2692    0.0175    0.0346    0.0714;...
    0.1738    0.0175    0.4308    0.2221    0.1538;...
    0.0467    0.0346    0.2221    0.1777    0.0504;...
    0.1005    0.0714    0.1538    0.0504    0.2497];
%
for is=1:length(scenarios)
    disp(sprintf('%3.0f->%32s, nfeatures=%2.0f nsamps=[%4.0f %4.0f]',...
        is,scenarios{is}.name,scenarios{is}.nfeat,scenarios{is}.nsamps));
end
islist=getinp('choice','d',[1 length(scenarios)],[1:length(scenarios)]);
if getinp('1 to override number of samples','d',[0 1],0);
    for isp=1:length(islist)
        nover(isp,:)=getinp(sprintf('number of samples for the two groups in scenario %2.0f',islist(isp)),...
            'd',[1 Inf],scenarios{islist(isp)}.nsamps);
    end
else
    nover=[];
end
scattervals=cell(0);
results=cell(0,0);
ou=cell(0,0);
for isp=1:length(islist)
    is=islist(isp);
    sc=scenarios{is};
    if ~isempty(nover)
        sc.nsamps=nover(isp,:);
    end
    for ig=1:2
        scattervals{isp,ig}=gnormcor(sc.covars(:,:,ig),sc.nsamps(ig));
    end
    disp(sprintf('made scatter values for scenario %3.0f->%32s, nfeatures=%2.0f nsamps=[%4.0f %4.0f]',...
        is,sc.name,sc.nfeat,sc.nsamps));
    for isigsize=1:length(sc.sigsize)
        sigsize=sc.sigsize(isigsize);
        disp(sprintf('  running sigsize=%7.3f',sigsize));
        %create samples by subtracting signal from group 1, adding to group 2
        samps=[scattervals{isp,1}-repmat(sigsize*sc.signal',1,sc.nsamps(1)),...
            scattervals{isp,2}+repmat(sigsize*sc.signal',1,sc.nsamps(2))];
        tags=[repmat(1,[1 sc.nsamps(1)]),repmat(2,[1 sc.nsamps(2)])];
        [results{isp,isigsize},ou{isp,isigsize}]=fisherdisc(samps,tags,opts);
    end %isigsize
    %
    % plot discriminants
    %
    figure;
    set(gcf,'Position',[100 150 1200 700]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',cat(2,'classifiers: ',sc.name));
    [nr,nc]=nicesubp(length(sc.sigsize),0.7);
    for isigsize=1:length(sc.sigsize)
        rj=results{isp,isigsize};
        subplot(nr,nc,isigsize);
        hp(1)=plot([1:sc.nfeat],sc.signal/sqrt(sum(sc.signal.^2)),'k-','LineWidth',1);hold on;
        hp(2)=plot([1:sc.nfeat],rj.discriminant,'b-','LineWidth',2);hold on;
        labels={'signal','classifier'};
        if opts.xvalid==1
            hp(3)=plot([1:sc.nfeat],rj.discriminant_jdebiased,'b-','LineWidth',1);hold on;
            labels{3}='debiased';
            for ifeat=1:sc.nfeat
                plot([ifeat ifeat],rj.discriminant_jdebiased(ifeat)+[-1 1]*rj.discriminant_jsem(ifeat),'b-','LineWidth',1);hold on;
            end
        end %opts.xvalid
        plot([0.5 sc.nfeat+0.5],[0 0],'k-');hold on;
        xlabel('feature');
        set(gca,'XTick',[1:sc.nfeat]);
        set(gca,'XTickLabel',strvcat('1',repmat(' ',sc.nfeat-2,1),sprintf('%3.0f',sc.nfeat)));
        set(gca,'XLim',[0.5 sc.nfeat+0.5]);
        ylabel('weight');
        set(gca,'YTick',[-1:0.5:1]);
        set(gca,'YLim',[-1 1]);
        legend(hp,labels,'Location','South','FontSize',7);
        title(cat(2,sc.name,sprintf(' sigsize=%4.1f',sc.sigsize(isigsize))),'FontSize',7);
    end
    %
    % if opts.classifiers>0 then plot fraction correct as function of
    % signal size, for direct calc and (if cross-validated) cross-validated calc
    %
    if ~isempty(opts.classifiers)
        figure;
        set(gcf,'Position',[100 200 800 600]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,'performance: ',sc.name));
        frac_corr=zeros(length(opts.classifiers),length(sc.sigsize),1+opts.xvalid);
        lstring=[];
        for ixv=0:opts.xvalid
            for ic=1:length(opts.classifiers)
                if (ixv==0)
                    tstring='raw';
                    cname=cat(2,'fc_',opts.classifiers{ic},'_cmatrix');
                else
                    tstring='xvalid';
                    cname=cat(2,'xv_',opts.classifiers{ic},'_cmatrix');
                end
                for isigsize=1:length(sc.sigsize)
                    cmat=getfield(results{isp,isigsize},cname);
                    frac_corr(ic,isigsize,1+ixv)=sum(diag(cmat))/sum(sc.nsamps);
                end
                lstring=strvcat(lstring,cat(2,opts.classifiers{ic},' ',tstring));
            end %ic
        end %ixv
        plot([1:length(sc.sigsize)],frac_corr(:,:,1),'-','LineWidth',2); hold on;
        if (opts.xvalid==1)
            plot([1:length(sc.sigsize)],frac_corr(:,:,2),':','LineWidth',2); hold on;
        end
        set(gca,'XTick',[1:length(sc.sigsize)]);
        set(gca,'XTickLabel',sc.sigsize);
        xlabel('sigsize');
        set(gca,'YLim',[0 1]);
        ylabel('frac corr');
        legend(lstring,'Location','SouthEast');
        title(sc.name);
    end %opts.classifiers
    %
    % plot shuffle analysis
    %
    nhist=40;
    hist_centers=([1:nhist]-0.5)/nhist;
    if ~(opts.nshuffle==0)
        figure;
        set(gcf,'Position',[100 150 1200 700]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,'shuffle variance analysis: ',sc.name));
        for isigsize=1:length(sc.sigsize)
            rj=results{isp,isigsize};
            rjs=rj.shuffle;
            subplot(nr,nc,isigsize);
            hdata=hist(rjs.ss_class./rjs.ss_total,hist_centers);hold on;
            hh=bar(hist_centers,hdata);
            set(hh,'FaceColor',[1 1 1]);
            set(gca,'XLim',[0 1]);
            xlabel('class/total');
            hd(1)=plot(rj.ss_class/rj.ss_total*[1 1],get(gca,'YLim'),'k-','LineWidth',2);hold on;
            p=sum(rjs.ss_class./rjs.ss_total<rj.ss_class/rj.ss_total)/length(rjs.ss_class);
            title(cat(2,sc.name,sprintf(' sigsize=%4.1f',sc.sigsize(isigsize))),'FontSize',7);
            hl=legend(hd,sprintf('p=%6.4f',p),'Location','NorthWest');
            set(hl,'FontSize',7);
        end
        if ~isempty(opts.classifiers)
            figure;
            set(gcf,'Position',[100 200  800 600]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',cat(2,'shuffle frac corr analysis: ',sc.name));
            frac_corr=zeros(length(opts.classifiers),length(sc.sigsize),2);
            lstring=[];
            for ish=1:2
                for ic=1:length(opts.classifiers)
                    if (ish==1) lstring=strvcat(lstring,cat(2,opts.classifiers{ic},' raw')); end
                    if (ish==2) lstring=strvcat(lstring,cat(2,opts.classifiers{ic},' shuffled')); end
                    for isigsize=1:length(sc.sigsize)
                        if (ish==1)
                            rj=results{isp,isigsize};
                            cmat=getfield(rj,cat(2,'fc_',opts.classifiers{ic},'_cmatrix'));
                            frac_corr(ic,isigsize,ish)=sum(diag(cmat))/sum(sc.nsamps);
                        end
                        if (ish==2)
                            rjs=rj.shuffle;
                            cmat_shuf=getfield(rjs,cat(2,'fc_',opts.classifiers{ic},'_cmatrix'));
                            frac_corr(ic,isigsize,ish)=sum(diag(sum(cmat_shuf,3)))/sum(sc.nsamps)/size(cmat_shuf,3);
                        end
                    end %isigsize
                end %ic
            end %ish
            plot([1:length(sc.sigsize)],frac_corr(:,:,1),'-','LineWidth',2); hold on;
            plot([1:length(sc.sigsize)],frac_corr(:,:,2),':','LineWidth',2); hold on;
            set(gca,'XTick',[1:length(sc.sigsize)]);
            set(gca,'XTickLabel',sc.sigsize);
            xlabel('sigsize');
            set(gca,'YLim',[0 1]);
            ylabel('frac corr');
            legend(lstring,'Location','SouthEast');
            title(sc.name);
        end %shuffle and classifiers
        clear r rj p hh hdata hd cmat cmat_shuf hl
    end %shuffle
end %is (scenario)
clear rj nr nc ic ig is ifeat ixv frac_corr lstring labels hp

%fisherdisc_demo: demo routine for Fisher discriminant
%
%  This uses spectral data from Goldfine's navigate/stop navigate
%  paradigm, it is derived from fisherdisc_test that uses artificial data.
%  See also:
%  FISHERDISC_TEST, FISHERDISC, GNORMCOR, SETFIELDS.
%
if ~exist('scenarios') scenarios=cell(0); end
% edit this to change the default frequency partitions
% each row is lower and upper frequency
if ~exist('fpart')
    fpart=[3.5 5.5;5.5 7.5;7.5 9.5;9.5 12.5;12.5 15.5;15.5 19.5;19.5 23.5;23.5 28.5;28.5 35.5;35.5 44.5;44.5 54.5];
end
disp('frequency partitions in Hz (this is the variable "fpart"):')
disp(fpart);
disp('if you want to change this, exit now and re-define fpart.');
%
% edit this to change the default options
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
% read the data
class_name={'Nav','StopNav'}; %each class (group) is a condition
scenario_name={'CP2','P3'}; %each scenario is a channel
subj_name='normal2';
datapath=getinp('data path','s',[],'../mcs/');
for is=1:2
    for ig=1:2
        file_name=cat(2,datapath,subj_name,class_name{ig},'Spectrum',scenario_name{is},'.mat');
        field_name=cat(2,'spectra',scenario_name{is},lower(class_name{ig}));
        disp(sprintf('attempting to load %s from %s',field_name,file_name));
        d{is,ig}=load(file_name);
        spectra{is,ig}=getfield(d{is,ig},field_name);
    end
end
freq_datasampling=d{1,1}.frequencyRecorded;
freqs=d{1,1}.f{1}*d{1,1}.frequencyRecorded;
freq_spectralsampling=freqs(2)-freqs(1);
mt_params=d{1,1}.params;
%
segsize=getinp('segment size (number of trials to consider contiguous','d',[0 Inf],2);
delrad=getinp('delrad (number of adjacent trials to consider dependent','d',[0 Inf],Inf);
%
% set up the scenarios
%
for is=1:2
    scenarios{is}.name=cat(2,scenario_name{is},': ',subj_name,' ',class_name{1},' vs ',class_name{2});
    scenarios{is}.sigsize=1; %just for compatibility with fisherdisc_test
    ifeat=0; %each feature is mean power in a non-empty partition of frequencies
    scenarios{is}.data=[];
    for ipart=1:size(fpart,1)
        fbins=intersect(find(freqs>=fpart(ipart,1)),find(freqs<fpart(ipart,2)));
        if length(fbins)>0
            ifeat=ifeat+1;
            for ig=1:2
                scenarios{is}.data{ig}(ifeat,:)=mean(spectra{is,ig}(fbins,:),1);
            end
        end %length(fbins)
    end %ipart
    for ig=1:2
        scenarios{is}.nsamps(1,ig)=size(scenarios{is}.data{ig},2);
    end
    scenarios{is}.nfeat=ipart;
    %difference of means, not necessariliy the discriminant
    dmean=mean(scenarios{is}.data{2},2)-mean(scenarios{is}.data{1},2);
    dmean=dmean./sqrt(sum(dmean.^2));
    scenarios{is}.signal=dmean;
    nsamps=sum(scenarios{is}.nsamps); %total number of samples
    %set up the segments
    scenarios{is}.segs=cell(0);
    for iseg=1:ceil(nsamps/segsize);
        scenarios{is}.segs{iseg}=[1+segsize*(iseg-1):min(segsize*iseg,nsamps)];
    end
end
%
%
for is=1:length(scenarios)
    disp(sprintf('%3.0f->%32s, nfeatures=%2.0f nsamps=[%4.0f %4.0f]',...
        is,scenarios{is}.name,scenarios{is}.nfeat,scenarios{is}.nsamps));
end
islist=getinp('choice','d',[1 length(scenarios)],[1:length(scenarios)]);
results=cell(0);
ou=cell(0);
%
% run and plot each requested scenario
%
for isp=1:length(islist)
    is=islist(isp);
    sc=scenarios{is};
    tags=[repmat(1,[1 sc.nsamps(1)]),repmat(2,[1 sc.nsamps(2)])];
    %
    % do the Fisher analysis
    %
    [results{isp},ou{isp}]=fisherdisc([sc.data{1},sc.data{2}],tags,...
        setfields(opts,{'segs','delrad'},{sc.segs,delrad}));
    %
    % plot discriminants
    %
    figure;
    set(gcf,'Position',[100 150 1300 900]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',sc.name);
    subplot(2,1,1);
    for isigsize=1:1 %for compatibility with fisherdisc_test (no "sigsize" here)
        rj=results{isp};
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
        %fancy argument to use bracketed ranges for tick labels
        set(gca,'XTickLabel',cellstr(reshape(sprintf('[%4.1f %4.1f]',fpart'),[11 sc.nfeat])'));
        set(gca,'XLim',[0.5 sc.nfeat+0.5]);
        set(gca,'FontSize',7);
        ylabel('weight','FontSize',10);
        set(gca,'YTick',[-1:0.5:1]);
        set(gca,'YLim',[-1 1]);
        legend(hp,labels,'Location','South','FontSize',7);
        title(cat(2,sc.name,sprintf(': segsize=%3.0f delrad=%3.0f',segsize,delrad)),'FontSize',10);
    end
    %
    % if opts.classifiers>0 then plot fraction correct 
    % for direct calc and (if cross-validated) cross-validated calc
    %
    if ~isempty(opts.classifiers)
        subplot(2,3,4);
        frac_corr=zeros(length(opts.classifiers),length(sc.sigsize),1+opts.xvalid);
        lstring=[];
        for ixv=0:2*opts.xvalid %2 if xvalid =1, to allow for plotting of xv_ and sg_
            for ic=1:length(opts.classifiers)
                if (ixv==0)
                    tstring='raw';
                    cname=cat(2,'fc_',opts.classifiers{ic},'_cmatrix');
                end
                if (ixv==1)
                    tstring='xvalid';
                    cname=cat(2,'xv_',opts.classifiers{ic},'_cmatrix');
                end
                if (ixv==2)
                    tstring='seg';
                    cname=cat(2,'sg_',opts.classifiers{ic},'_cmatrix');
                end
                for isigsize=1:1 %for compatibility
                    cmat=getfield(results{isp},cname);
                    frac_corr(ic,1,1+ixv)=sum(diag(cmat))/sum(sc.nsamps);
                end
                lstring=strvcat(lstring,cat(2,opts.classifiers{ic},' ',tstring));
            end %ic
        end %ixv
        plot([0 1],frac_corr(:,[1 1],1),'-','LineWidth',2); hold on;
        if (opts.xvalid==1)
            plot([0 1],frac_corr(:,[1 1],2),'-.','LineWidth',2); hold on;
            plot([0 1],frac_corr(:,[1 1],3),':','LineWidth',2); hold on;
        end
        plot([0 1],[0.5 0.5],'k-');
        set(gca,'XTick',[0 1]);
        set(gca,'XTickLabel',' ');
        set(gca,'XLim',[0 1]);
        set(gca,'YLim',[0 1]);
        ylabel('frac corr');
        hl=legend(lstring,'Location','SouthEast');
        set(hl,'FontSize',7);
        title('raw vs. x-validated','FontSize',7);
    end %opts.classifiers
    %
    % plot shuffle analysis
    %
    nhist=40;
    hist_centers=([1:nhist]-0.5)/nhist;
    if ~(opts.nshuffle==0)
        subplot(2,3,6);
        for isigsize=1:1 %for compatibility
            rj=results{isp};
            rjs=rj.shuffle;
            %
            hd=[];
            hdata=hist(rjs.ss_class./rjs.ss_total,hist_centers);
            hh=bar(hist_centers,hdata/sum(hdata));hold on;
            hd(1)=hh;
            set(hh,'FaceColor',[1 1 1]);
            set(gca,'XLim',[0 1]);
            p_shuffle=sum(rjs.ss_class./rjs.ss_total<rj.ss_class/rj.ss_total)/length(rjs.ss_class);
            hstring=sprintf('p=%6.4f (shuffle)',p_shuffle);
            %
            if (isfield(rj,'segflip'));
                rjs=rj.segflip;
                hdata=hist(rjs.ss_class./rjs.ss_total,hist_centers);
                hh=bar(hist_centers,hdata/sum(hdata));hold on;
                set(hh,'FaceColor',[1 1 1]);
                set(hh,'LineStyle',':');
                hd(2)=hh;
                p_segflip=sum(rjs.ss_class./rjs.ss_total<rj.ss_class/rj.ss_total)/length(rjs.ss_class);
                hstring=strvcat(hstring,sprintf('p=%6.4f (segflip)',p_segflip));
            end
            hd(end+1)=plot(rj.ss_class/rj.ss_total*[1 1],get(gca,'YLim'),'k-','LineWidth',2);hold on;
            hstring=strvcat(hstring,'data');
            xlabel('class/total');
            ylabel('fraction');
            %
            hl=legend(hd,hstring,'Location','NorthWest');
            set(hl,'FontSize',7);
            title('variance explained','FontSize',7);
        end
        if ~isempty(opts.classifiers)
            subplot(2,3,5);
            frac_corr=zeros(length(opts.classifiers),length(sc.sigsize),2);
            lstring=[];
            ifflip=isfield(rj,'segflip');
            for ish=1:2+ifflip %allow for plotting of nflips
                for ic=1:length(opts.classifiers)
                    if (ish==1) lstring=strvcat(lstring,cat(2,opts.classifiers{ic},' raw')); end
                    if (ish==2) lstring=strvcat(lstring,cat(2,opts.classifiers{ic},' shuffled')); end
                    if (ish==3) lstring=strvcat(lstring,cat(2,opts.classifiers{ic},' segflip')); end
                    for isigsize=1:1 %for compatibility
                        if (ish==1)
                            rj=results{isp};
                            cmat=getfield(rj,cat(2,'fc_',opts.classifiers{ic},'_cmatrix'));
                            frac_corr(ic,1,ish)=sum(diag(cmat))/sum(sc.nsamps);
                        end
                        if (ish==2)
                            rjs=rj.shuffle;
                            cmat_shuf=getfield(rjs,cat(2,'fc_',opts.classifiers{ic},'_cmatrix'));
                            frac_corr(ic,1,ish)=sum(diag(sum(cmat_shuf,3)))/sum(sc.nsamps)/size(cmat_shuf,3);
                        end
                        if (ish==3)
                            rjs=rj.segflip;
                            cmat_flip=getfield(rjs,cat(2,'fc_',opts.classifiers{ic},'_cmatrix'));
                            frac_corr(ic,1,ish)=sum(diag(sum(cmat_flip,3)))/sum(sc.nsamps)/size(cmat_flip,3);
                        end
                    end %isigsize
                end %ic
            end %ish
            plot([0 1],frac_corr(:,[1 1],1),'-','LineWidth',2); hold on;
            plot([0 1],frac_corr(:,[1 1],2),'-.','LineWidth',2); hold on;
            if (ifflip)
                plot([0 1],frac_corr(:,[1 1],3),':','LineWidth',2); hold on;
            end
            plot([0 1],[0.5 0.5],'k-');
            set(gca,'XTick',[0 1]);
            set(gca,'XTickLabel',' ');
            set(gca,'XLim',[0 1]);
            set(gca,'YLim',[0 1]);
            ylabel('frac corr');
            hl=legend(lstring,'Location','SouthEast');
            set(hl,'FontSize',7);
            title('raw vs. shuffled','FontSize',7);
        end %shuffle and classifiers
        clear r rjs  hh hdata hd cmat cmat_shuf cmat_flip hl ish isp ic ifflip isigsize hstring
        clear p_shuffle p_segflip
    end %shuffle
end %is (scenario)
clear rj nr nc ic ig is ifeat ixv frac_corr lstring labels hp ipart tstring

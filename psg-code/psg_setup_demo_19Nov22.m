%psg_setup_demo
% demonstrate setup for perceptual space geometry
%
% file writing will need to be changed to make use of actual stimulus file
% names (and instances of them) -- here, 'type01', etc, is used
%
% to do:
%  * fix colormap of context map: black ref, white for not used, other
%  colors for stimuli tested in multiple contexts
%  * make analysis of context data a function, which includes statistical table of context data
%
%  18Nov22: add option for frozen random numbers
%  18Nov22: add setseq: options for stimulus sets across sessions
%  19Nov22: if multiplicative scheme requested, only redo the sessions if the number of sessions needs to be adjusted
%  19Nov22: allow for a list of actions (plots, stats, write files)
%
%   See also:  PSG_SESSCONFIG_MAKE, PSG_COND_CREATE, PSG_COND_WRITE, PSG_COND_CREATE, NICESUBP, ZPAD
%    PSG_SESSION_STATS, PSG_TRIAD_STATS, PSG_SPOKES_SETUP, PRIMROOT, PSG_SESS_PERM.
%
if ~exist('opts_psg')
    opts_psg=struct;
end
opts_psg=psg_defopts(opts_psg);
%
if ~exist('zpad_typeno') zpad_typeno=2; end %padding for type number (unique ID for stimulus, 1 to nstims)
%
opts_psg=filldefault(opts_psg,'if_log',1);
if ~exist('psg_alt_setups') %alternative setups
    psg_alt_setups{1}.nstims=37;
    psg_alt_setups{1}.ncompares=8;
    psg_alt_setups{1}.novlp=2;
    psg_alt_setups{1}.nsess=5;
    %
    psg_alt_setups{2}.nstims=25;
    psg_alt_setups{2}.ncompares=8;
    psg_alt_setups{2}.novlp=2;
    psg_alt_setups{2}.nsess=10;
    %
    psg_alt_setups{3}.nstims=19;
    psg_alt_setups{3}.ncompares=8;
    psg_alt_setups{3}.novlp=2;
    psg_alt_setups{3}.nsess=16;
    %
    psg_alt_setups{4}.nstims=29;
    psg_alt_setups{4}.ncompares=10;
    psg_alt_setups{4}.novlp=3;
    psg_alt_setups{4}.nsess=8;
end
disp(sprintf('%1.0f->setup with %3.0f stimuli, %3.0f comparison stimuli per trial, overlap %3.0f; %3.0f sessions (default)',...
    0,opts_psg.cond_nstims,opts_psg.cond_ncompares,opts_psg.cond_novlp,opts_psg.cond_nsess));
for k=1:length(psg_alt_setups)
    disp(sprintf('%1.0f->setup with %3.0f stimuli, %3.0f comparison stimuli per trial, overlap %3.0f',...
        k,psg_alt_setups{k}.nstims,psg_alt_setups{k}.ncompares,psg_alt_setups{k}.novlp));
end
k=getinp('choice (negative to modify)','d',length(psg_alt_setups)*[-1 1],0);
if k~=0
    opts_psg.cond_nstims=psg_alt_setups{abs(k)}.nstims;
    opts_psg.cond_ncompares=psg_alt_setups{abs(k)}.ncompares;
    opts_psg.cond_novlp=psg_alt_setups{abs(k)}.novlp;
    opts_psg.cond_nsess=psg_alt_setups{abs(k)}.nsess;
end
if k<0
    opts_psg.cond_nstims=getinp('nstims','d',[1 1000],opts_psg.cond_nstims);
    opts_psg.cond_ncompares=getinp('ncompares','d',[1 1000],opts_psg.cond_ncompares);
    opts_psg.cond_novlp=getinp('novlp','d',[1 1000],opts_psg.cond_novlp);
end
%
opts_psg.cond_nsess=getinp('number of sessions','d',[1 100],opts_psg.cond_nsess);
%
for k=1:length(opts_psg.refseq_labels)
    disp(sprintf('%1.0f->method for choosing stimuli in overlap: %s',k,opts_psg.refseq_labels{k}));
end
opts_psg.refseq=getinp('choice','d',[1 length(opts_psg.refseq_labels)],opts_psg.refseq);
%
for k=1:length(opts_psg.setseq_labels)
    disp(sprintf('%1.0f->method for sequencing stimulus sets across sessions: %s',k,opts_psg.setseq_labels{k}));
end
opts_psg.setseq=getinp('choice','d',[1 length(opts_psg.setseq_labels)],opts_psg.setseq);
%
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1]);
%
if (if_frozen~=0)
    rng('default');
    if (if_frozen<0)
        rand(1,abs(if_frozen));
    end
else
    rng('shuffle');
end
%
[sessions,sessions_sorted,opts_used,psg_desc]=psg_sessconfig_make(opts_psg);
opts_psg.cond_desc=psg_desc;
%
disp(sprintf('Analyzing setup with %s',opts_psg.cond_desc));
ntrials=size(sessions,1);
%
%accumulate and display statistics of the configuration
%
stats=psg_session_stats(sessions,setfield(opts_psg,'if_log',1));
%
%do we want to permute the stimulus number via a multiplicative scheme?
%
[primroots,nrp_phi]=primroot(opts_psg.cond_nstims);
if length(primroots)==0
    disp(sprintf('no primitive roots mod %3.0f, so no option for a multiplicative scheme',opts_psg.cond_nstims));
else
    disp(sprintf('primitive roots mod %3.0f:',opts_psg.cond_nstims));
    disp(primroots);
    if_mult=getinp('1 to apply a multiplicative scheme','d',[0 1]);
    if (if_mult)
        disp('available primitive roots')
        disp(primroots);
        prim=NaN;
        while (~ismember(prim,primroots))
            prim=getinp('a primitive root','d',primroots([1 end]),primroots(1));
        end
        disp(sprintf('period is %3.0f',nrp_phi));
        if (nrp_phi~=opts_psg.cond_nsess)
            if_change_nsess=getinp(sprintf('1 to change to have %3.0f sessions instead of %3.0f',nrp_phi,opts_psg.cond_nsess),'d',[0 1]);
        else
            if_change_nsess=0;
        end
        if if_change_nsess
            opts_psg.cond_nsess=nrp_phi;
            [sessions,sessions_sorted,opts_used,psg_desc]=psg_sessconfig_make(opts_psg);
            opts_psg.cond_desc=psg_desc;
            stats=psg_session_stats(sessions,setfield(opts_psg,'if_log',1)); %redo statistics
        end
        %apply permutation based on primitive roots
        pmode=struct;
        pmode.type='primroot';
        pmode.prim=prim;
        [sessions,sessions_sorted,opts_used,psg_desc,warnings,perms_used]=...
            psg_sess_perm(sessions,sessions_sorted,pmode,setfield(opts_psg,'if_log',1));
        opts_psg.cond_desc=psg_desc;
        stats=psg_session_stats(sessions,setfield(opts_psg,'if_log',1));
    end
end
%
% plot and save if requested
%
choice=0;
sstrings={'session order','sorted order'};
ystrings={'trial in session','sorted trial'};
while (choice>=0)
    disp('-1->done')
    disp(' 0->write session files')
    disp(' 1->show heatmap of conditions: which stimuli are where on each trial');
    disp(' 2->show heatmap of conditions: roles of stimuli on each trial');
    disp(' 3->show heatmap of conditions: roles of stimuli w.r.t. context');
    disp(' 4->compute and show triad statistics');
    disp(' 5->compute and show triad statistics, including cumulative statistics');
    choices=getinp('one or more choices (-1 at end to quit)','d',[-1 5]);
    for ichoice=1:length(choices)
        choice=choices(ichoice);
        switch choice
            case 0 %write cond files for each session
                filename_base=getinp('file name base (and path), _sess[#].csv will be appended','s',[]);
                disp('options for stimulus example re-use')
                for k=1:length(opts_psg.example_infix_labels)
                    disp(sprintf('%1.0f->%s',k,opts_psg.example_infix_labels{k}));
                end
                opts_psg.example_infix_mode=getinp('mode','d',[1 length(opts_psg.example_infix_labels)],opts_psg.example_infix_mode);
                %
                typenames=cell(opts_psg.cond_nstims,1);
                for istim=1:opts_psg.cond_nstims
                    typenames{istim}=cat(2,'type',zpad(istim,zpad_typeno));
                end
                [session_cells,perms_used,examps_used]=psg_cond_create(sessions,typenames,opts_psg);
                for isess=1:opts_psg.cond_nsess
                    filename=cat(2,filename_base,'_sess',zpad(isess,opts_psg.sess_zpad));
                    psg_cond_write(filename,session_cells{isess},setfield(opts_psg,'if_log',1));
                end
            case {4,5} %compute and show triad statistics
                if (choice==4)
                    if_cumulative=0;
                else
                    if_cumulative=1;
                end
                triad_stats=psg_triad_stats(sessions,setfields(opts_psg,{'if_log','if_cumulative'},{1,if_cumulative}));
            case 1 %show which stimuli are where on each trial, for session-order and sorted-order
                for iorder=1:2
                    switch iorder
                        case 1
                            sess_use=sessions;
                        case 2
                            sess_use=sessions_sorted;
                    end
                    tstring=cat(2,'config, ',sstrings{iorder},': ',opts_psg.cond_desc,sprintf(' (%3.0f trials)',ntrials));
                    figure;
                    set(gcf,'Position',[50 50 1200 750]);
                    set(gcf,'NumberTitle','off');
                    set(gcf,'Name',tstring);
                    [nr,nc]=nicesubp(opts_psg.cond_nsess,0.7);
                    for isess=1:opts_psg.cond_nsess
                        subplot(nr,nc,isess);
                        imagesc(sess_use(:,:,isess));
                        ylabel(ystrings{iorder});
                        set(gca,'YTick',[1 ntrials]);
                        set(gca,'YTickLabel',[1 ntrials]);
                        set(gca,'XTick',[1:opts_psg.cond_ncompares+1]);
                        set(gca,'XTickLabel',strvcat('ref',num2str([1:opts_psg.cond_ncompares]')));
                        title(sprintf('session %2.0f',isess));
                    end
                    colormap jet;
                    axes('Position',[0.02,0.02,0.01,0.01]); %for text
                    text(0,0,tstring,'Interpreter','none');
                    axis off;
                end % iorder
            case 2 %show roles of stimuli on each trial, for session-order and sorted-order
                for iorder=1:2
                    switch iorder
                        case 1
                            sess_use=sessions;
                        case 2
                            sess_use=sessions_sorted;
                    end
                    tstring=cat(2,'roles, ',sstrings{iorder},': ',opts_psg.cond_desc,sprintf(' (%3.0f trials)',ntrials));
                    figure;
                    set(gcf,'Position',[50 50 1200 750]);
                    set(gcf,'NumberTitle','off');
                    set(gcf,'Name',tstring);
                    [nr,nc]=nicesubp(opts_psg.cond_nsess,0.7);
                    for isess=1:opts_psg.cond_nsess
                        roles=zeros(ntrials,opts_psg.cond_nstims);
                        for itrial=1:ntrials
                            roles(itrial,sess_use(itrial,1,isess))=1; % reference
                            roles(itrial,sess_use(itrial,1+[1:opts_psg.cond_ncompares],isess))=2; %comparison       
                        end
                        subplot(nr,nc,isess);
                        imagesc(roles);
                        ylabel(ystrings{iorder});
                        set(gca,'YTick',[1 ntrials]);
                        set(gca,'YTickLabel',[1 ntrials]);
                        xlabel('stimuli')
                        set(gca,'XTick',[1 opts_psg.cond_nstims]);
                    end
                    colormap hot;
                    axes('Position',[0.02,0.02,0.01,0.01]); %for text
                    text(0,0,tstring,'Interpreter','none');
                    axis off;
                end %iorder
            case 3 %show roles of stimuli on each trial w.r.t. context
                %color= role (comparison, each context, or not used for context
                %since the first step is to pull out the trials for each stimulus as reference
                %the diagonal (the reference stimulus) doesn't depend on trial order vs sorted order
                %off-diagonal elements may change their contexts (i.e., color tag) but should retain their groupings by color
                for iorder=1:2
                    switch iorder
                        case 1
                            sess_use=sessions;
                        case 2
                            sess_use=sessions_sorted;
                    end
                    tstring=cat(2,'contexts, ',sstrings{iorder},': ',opts_psg.cond_desc,sprintf(' (%3.0f trials)',ntrials));
                    figure;
                    set(gcf,'Position',[50 50 1200 750]);
                    set(gcf,'NumberTitle','off');
                    set(gcf,'Name',tstring);
                    [nr,nc]=nicesubp(opts_psg.cond_nsess,0.7);
                    for isess=1:opts_psg.cond_nsess
                        indisj=0;% flag for nondisjoint contexts
                        contexts=zeros(opts_psg.cond_nstims,opts_psg.cond_nstims);
                        for istim=1:opts_psg.cond_nstims
                            contexts(istim,istim)=-1; %flag for reference
                            trials=find(sess_use(:,1,isess)==istim);
                            configs=sess_use(trials,1+[1:opts_psg.cond_ncompares],isess); %configurations of comparison stimuli for this reference stimulus
                            %assemble an array: first two columns are the sorted pairs of comparison stimuli; third column is the trial in configs
                            npairs=nchoosek(opts_psg.cond_ncompares,2);
                            nconfigs=size(configs,1);
                            pair_dict=zeros(npairs*nconfigs,3);
                            for iconfig=1:nconfigs
                                pair_dict((iconfig-1)*npairs+[1:npairs],:)=[nchoosek(sort(configs(iconfig,:)),2),repmat(iconfig,[npairs,1])];
                            end
                            pair_dict=sortrows(pair_dict); 
                            %now find sets of two or more rows of pair_dict that agree in first two elements
                            %and create an array (blocks), in which blocks(:,[1 2]) are the unique pair, and blocks(:,2+[1:nconfigs] are the trials in which they occur
                            upairs=unique(pair_dict(:,[1 2]),'rows'); %unique pairs
                            blocks=zeros(0,2+nconfigs); %first two cols are the pair, rest are the trials that include it
                            for iup=1:size(upairs,1)
                                matches=find(all(upairs(iup,:)==pair_dict(:,[1 2]),2));
                                if length(matches)>1
                                    newblock=[upairs(iup,:) pair_dict(matches,3)'];
                                    blocks(end+1,1:2+length(matches))=newblock;
                                end
                            end
                            blocks=sortrows(blocks,2+[1:nconfigs]); %now blocks are sorted by unique sets of contexts (trials) in which they occur
                            ucontexts=unique(blocks(:,2+[1:nconfigs]),'rows');
                            %for each unique set of contexts that contain a common pair, add a tag to the stimuli in contexts(istim,:) for that pair
                            for iuc=1:size(ucontexts,1)
                                pairs_involved=find(all(blocks(:,2+[1:nconfigs])==ucontexts(iuc,:),2));
                                for ip=1:length(pairs_involved)
                                    stims_involved=blocks(pairs_involved(ip),[1 2]);
                                    if any(contexts(istim,stims_involved)>0 & contexts(istim,stims_involved)~=iuc)
                                        if (indisj==0)
                                            disp(sprintf('nondisjoint contexts found, session %2.0f reference stimulus %3.0f',isess,istim));
                                            indisj=1;
                                        end
                                    else
                                        contexts(istim,stims_involved)=iuc;
                                    end %nondisjoint contexts?
                                end %pair
                            end %iuc: unique context
                        end %istim
                        subplot(nr,nc,isess);
                        imagesc(contexts);
                        ylabel(cat(2,'stimuli (',ystrings{iorder},')'));
                        set(gca,'YTick',[1 opts_psg.cond_nstims]);
                        set(gca,'YTickLabel',[1 opts_psg.cond_nstims]);
                        xlabel('stimuli')
                        set(gca,'XTick',[1 opts_psg.cond_nstims]);
                        axis square;
                    end
                    colormap hot;
                    axes('Position',[0.02,0.02,0.01,0.01]); %for text
                    text(0,0,tstring,'Interpreter','none');
                    axis off;
                end %iorder
            otherwise
                disp('done.');
        end %switch
    end %ichoice
end %choice>=0

function opts_psg_used=psg_sess_process(sessions,sessions_sorted,opts_psg) 
%
% opts_psg_used=psg_sess_process(sessions,sessions_sorted,opts_psg) is
% menu-driven routine to write stimulus files, plot, or calculate statistics for psg package
%
% opts_psg: options structure, see psg_defopts, can be omitted
%    if omitted:
%      cond_ncompares, cond_nsess, cond_nstims are determined from sess
%      cond_novlp, refseq, setseq set to NaN,
%      cond_desc is created (typically created in psg_sessconfig_make)
% sessions: dim 1 is trial, dim 2 is stimulus number (1 to nstims), dim 3 is session 
% sessions_sorted: as in sessions, but trials are ordered by reference stimulus
%     and the positions of the comparison stimuli are sequential
%
% Note that if cond files are creaed, the type names are generic (type01, typeo2, ...); 
%   this is in contreast to psg_spokes_setup, which calls PSG_COND_CREATE uses with custom names 
%   that specify the btc coordinates
%
% opts_psg_used: opts_psg, with unused options filled by default
%
%   See also:  PSG_SESSCONFIG_MAKE, PSG_COND_CREATE, PSG_COND_WRITE,NICESUBP, ZPAD
%    PSG_SESSION_STATS, PSG_TRIAD_STATS, PSG_SPOKES_SETUP, PRIMROOT, PSG_SESS_PERM, PSG_SESS_DEMO,
%    PSG_UMI_STATS, PSG_QUAD_STATS, PSG_TENT_STATS.
%
undets={'cond_novlp','refseq','setseq'};
if (nargin<=2)
    opts_psg=struct;
end
%fill in values from sessions
ntrials=size(sessions,1);
ncompares=size(sessions,2)-1;
nsess=size(sessions,3);
nstims=length(unique(sessions(:)));
%
if ~isfield(opts_psg,'cond_ncompares')
    opts_psg.cond_ncompares=ncompares;
end
if ~isfield(opts_psg,'cond_nsess')
    opts_psg.cond_nsess=nsess;
end
if ~isfield(opts_psg,'cond_nstims')
    opts_psg.cond_nstims=nstims;
end
if opts_psg.cond_ncompares~=ncompares
    warning('mismatch between number of comparison stimuli specified in options (%3.0f) and found in session array (%3.0f)',opts_psg.cond_ncompares,ncompares);
end
if opts_psg.cond_nsess~=nsess
    warning('mismatch between number of sessions specified in options (%3.0f) and found in session array (%3.0f)',opts_psg.cond_nsess,nsess);
end
if opts_psg.cond_nstims~=nstims
    warning('mismatch between number of stimuli specified in options (%3.0f) and found in session array (%3.0f)',opts_psg.cond_nstims,nstims);
end
opts_psg=filldefault(opts_psg,'cond_desc',sprintf('nstims=%2.0f, ncompares=%2.0f, nsess=%2.0f',nstims,ncompares,nsess));
opts_psg=filldefault(opts_psg,'if_log',1);
%
for iundet=1:length(undets)
    opts_psg=filldefault(opts_psg,undets{iundet},NaN);
end
opts_psg=psg_defopts(opts_psg);  %fill in other values
opts_psg_used=opts_psg;
%
% write files, make plots, and calculate statistics
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
    disp(' 6->compute and show ultrametric statistics');
    disp(' 7->compute and show ultrametric statistics, including cumulative statistics');
    disp(' 8->compute and show quad statistics');
    disp(' 9->compute and show quad statistics, including cumulative statistics')
    disp('10->compute and show tent statistics');
    disp('11->compute and show tent statistics, including cumulative statistics');
;
%
    choices=getinp('one or more choices (-1 at end to quit)','d',[-1 11]);
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
                    typenames{istim}=cat(2,'type',zpad(istim,opts_psg.typeno_zpad));
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
            case {6,7} %compute and show ultrametric statistics
                if (choice==6)
                    if_cumulative=0;
                else
                    if_cumulative=1;
                end
                umi_stats=psg_umi_stats(sessions,setfields(opts_psg,{'if_log','if_cumulative'},{1,if_cumulative}));
            case {8,9} %compute and show quad statistics
                if (choice==8)
                    if_cumulative=0;
                else
                    if_cumulative=1;
                end
                quad_stats=psg_quad_stats(sessions,setfields(opts_psg,{'if_log','if_cumulative'},{1,if_cumulative}));
            case {10,11} %compute and show quad statistics
                if (choice==10)
                    if_cumulative=0;
                else
                    if_cumulative=1;
                end
                tent_stats=psg_tent_stats(sessions,setfields(opts_psg,{'if_log','if_cumulative'},{1,if_cumulative}));
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
return
end

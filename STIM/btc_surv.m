% btc_surv
% surveys which combinations of parameters have been implemented in btc_augcoords
%
% creates a structure pmethsum, that summarizes pairwise methods,
% indicating available methods and extreme values that can be reached
%
%    See also: BTC_TEST, BTC_DEFINE, BTC_CORRS2VEC, BTC_VEC2CORRS, BTC_VEC2LETCODE, BTC_LETCODE2VEC,
%      BTC_EXPTNAME, BTC_COORKINDS, BTC_AUGCOORDS.
%
if ~exist('opts') 
    opts=[];
end
opts=filldefault(opts,'ifshow',1);
opts=filldefault(opts,'ifshow_fig',0);
opts=filldefault(opts,'ordernames',{'gamma','beta','theta','alpha'});
%
aug_opts=[];
%
disp(opts)
disp(' exit any time and re-define opts if desired.');
%
dict=btc_define(opts);
dict_def=btc_define();
%
nparlist=getinp('number of parameters, +n  or range for all, -n to specify parameter sets','d',[-10 3],[0 1 2]);
if (any(nparlist<0))
    npar=-nparlist(1);
    nlist=getinp('number of sets of parameters','d',[1 nchoosek(length(dict.codel),npar)],1);
    coordnums_list=zeros(nlist,-nparlist(1));
    for icomb=1:size(coordnums_list,1)
        coordnums_list(icomb,:)=getinp(sprintf('which %2.0f coordinate(s) to manipulate for set %5.0f',npar,icomb),'d',[1 length(dict.checks)]);
    end
    clear nlist
end
ifprobe=getinp('1 to probe pairwise methods for range of permissible values','d',[0 1]);
if (ifprobe>0)
    dp=getinp('resolution of probe (.3 just to check for Pickard violations, 0.05 standard res, 0.02 fine)','f',[0 .5],0.05);
    dpend=getinp('fractional reduction of probe at extremes to avoid "pathologic" correls=1','f',[0 0.1],0.01);
    paramvals=[-1:dp:1];
    paramvals=paramvals-mean(paramvals);
    paramvals(1)=paramvals(1)+dpend*dp;
    paramvals(end)=paramvals(end)-dpend*dp;
    pmethsum=[];
    %
    ifsmin=getinp('1 for simple exploration of entropy in neighborhood of extremes of permissible values (2 for detailed output)','d',[0 2]);
    dpfrac=getinp('fraction of probe step size','f',[0 1],dpend);
    ddp=dp*dpfrac;
end
disp(' experiment   method(s).....');
for nparptr=1:length(nparlist)
    if (nparlist(1)>=0)
        npar=nparlist(nparptr);
        coordnums_list=nchoosek([1:length(dict.codel)],npar);
    end
    for icomb=1:size(coordnums_list,1)
        spec=[];
        for ivar=coordnums_list(icomb,:)
            spec=setfield(spec,dict.codel(ivar),0);
        end
        if (isempty(spec))
            spec=struct([]);
        end
        results=btc_augcoords(spec,dict,aug_opts);
        outstring=sprintf('%12s',results.exptname);
        for im=1:length(results.method)
            outstring=cat(2,outstring,'   ',sprintf('%s %s',results.method{im}.name,results.method{im}.variant_lab));
        end
        disp(outstring);
        %
        % sweep param1 and param2, and see which ones return a method in which
        % probabilities are normalized (ok_norm==1) and >=0 (ok_probs==1),
        % and plot the results
        %
        nmethods=length(results.method);
        if (npar==2) & (ifprobe==1) & (nmethods>0)
            method_bare=cell(1,nmethods);
            for im=1:nmethods
                method_bare{im}.name=results.method{im}.name;
                method_bare{im}.variant_num=results.method{im}.variant_num;
                method_bare{im}.variant_lab=results.method{im}.variant_lab;
            end
            for iv=1:2
                icn(iv)=dict.codel(coordnums_list(icomb,iv));
            end
            ic=char(icn);
            %
            pmethsum=setfield(pmethsum,results.exptname,method_bare);
            entropies=repmat(-Inf,[length(paramvals),length(paramvals),length(method_bare)]);
            extremes=zeros(8,2,nmethods); %extreme values on x axis, in quadrant 1, on y axis, in quadrant 2, etc 
            for p1val=1:length(paramvals)
                xy(1,1)=paramvals(p1val);
                for p2val=1:length(paramvals)
                    xy(1,2)=paramvals(p2val);
                    specv=spec;
                    for iv=1:2
                        specv=setfield(specv,ic(iv),xy(iv));
                    end
                    rv=btc_augcoords(specv,dict,aug_opts);
                    for im=1:nmethods
                        if (rv.method{im}.ok_norm==1) & (rv.method{im}.ok_probs==1)
                            entropies(p1val,p2val,im)=rv.method{im}.entropy_mrf;
                            if (xy(1)>extremes(1,1,im) & abs(xy(2))<=dp/2) extremes(1,1,im)=xy(1); end
                            if (xy(2)>extremes(3,2,im) & abs(xy(1))<=dp/2) extremes(3,2,im)=xy(2); end
                            if (xy(1)<extremes(5,1,im) & abs(xy(2))<=dp/2) extremes(5,1,im)=xy(1); end
                            if (xy(2)<extremes(7,2,im) & abs(xy(1))<=dp/2) extremes(7,2,im)=xy(2); end
                            r2=sum(xy.^2);
                            if (xy(1)>dp/4 & xy(2)>dp/4 & (abs(xy(1)-xy(2))<dp/2) & r2>sum(extremes(2,:,im).^2)) extremes(2,:,im)=xy; end
                            if (xy(1)<(-dp/4) & xy(2)>dp/4 & (abs(xy(1)+xy(2))<dp/2) & r2>sum(extremes(4,:,im).^2)) extremes(4,:,im)=xy; end
                            if (xy(1)<(-dp/4) & xy(2)<(-dp/4) & (abs(xy(1)-xy(2))<dp/2) & r2>sum(extremes(6,:,im).^2)) extremes(6,:,im)=xy; end
                            if (xy(1)>dp/4 & xy(2)<(-dp/4) & (abs(xy(1)+xy(2))<dp/2) & r2>sum(extremes(8,:,im).^2)) extremes(8,:,im)=xy; end
                        end
                    end
                end %p2val
            end %p1val
            clear xy r2
            for im=1:nmethods
                method_bare{im}.extremes=extremes(:,:,im);
            end
            if (ifsmin>0) %simple exploration, standard Pickard methods only
                for im=1:nmethods
                    entropy_base=zeros(9,1); %origin and 8 extremes
                    entropy_pert=-ones(9,length(dict.codel),2); % origin and 8 extremes, variable, isign
                    vec_base=zeros(9,length(dict.codel));
                    disp(sprintf(' simple exploration of entropy in neighborhood of extrema for method %s %2.0f %s with delta coord of %8.5f',...,
                        method_bare{im}.name,method_bare{im}.variant_num,method_bare{im}.variant_lab,ddp));
                    for iex=0:8
                        specv=spec;
                        if (iex>0)
                            for iv=1:2
                                specv=setfield(specv,ic(iv),extremes(iex,iv,im));
                            end
                        end
                        rm_base=btc_augcoords(specv,dict,aug_opts);
                        vec_base(1+iex,:)=rm_base.method{im}.vec;
                        entropy_base(1+iex)=rm_base.method{im}.entropy_mrf;
                        msg_start=sprintf('  extr pt %2.0f (%1s=%8.5f %1s=%8.5f): entropy %9.6f',...
                            iex,ic(1),getfield(specv,ic(1)),ic(2),getfield(specv,ic(2)),entropy_base(1+iex));
                        if (ifsmin>1) disp(msg_start); end
                        ivelist=[1:length(dict.codel)];
                        ivelist=setdiff(ivelist,coordnums_list(icomb,:));
                        if (strcmp(method_bare{im}.name,'Pickard'))
                            %remove Pickard values also
                            %Pickard relations:
                            picklist=find(dict.inpickard(method_bare{im}.variant_num,:)==1);
                            %Pickard 1: conditional independence in NW and SE use variables bcetv: 
                            %t=v=gb+gc-ge, gt=gv=g^2+bc-e;
                            %Pickard 2: conditional independence in NW and SE use variables bcetv: 
                            %u=w=gb+gc-gd, gu=gw=g^2+bc-d;
                        else
                            picklist=[];
                        end
                        ivelist=setdiff(ivelist,picklist);
                        for ive=ivelist
                            for isign=[1 2]
                                vec_pert=vec_base(1+iex,:);
                                vec_pert(ive)=vec_pert(ive)+(isign*2-3)*ddp;
                                %disp(vec_pert);
                                corrs_pert=getcorrs_p2x2(getp2x2_corrs(btc_vec2corrs(vec_pert,dict_def)),1); %note no dict to ensure getp2x2 works, and warning
                                % if probs are out of range, then entropy_pert remains at -1
                                if (corrs_pert.ok==1)
                                    entropy_pert(1+iex,ive,isign)=corrs_pert.entropy;
                                end
                            end % isign
                            max_pert=max(entropy_pert(1+iex,ive,:));
                            min_pert=min(entropy_pert(1+iex,ive,:));
                            if (max_pert==-1)
                                msg_expand=sprintf('%2s: no room to move',dict.codel(ive));
                                msg_short=sprintf('%2s->noro',dict.codel(ive));
                            elseif entropy_base(1+iex)>max_pert
                                if (min_pert==-1)
                                    msg_expand=sprintf('%2s: passes, one-sided test',dict.codel(ive));
                                    msg_short=sprintf('%2s->1sok',dict.codel(ive));
                                else
                                    msg_expand=sprintf('%2s: passes, two-sided test',dict.codel(ive));
                                    msg_short=sprintf('%2s->2sok',dict.codel(ive));
                                end
                            else
                                msg_expand=sprintf('%2s: FAILS.  Extreme point: %8.5f, neighbors: %8.5f %8.5f',dict.codel(ive),entropy_base(1+iex),squeeze(entropy_pert(1+iex,ive,:)));
                                msg_short=sprintf('%2s->FAIL',dict.codel(ive));
                            end
                            if (ifsmin>1)
                                disp(msg_expand);
                            else
                                msg_start=cat(2,msg_start,' ',msg_short);
                            end
                        end %ive
                        if (ifsmin==1)
                            disp(msg_start);
                        end
                    end %iex
                    method_bare{im}.entropy_base=entropy_base;
                    method_bare{im}.entropy_pert=entropy_pert;
                    method_bare{im}.vec=vec_base;
                end %im
                clear entropy_base entropy_pert iex ive ivelist max_pert min_pert msg_expand msg_short msg_start
                clear vec_base vec_pert isign rm_base picklist
            end %ifsmin
            pmethsum=setfield(pmethsum,results.exptname,method_bare);
            %
            p1name=cat(2,dict.codel(coordnums_list(icomb,1)),':',dict.name{coordnums_list(icomb,1)});
            p2name=cat(2,dict.codel(coordnums_list(icomb,2)),':',dict.name{coordnums_list(icomb,2)});
            figure;
            set(gcf,'Position',[100 100 1400 700]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',cat(2,results.exptname,' -> ',p1name,' vs ',p2name));
            %make the bottom of the colormap black, and map -Inf to it
            cm=colormap;
            cbottom=-1/(size(cm,1)-1);
            cm(1,:)=0;
            colormap(cm);
            for im=1:nmethods %allow for 4 methods
                subplot(2,4,im); %top row shows domain, bottom row will show examples
                imagesc(paramvals,paramvals,entropies(:,:,im)',[cbottom 1]);
                set(gca,'YDir','Normal');
                xlabel(p1name);
                ylabel(p2name);
                axis square;
                set(gca,'XLim',(1+dp/2)*[-1 1]);
                set(gca,'YLim',(1+dp/2)*[-1 1]);
                title(cat(2,method_bare{im}.name,' ',method_bare{im}.variant_lab));
                hold on;
                hp=plot([0 0],[-1 1],'k-');
                set(hp,'Color',[0.25 0.25 0.25]);
                hp=plot([-1 1],[0 0],'k-');
                set(hp,'Color',[0.25 0.25 0.25]);
            end
            clear hp p1name p2name p1val p2val cm cbottom rv method_bare specv extremes iv ic icn entropies
        end
    end
end
clear spec nmethods results outstring im icomb ivar npar nparlist nparptr coordnums_list

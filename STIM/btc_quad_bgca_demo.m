%btc_quad_bgca_demo: demonstrate statistics and examples of textures
%parameterized by the binary pairwise statistics b,g,c,a
%
% Looks at induced correlations and alternate ways of generating the
% textures. Defines variables for generation of maps and creation of a resource.
%
% In contrast to btc_quad_bcde_demo:
%   *separate scales for the varied parameters, in ratio given by vmax_bgca
%   *induced correls also calculated for gamma
%   *several compensation methods, for each way of splitting bgca into pairs
%      but none are perfect, sincd bc->a and g->a are ignored, and all
%      ignore g,{b,c}->{d,e}
% 
%   See also:  BTC_DEFINE, DONUT_METRO, BTC_MIX2TEX_DEMO, BTC_METRO_DEMO,
% BTC_METRO_MIX_DEMO, DONUT_METRO, BTC_AUGCOORDS, BTC_MAKEMAPS, BTC_QUAD_BCDE_DEMO.
%
% in bgca, max values for g is 0.4, b,c is 0.6, a is 1.0
%
vmax_bgca=struct;
vmax_bgca.g=0.4;
vmax_bgca.b=0.6;
vmax_bgca.c=0.6;
vmax_bgca.a=1.0;
vmax_bgca.ref=0.6; %b and c as reference
%
if ~exist('vals_probed') vals_probed=[0:.05:.5]; end %values of pairwise correls in component textures to look at induced correlations
mix_scenarios={{'ag','bc'},{'bg','ca'},{'cg','ba'},{'ag','bc'},{'bg','ca'},{'cg','ba'}};
if_compensates=[0 0 0 1 2 3]; %whether to compensate for induced correlations
%
dict=btc_define;
codel=dict.codel;
nbtc=length(codel);
%
n4=sum(dict.order==2); %number of binary params
quad_signs=1-2*int2nary([0:2^n4-1]'); %all possible signs
quad_lets='bgca';
quad_strings=cell(1,2^n4);
for iq=1:2^n4
    quad_strings{iq}=[];
    for iql=1:n4
        if quad_signs(iq,iql)>0
            qs='+';
        elseif quad_signs(iq,iql)<0
            qs='-';
        else
            qs=0; %shouldn't happen
        end
        quad_strings{iq}=cat(2,quad_strings{iq},qs,quad_lets(iql),' ');
    end
end
%
params='bcdeag'; %parameters to look at induced correlations
np=length(params);
param_ptrs=zeros(1,np);
for ip=1:np
    param_ptrs(ip)=find(codel==params(ip));
end
%
nmix=length(mix_scenarios);
%
inducor=zeros(length(vals_probed),2^n4,nmix,np);
specs=cell(length(vals_probed),2^n4,nmix);
%
mix_scenario_names=cell(1,nmix);
for imx=1:nmix
    mix_scenario=mix_scenarios{imx};
    ncomponents=length(mix_scenarios{imx});
    mix_scenario_names{imx}=[];
    for ic=1:ncomponents
        mix_scenario_names{imx}=cat(2,mix_scenario_names{imx},'[',mix_scenario{ic},'] ');
    end
    if if_compensates(imx)>0
        mix_scenario_names{imx}=cat(2,mix_scenario_names{imx},sprintf(' compensation, method %1.0f',if_compensates(imx)));
    end
    disp(sprintf('doing mix scenario %1.0f:',imx));
    disp(mix_scenario_names{imx})
    for iq=1:2^n4
        for ivp=1:length(vals_probed)
            quad_val=vals_probed(ivp)*quad_signs(iq,:);
            specs{ivp,iq,imx}=cell(1,ncomponents);
            vecs=zeros(ncomponents,nbtc);
            quad_val_use_unscaled=repmat(quad_val,[ncomponents 1]);
            quad_val_use=quad_val_use_unscaled;
            for ic=1:ncomponents
                for ibtc=1:length(mix_scenario{ic})
                    let=mix_scenario{ic}(ibtc);
                    quad_val_use(ic,find(quad_lets==let))=quad_val_use_unscaled(ic,find(quad_lets==let))*vmax_bgca.(let)/vmax_bgca.ref;
                end
            end
            switch if_compensates(imx)
                case 0
                case 1 %compensate b and c (positions 1 and 3 in bgca) in second set from induced correls of g (position 2 in bgca) in first set
                   quad_val_use(2,[1 3])=quad_val_use(2,[1 3])-quad_val_use(1,2).^2;
                case 2 %compensate a (position 4 in bgca) in second set from induced correls of b (position 1 in bgca) in first set
                   quad_val_use(2,4)=quad_val_use(2,4)-quad_val_use(1,1).^2;
                   quad_val_use(2,3)=quad_val_use(2,3)-quad_val_use(1,2).^2; %compensate c for induced by g
                 case 3 %compensate a (position 4 in bgca) in second set from induced correls of c (position 3 in bgca) in first set
                   quad_val_use(2,4)=quad_val_use(2,4)-quad_val_use(1,3).^2;
                   quad_val_use(2,1)=quad_val_use(2,1)-quad_val_use(1,2).^2; %compensate b for induced by g
            end
            for ic=1:ncomponents
                %create the structures for each component
                s=struct;
                for ibtc=1:length(mix_scenario{ic})
                    let=mix_scenario{ic}(ibtc);
                    s.(let)=quad_val_use(ic,find(quad_lets==let));
                end
                specs{ivp,iq,imx}{ic}=s;
                aug=btc_augcoords(specs{ivp,iq,imx}{ic},dict);
                vecs(ic,:)=aug.method{1}.vec; %statistics of each component
            end
            inducor(ivp,iq,imx,:)=mean(vecs(:,param_ptrs));
        end %value probed
    end
end %imx
%plot
[nrows,ncols]=nicesubp(2^n4);
for imx=1:nmix
    mix_scenario=mix_scenarios{imx};
    nmix=length(mix_scenarios{imx});
    tstring=cat(2,'final corrs for ',mix_scenario_names{imx});
    figure;
    set(gcf,'Position',[50 50 1500 850]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',tstring)
    %
    for iq=1:2^n4
        subplot(nrows,ncols,iq);
        data=reshape(inducor(:,iq,imx,:),[length(vals_probed) np]);
        plot(vals_probed,data,'LineWidth',2);
        h=get(gca,'Children');
        set(h(1),'Marker','s');
        set(h(2),'Marker','<');
        set(h(3),'Marker','>');
        set(h(4),'Marker','o');
        set(h(5),'Marker','x');
        set(h(6),'Marker','.');
        set(gca,'XLim',[min(vals_probed) max(vals_probed)]);
        set(gca,'YLim',[-.5 .5]);
        xlabel('component b');
        legend(params','FontSize',7,'Location','NorthWest');
        title(quad_strings{iq});
    end
    %
    axes('Position',[0.01,0.03,0.01,0.01]); %for text
    text(0,0,tstring,'Interpreter','none');
    axis off;
end

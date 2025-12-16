%btc_quad_bcde_demo: demonstrate statistics and examples of textures
%parameterized by the binary pairwise statistics b,c,d,e
%
% Looks at induced correlations and alternate ways of generating the
% textures. Defines variables for generation of maps and creation of a resource.
% 
%   See also:  BTC_DEFINE, DONUT_METRO, BTC_MIX2TEX_DEMO, BTC_METRO_DEMO,
% BTC_METRO_MIX_DEMO, DONUT_METRO, BTC_AUGCOORDS, BTC_MAKEMAPS, BTC_QUAD_BCDE_MAKE,
% BTC_QUAD_BGCA_DEMO.
%
if ~exist('vals_probed') vals_probed=[0:.05:.5]; end %values in component textures to look at induced correlations
mix_scenarios={{'bc','de'},{'bd','ce'},{'b','c','d','e'},{'bc','de'},{'bd','be','cd','ce'}};
if_compensates=[0 0 0 1 0]; %whether to compensate for induced correlations
%
dict=btc_define;
codel=dict.codel;
nbtc=length(codel);
%
n4=sum(dict.order==2); %number of binary params
quad_signs=1-2*int2nary([0:2^n4-1]'); %all possible signs
quad_lets='bcde';
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
params='bcdea'; %parameters to look at induced correlations
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
    if if_compensates(imx)==1
        mix_scenario_names{imx}=cat(2,mix_scenario_names{imx},' compensated');
    end
    disp(sprintf('doing mix scenario %1.0f:',imx));
    disp(mix_scenario_names{imx})
    for iq=1:2^n4
        for ivp=1:length(vals_probed)
            quad_val=vals_probed(ivp)*quad_signs(iq,:);
            specs{ivp,iq,imx}=cell(1,ncomponents);
            vecs=zeros(ncomponents,nbtc);
            quad_val_use=repmat(quad_val,[ncomponents 1]);
            if if_compensates(imx)==1 %for {bc}{de}, change the second component to compensate for d and e induced by first component
                quad_val_use(2,[3 4])=quad_val_use(2,[3 4])-quad_val_use(1,1)*quad_val_use(1,2);
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
        set(gca,'XLim',[min(vals_probed) max(vals_probed)]);
        set(gca,'YLim',[-.5 .5]);
        legend(params','FontSize',7,'Location','NorthWest');
        title(quad_strings{iq});
    end
    %
    axes('Position',[0.01,0.03,0.01,0.01]); %for text
    text(0,0,tstring,'Interpreter','none');
    axis off;
end

%btc_quad_bcde_make: create examples and a resource for maps 
%parameterized by the binary pairwise statistics tuvw
% 
%   See also:  BTC_DEFINE, DONUT_METRO, BTC_MIX2TEX_DEMO, BTC_METRO_DEMO,
% BTC_METRO_MIX_DEMO, DONUT_METRO, BTC_AUGCOORDS, BTC_MAKEMAPS, BTC_QUAD_BCDE_DEMO
% PSG_SPEC2FILENAME, BTC_LETCODE2VEC, PSG_SPOKES_SETUP, BTC_QUAD_BCDE_MAKE.
%
fn_prompt='file name, typically btc_quad_tuvw*.mat';
%
if ~exist('ndigs') ndigs=4; end %number of digits in file name
mix_scenarios={{'tv','uw'}}; %only one scenario
if_compensates=[0]; %whether to compensate for induced correlations
%
dict=btc_define;
codel=dict.codel;
nbtc=length(codel);
%
opts_stn=struct; %for psg_spec2filename
%
n4=sum(dict.order==2); %number of binary params
quad_signs=1-2*int2nary([0:2^n4-1]'); %all possible signs
quad_lets='tuvw';
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
params='tuvw'; %parameters to look at induced correlations
np=length(params);
param_ptrs=zeros(1,np);
for ip=1:np
    param_ptrs(ip)=find(codel==params(ip));
end
%
nmix=length(mix_scenarios);
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
end
%
%read or make texture examples
%
if_read=getinp('0 to read a texture resource, 1 to create','d',[0 1]);
if if_read==0
    fn=getinp(fn_prompt,'s',[]);
    load(fn,'r');
else
    for imx=1:nmix
        disp(sprintf('%1.0f-> mix scenario %s',imx,mix_scenario_names{imx}))
    end
    mix_choice=getinp('choice','d',[1 nmix]);
    r=struct;
    mix_scenario=mix_scenarios{mix_choice};
    r.mix_scenario=mix_scenario;
    r.if_compensate=if_compensates(mix_choice);
    r.mix_scenario_name=mix_scenario_names{mix_choice};
    ncomponents=length(mix_scenario);
    r.ncomponents=ncomponents;
    tuvw_target=getinp('target value for t,u,v,w in final texture','f',[0 1/ncomponents]);
    r.tuvw_target=tuvw_target;
    if ~exist('size_recur') size_recur=256; end %size of map to generate via recursion (only middle is used for Metropolis)
    if ~exist('size_metro') size_metro=size_recur-24; end %size of map to run Metropolis algorithm on, and for stats
    if ~exist('size_show') size_show=32; end %size of map sample to show by this routine
    size_recur=getinp('size of map to generate via recursion','d',[64 1024],size_recur);
    size_metro=getinp('size of map to run Metropolis donut algorithm on','d',[40 size_recur],min(size_metro,size_recur-24));
    select_metro=[1:size_metro]+round((size_recur-size_metro)/2);
    r.map_size_final=size_metro;
    %
    %setup for donut algorithm
    if exist('auxopts')
        auxopts=btc_auxopts(auxopts);
    else
        auxopts=btc_auxopts;
    end
    auxopts.metro_opts.numiters=getinp('number of Metropolis iterations','d',[0 10^6]);
    auxopts.metro_opts.sampfreq_map=0; %frequency to save sampled maps during Metropolis, also used to calculate statistics
    auxopts.metro_opts.sampfreq_stat=0; %never do stats during Metropolis; here we calculate statistics are based on sampled maps
    auxopts.metro_opts.showfreq_map=0; %frequency to show maps during Metropolis
    auxopts.metro_opts.nf_frac=(1/2^8)*getinp('target flip fraction*256','f',[0 2^8],auxopts.metro_opts.nf_frac*256);
    auxopts.metro_opts.nf_dist=getinp('nf_dist (0->const, 1->unif, 2->Bern, 3->Pois, 4->exp)','d',[0 4],auxopts.metro_opts.nf_dist);
    donut=[];
    donut.name='donut';
    donut.matrix=[1 1 1;1 0 1; 1 1 1];
    donut=glider_addcoords(donut);
    metro_opts=auxopts.metro_opts;
    %
    opts_component=struct;
    opts_component.area=[size_recur,size_recur];
    %
    r.opts_component=opts_component;
    r.opts_metro=metro_opts;
    %
    %these fields should facitate creation of sa datastructure and map files for an experiment
    r.spec_final=cell(2^n4,1); %specify b,c,d,e
    r.btc_specoords=cell(2^n4,1); %has NaNs other than for b,c,d,e (coordinates specified)
    r.btc_augcoords=cell(2^n4,1); %includes a and zeros (actual coords)
    r.filebase=cell(2^n4,1); %name base for map filenames, e.g., 'bp0200cm0200dp0140em0140'
    %
    r.stats_init_ideal=zeros(2^n4,nbtc,ncomponents);
    r.stats_final_ideal=zeros(2^n4,nbtc);
    r.stats_init=zeros(2^n4,nbtc,ncomponents);
    r.stats_final=zeros(2^n4,nbtc,ncomponents);
    r.maps_init=cell(1,2^n4);
    r.maps_final=cell(1,2^n4);
    r.opts_metro_used=cell(1,2^n4);   
    for iq=1:2^n4
        map_components=zeros(size_recur,size_recur,ncomponents);
        %create specifications for the components
        % ncomponents * target in final, then add compensation
        comp_vals=repmat(tuvw_target*ncomponents*quad_signs(iq,:),[ncomponents 1]);
        if if_compensates(mix_choice)
            comp_vals(2,[3 4])=comp_vals(2,[3 4])-comp_vals(1,1)*comp_vals(1,2); %in [bc][de], adjust de in comopnent 2 by the induced contribution from bc
        end
        disp(' ');
        disp(sprintf('creating map type %2.0f (%s) with target value %7.3f, using mixing scenario %s',...
            iq,quad_strings{iq},tuvw_target,mix_scenario_names{mix_choice}));
        for ic=1:ncomponents
            s_string=[];
            s=struct;
            for ibtc=1:length(mix_scenario{ic})
                let=mix_scenario{ic}(ibtc);
                s.(let)=comp_vals(ic,find(quad_lets==let));
                s_string=cat(2,s_string,sprintf('%s=%5.2f ',let,s.(let)));
            end
            aug=btc_augcoords(s,dict);
            r.stats_init_ideal(iq,:,ic)=aug.method{1}.vec;
            map_components(:,:,ic)=btc_makemaps(aug.method{1},opts_component,dict);
            counts=btc_map2counts(map_components(:,:,ic));
            r.stats_init(iq,:,ic)=btc_corrs2vec(getcorrs_p2x2(counts/sum(counts(:))))';
            disp(sprintf('                         component %2.0f: map made for spec :   %s',ic,s_string));
        end
        %
        r.stats_final_ideal(iq,:)=mean(r.stats_init_ideal(iq,:,:),3);
        %
        spec_final=struct;
        for iql=1:length(quad_lets)
            spec_final.(quad_lets(iql))=r.stats_final_ideal(iq,param_ptrs(iql));
        end
        r.spec_final{iq}=spec_final; %specify tuvw
        r.btc_specoords{iq}=btc_letcode2vec(spec_final,dict); %has NaNs other than tuvw 
        r.btc_augcoords{iq}=r.stats_final_ideal; %includes a and zeros
        r.filebase{iq}=psg_spec2filename(spec_final,opts_stn); %name base for map filenames, e.g., 'bp0200cm0200dp0200em0200'
        %
        map_init=map_components(select_metro,select_metro,:);
        r.maps_init{iq}=map_init;
        %
        [map_donut,metro_samp,opts_metro_used]=donut_metro(2,donut,map_init,metro_opts); %mix via Metropolis
        %
        r.maps_final{iq}=map_donut;
        r.opts_metro_used{iq}=opts_metro_used;
        %statistics
        for ic=1:ncomponents
            counts=btc_map2counts(map_donut(:,:,ic));
            r.stats_final(iq,:,ic)=btc_corrs2vec(getcorrs_p2x2(counts/sum(counts(:))))';
        end
    end %iq
    r.stats_final_avg=mean(r.stats_final,3); %average statsitics after mixing
    r.stats_final_max_component_diff=max(r.stats_final,[],3)-min(r.stats_final,[],3); %max dev of any component from another cdomopnent
    r.stats_final_max_dev_ideal=max(abs(r.stats_final-repmat(r.stats_final_ideal,[1 1 ncomponents])),[],3); %max def of any component from ideal
    r.stats_final_avg_dev_ideal=abs(r.stats_final_avg-r.stats_final_ideal); %max def of any component from ideal
    if getinp('write a resource file','d',[0 1],1);
        fn_write=getinp(fn_prompt,'s',[]);
        save(fn_write,'r');
    end
end
%plot
map_showsize=16;
[nrows,ncols]=nicesubp(2^n4);
while(map_showsize>0)
    map_showsize=getinp('size of map to show, 0 to end','d',[0 r.map_size_final],16);
    if map_showsize>0
       start_locs=getinp('start locations','d',[1 r.map_size_final-map_showsize+1],[1 1]);
        if length(start_locs)==1
            start_locs=[start_locs start_locs];
        end
        map_sel_x=[0:map_showsize-1]+start_locs(1);
        map_sel_y=[0:map_showsize-1]+start_locs(2);
        for icomp=1:r.ncomponents
            figure;
            tstring=sprintf('mixing scenario %s, tuvw_target %7.3f, component %1.0f, locs [%4.0f %4.0f]',...
                r.mix_scenario_name,r.tuvw_target,icomp,start_locs);
            set(gcf,'Position',[100 100 1200 800]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',tstring);
            for iq=1:2^n4
                subplot(nrows,ncols,iq);
                imagesc(r.maps_final{iq}(map_sel_x,map_sel_y,icomp),[0 1]);       
                title(r.filebase{iq});
                axis square;
                colormap('gray');
                set(gca,'XTick',[]);
                set(gca,'YTick',[]);
            end
            %
            axes('Position',[0.01,0.03,0.01,0.01]); %for text
            text(0,0,tstring,'Interpreter','none');
            axis off;
        end
    end
end
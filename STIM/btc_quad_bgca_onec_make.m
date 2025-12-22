%btc_quad_bgca_onec_make: create examples and a resource for maps 
%parameterized by the binary pairwise statistics g,b,c,a
%using a single-component strategy, rather than mixing
%
% See btc_quad_onecomponent_notes.docx for details.
% Strategy: make a teture with defined beta-cardinals and alpha,  which is Pickard,
% and make g nonzero by flipping |g| of the white checks if g>0, or |g| of the black checks if g<0.
% This will change the betas and alpha, so the initial betas are chosen with that in mind.
%
%%%%%%ideal values of alpha and gamma not done yet
%
%   See also:  BTC_DEFINE, DONUT_METRO, BTC_MIX2TEX_DEMO, BTC_METRO_DEMO,
% BTC_METRO_MIX_DEMO, DONUT_METRO, BTC_AUGCOORDS, BTC_MAKEMAPS, BTC_QUAD_BCDE_DEMO
% PSG_SPEC2FILENAME, BTC_LETCODE2VEC, PSG_SPOKES_SETUP, BTC_QUAD_BGCA_MAKE, BTC_ALPHARANGE.
%
fn_prompt='file name, typically btc_quad_bgca_onec*.mat';
%
if ~exist('alpharange_tol') alpharange_tol=10^-6; end % tolerance on alpha range
if ~exist('ndigs') ndigs=4; end %number of digits in file name
vmax_bgca=struct;
vmax_bgca.g=0.4;
vmax_bgca.b=0.6;
vmax_bgca.c=0.6;
vmax_bgca.a=1.0;
vmax_bgca.ref=0.6; %b and c as reference
%
dict=btc_define;
codel=dict.codel;
nbtc=length(codel);
%
opts_stn=struct; %for psg_spec2filename
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
params=quad_lets; %parameters that go into spec_final
np=length(params);
param_ptrs=zeros(1,np);
for ip=1:np
    param_ptrs(ip)=find(codel==params(ip));
end
%
%
%read or make texture examples
%
if_read=getinp('0 to read a texture resource, 1 to create','d',[0 1]);
if if_read==0
    fn=getinp(fn_prompt,'s',[]);
    load(fn,'r');
else
    r=struct;
    %for compatibility with btc_quad_bgca_make
    r.mix_scenario=[];
    r.if_compensate=0;
    r.mix_scenario_name='one component';
    ncomponents=1;
    r.ncomponents=ncomponents;
    bc_target=getinp('target value for b and c in final texture','f',[0 1/ncomponents]);
    r.bc_target=bc_target;
    r.aneg_option='n/a';
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
    r.spec_final=cell(2^n4,1); %specify coordinates b,g,c,a
    r.btc_specoords=cell(2^n4,1); %has NaNs other than for coordinates specified
    r.btc_augcoords=cell(2^n4,1); %actual coords after augmentation
    r.filebase=cell(2^n4,1); %name base for map filenames, e.g., 'bp0200cm0200dp0140em0140'
    %
    r.stats_init_ideal=zeros(2^n4,nbtc,ncomponents);
    r.stats_final_ideal=zeros(2^n4,nbtc);
    r.stats_init=zeros(2^n4,nbtc,ncomponents);
    r.stats_final=zeros(2^n4,nbtc,ncomponents);
    r.maps_init=cell(1,2^n4);
    r.maps_final=cell(1,2^n4);
    r.opts_metro_used=cell(1,2^n4);
    r.alpharange_tol=alpharange_tol;
    for iq=1:2^n4
        disp('**********************');
        disp(sprintf('  beginning map type %2.0f',iq));
        map_components=zeros(size_recur,size_recur,ncomponents);
        %create specifications for the components
        % ncomponents * target in final, then add compensation
        comp_vals=zeros(1,length(quad_lets));
        for iql=1:length(quad_lets)
            let=quad_lets(iql);
            comp_vals(1,iql)=bc_target*quad_signs(iq,iql)*vmax_bgca.(let)/vmax_bgca.ref;
        end
        disp(sprintf('creating map type %2.0f (%s) with target b,c value %7.3f, using one-component strategy',...
            iq,quad_strings{iq},bc_target));
        disp('initial component values (bgca)')
        disp(comp_vals)
        %adjust second-order stats since gamma will be incremented
        g=comp_vals(1,find(quad_lets=='g'));
        b_targ=comp_vals(1,find(quad_lets=='b'));
        c_targ=comp_vals(1,find(quad_lets=='c'));
        a_trial=comp_vals(1,find(quad_lets=='a'));
        b_use=(b_targ-g^2)/(1-abs(g))^2;
        c_use=(c_targ-g^2)/(1-abs(g))^2;
        %determine range of possible a
        s_trial=struct; %temporary choice to determine feasible range
        s_trial.b=b_use;
        s_trial.c=c_use;
        aug_trial=btc_augcoords(s_trial,dict);
        p2x2_trial=aug_trial.method{1}.p2x2;
        [amin,amax,p2x2_extremes]=btc_alpharange(p2x2_trial);
        amin=amin+alpharange_tol; %add some tolerance to ensure probs >0
        amax=amax-alpharange_tol;
        if_adj_a=0;
        if a_trial>amax
            if_adj_a=1;
            comp_vals(1,find(quad_lets=='a'))=amax;
            disp(sprintf('in single component, a was %7.3f, feasible range [%7.3f %7.3f], set to %7.3f',...
                a_trial,amin,amax,amax));
        end
        if a_trial<amin
            if_adj_a=1;
            comp_vals(1,find(quad_lets=='a'))=amin;
            disp(sprintf('in single component, a was %7.3f, feasible range [%7.3f %7.3f], set to %7.3f',...
                a_trial,amin,amax,amin));
        end
        if if_adj_a==0
            disp(sprintf('in single component, a was %7.3f, feasible range [%7.3f %7.3f]',a_trial,amin,amax));
        end
        a_use=comp_vals(1,find(quad_lets=='a'));
        s_use=s_trial;
        s_use.a=a_use; %a after restriction to feasible range
        %now make the one-component map
        aug=btc_augcoords(s_use,dict);
        map_components_pregamma=btc_makemaps(aug.method{1},opts_component,dict);
        r.stats_init_ideal(iq,:,1)=aug.method{1}.vec;
        %adjust final ideal stats
        r.stats_final_ideal(iq,1,1)=g;
        r.stats_final_ideal(iq,[2 3 4 5],1)=r.stats_init_ideal(iq,[2 3 4 5],1)*(1-abs(g))^2+g^2;
%        r.stats_final_ideal(iq,[6,7,8 9],1)=NaN; %don't yet know the ideal thetas
%        r.stats_final_ideal(iq,10,1)=NaN; %don't yet know the ideal alphas
        %adjust map
        mask=double(rand(size_recur,size_recur)<abs(g)); %|g| of these are white
        if g>0 %increase gamma by flipping fraction |g| of the black checks to white
            mc_flip=find(map_components_pregamma(:)==0);
            map_components_pregamma(mc_flip)=mask(mc_flip);
        end
        if g<0 %decrease gamma by flipping fraction |g| of the white checks to black
            mc_flip=find(map_components_pregamma(:)==1);
            map_components_pregamma(mc_flip)=1-mask(mc_flip);
        end
        map_components(:,:,1)=map_components_pregamma;
        %
        spec_final=struct;
        for iql=1:length(quad_lets)
            spec_final.(quad_lets(iql))=r.stats_final_ideal(iq,param_ptrs(iql));
        end
        r.spec_final{iq}=spec_final; %specify b,c,g,a
        r.btc_specoords{iq}=btc_letcode2vec(spec_final,dict); %has NaNs other than for b,c,g,a 
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
        %%%%kludge since ideal a not yet calculated
        %%%%remove these lines when value installed
        r.stats_final_ideal(iq,10,ic)=r.stats_final(iq,10,ic);
        r.spec_final{iq}.a=r.stats_final(iq,10,ic);
        r.filebase{iq}=psg_spec2filename(r.spec_final{iq},opts_stn); %name base for map filenames, e.g., 'bp0200cm0200dp0200em0200'
        %%%
    end %iq
    r.stats_final_avg=mean(r.stats_final,3); %average statistics after mixing
    r.stats_final_max_component_diff=max(r.stats_final,[],3)-min(r.stats_final,[],3); %max dev of any component from another cdomopnent
    r.stats_final_max_dev_ideal=max(abs(r.stats_final-repmat(r.stats_final_ideal,[1 1 ncomponents])),[],3); %max def of any component from ideal
    r.stats_final_avg_dev_ideal=abs(r.stats_final_avg-r.stats_final_ideal); %max def of any component from ideal
    if getinp('1 to write a resource file','d',[0 1],1);
        fn_write=getinp(fn_prompt,'s',[]);
        save(fn_write,'r');
    end
end
%plot
map_showsize=16;
[nrows,ncols]=nicesubp(2^n4);
while(map_showsize>0)
    disp('**********Notice:  ideal a is empiric not exact**************');
    map_showsize=getinp('size of map to show, 0 to end','d',[0 r.map_size_final],16);
    if map_showsize>0
       start_locs=getinp('start locations','d',[1 r.map_size_final-map_showsize+1],[1 1]);
        if length(start_locs)==1
            start_locs=[start_locs start_locs];
        end
        map_sel_x=[0:map_showsize-1]+start_locs(1);
        map_sel_y=[0:map_showsize-1]+start_locs(2);
        for icomp=1:1
            figure;
            tstring=sprintf('one-component strategy, bc_target %7.3f, locs [%4.0f %4.0f]',...
                r.bc_target,start_locs);
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
%btc_quad_bgca_make: create examples and a resource for maps 
%parameterized by the binary pairwise statistics g,b,c,a
% 
% In contrast to btc_quad_bcde_make:
%   *separate scales for the varied parameters, in ratio given by vmax_bgca
%   *induced correls also calculated for gamma
%   *several compensation methods, for each way of splitting bgca into pairs
%      but none are perfect, sincd bc->a and g->a are ignored, and all
%      ignore g,{b,c}->{d,e}
%   *params set equal to quad_lets, not the induced-correlation params
%   *also a mix scenario of triplets in each componenent, using btc_alpharange to determine
%    feasible range for a given (b,g) or (c,g); this seems to yield the largest range for a
%
%   See also:  BTC_DEFINE, DONUT_METRO, BTC_MIX2TEX_DEMO, BTC_METRO_DEMO,
% BTC_METRO_MIX_DEMO, DONUT_METRO, BTC_AUGCOORDS, BTC_MAKEMAPS, BTC_QUAD_BCDE_DEMO
% PSG_SPEC2FILENAME, BTC_LETCODE2VEC, PSG_SPOKES_SETUP, BTC_QUAD_BCDE_MAKE,
% BTC_QUAD_BGCA_DEMO, BTC_ALPHARANGE, BTC_QUAD_BGCA_ONEC_DEMO.
%
fn_prompt='file name, typically btc_quad_bgca*.mat';
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
mix_scenarios={{'ag','bc'},{'bg','ca'},{'cg','ba'},{'ag','bc'},{'bg','ca'},{'cg','ba'},{'abg','acg'}};
%{abg}{acg} uses btc_alpharange to find the smallest feasible value for a if its target is negative
% 
a_paired_list={'g','c','b','g','c','b','z'}; %what a is paired with
a_component_list=[1 2 2 1 2 2 NaN];%which component has a
if_compensates=[0 0 0 1 2 3 0]; %whether to compensate for induced correlations
% (always compensate for scenario 7, 'abg' amd 'acg'
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
    if if_compensates(imx)>=1
        mix_scenario_names{imx}=cat(2,mix_scenario_names{imx},sprintf(' compensation, method %1.0f',if_compensates(imx)));
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
    bc_target=getinp('target value for b and c in final texture','f',[0 1/ncomponents]);
    r.bc_target=bc_target;
    a_pair_def=a_paired_list{mix_choice};
    if ismember(a_pair_def,'bc')
        aneg_option=getinp('1 to force negative a to zero, -1 to force negative a to lowest feasible value','d',[-1 1],-1);
    elseif a_pair_def~='z'
        aneg_option=getinp('1 to force negative a to zero','d',[0 1],1);
    else
        aneg_option=0;
    end
    r.aneg_option=aneg_option;
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
        comp_vals_unscaled=repmat(bc_target*ncomponents*quad_signs(iq,:),[ncomponents 1]);
        comp_vals=comp_vals_unscaled;
        for iql=1:length(quad_lets)
            let=quad_lets(iql);
            comp_vals(:,iql)=comp_vals(:,iql)*vmax_bgca.(let)/vmax_bgca.ref;
        end
        if a_paired_list{mix_choice}=='z' %a and g are in both components
            comp_vals(:,find(quad_lets=='a'))=comp_vals(:,find(quad_lets=='a'))/ncomponents;
            comp_vals(:,find(quad_lets=='g'))=comp_vals(:,find(quad_lets=='g'))/ncomponents;
            comp_vals(1,find(quad_lets=='b'))=comp_vals(1,find(quad_lets=='b'))-comp_vals(2,find(quad_lets=='g')).^2;
            comp_vals(2,find(quad_lets=='b'))=0;
            comp_vals(1,find(quad_lets=='c'))=0;
            comp_vals(2,find(quad_lets=='c'))=comp_vals(2,find(quad_lets=='c'))-comp_vals(1,find(quad_lets=='g')).^2;
            %use btc_alpharange to restrict a to a feasible value
            disp('initial component values (bgca) prior to checking feasibility')
            disp(comp_vals)
            for ic=1:ncomponents
                s=struct;
                for ibtc=1:length(mix_scenario{ic})
                    let=mix_scenario{ic}(ibtc);
                    if (let~='a')
                        s.(let)=comp_vals(ic,find(quad_lets==let));
                    end
                end               
                aug_trial=btc_augcoords(s,dict);
                p2x2_trial=aug_trial.method{1}.p2x2;
                [amin,amax,p2x2_extremes]=btc_alpharange(p2x2_trial);
                amin=amin+alpharange_tol; %add some tolerance to ensure probs >0
                amax=amax-alpharange_tol;
                cv_prev=comp_vals(ic,find(quad_lets=='a'));
                if cv_prev>amax
                    comp_vals(ic,find(quad_lets=='a'))=amax;
                    disp(sprintf('for component %1.0f, a was %7.3f, feasible range [%7.3f %7.3f], set to %7.3f',...
                        ic,cv_prev,amin,amax,amax));
                end
                if cv_prev<amin
                    comp_vals(ic,find(quad_lets=='a'))=amin;
                    disp(sprintf('for component %1.0f, a was %7.3f, feasible range [%7.3f %7.3f], set to %7.3f',...
                        ic,cv_prev,amin,amax,amin));
                end
            end %ic
        end %scenario with a and g in both components
        disp(sprintf('creating map type %2.0f (%s) with target b,c value %7.3f, using mixing scenario %s',...
            iq,quad_strings{iq},bc_target,mix_scenario_names{mix_choice}));
        if if_compensates(mix_choice)>0
            disp(' uncompensated spec values for bgca')
            disp(comp_vals)
        end
        switch if_compensates(mix_choice)
            case 0
            case 1 %{'ag','bc'}, compensate b and c (positions 1 and 3 in bgca) in second set from induced correls of g (position 2 in bgca) in first set
               comp_vals(2,[1 3])=comp_vals(2,[1 3])-comp_vals(1,2).^2;
            case 2 %{'bg','ca'}, compensate a (position 4 in bgca) in second set from induced correls of b (position 1 in bgca) in first set
               comp_vals(2,4)=comp_vals(2,4)-comp_vals(1,1).^2;
               comp_vals(2,3)=comp_vals(2,3)-comp_vals(1,2).^2; %compensate c for induced by g
            case 3 %,{'cg','ba'}, compensate a (position 4 in bgca) in second set from induced correls of c (position 3 in bgca) in first set
               comp_vals(2,4)=comp_vals(2,4)-comp_vals(1,3).^2;
               comp_vals(2,1)=comp_vals(2,1)-comp_vals(1,2).^2; %compensate b for induced by g
        end
        disp('component values for bgca');
        disp(comp_vals);
        %
        %process aneg_option
        %repeat calculation of mix_choice, to allow for adaptable choice
        a_component=a_component_list(mix_choice);
        a_pair=a_paired_list{mix_choice};
        let_a_pair=find(quad_lets==a_pair); %points to the letter paired with a
        let_a=find(quad_lets=='a');
         switch aneg_option
            case 0
            case 1
                comp_vals(a_component,let_a)=max(comp_vals(a_component,let_a),0);
                disp('after possible increase of a to 0')
                disp(comp_vals);
            case -1
                %what is a paired with?
                switch a_pair
                    case {'b','c'}
                        a_min=2*abs(comp_vals(a_component,let_a_pair))-1; %alpha >=2*|beta|-1
                    otherwise
                        warning('inconsist aneg_option');
                end
                comp_vals(a_component,let_a)=max(a_min,comp_vals(a_component,let_a));
                disp('after possible increase of a to smallest feasible value')
                disp(comp_vals);
        end
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
        r.spec_final{iq}=spec_final; %specify b,c,g,a
        r.btc_specoords{iq}=btc_letcode2vec(spec_final,dict); %has NaNs other than for b,g,c,a
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
            tstring=sprintf('mixing scenario %s, bc_target %7.3f, component %1.0f, locs [%4.0f %4.0f], aneg option %1.0f',...
                r.mix_scenario_name,r.bc_target,icomp,start_locs,r.aneg_option);
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
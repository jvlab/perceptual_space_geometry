%psg_raystats_summ: summarize ray angle statistics (see psg_visualize_demo for plotting)
%
% tabulates, for multiple datasets
%angles between positive and negative rays on each axis (using psg_rayfit with bid=1)
%angles between best-fitting line for each pair of axes (using psg_rayfit with bid=0)
%multipliers (gains) on each axis
%
% reads choice files to find critical jitter, and then
% includes confidence limits based on critical jitter
%
%  See also: PSG_GET_COORDSETS, PSG_READ_COORDDATA, PSG_FINDRAYS, PSG_DEFOPTS, BTC_DEFINE,
%  PSG_RAYFIT, PSG_RAYANGLES, PSG_RAYMULTS, PSG_RAYSTATS, PSG_LLJIT_CRIT.
%
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays and psg_get_coordsets
if ~exist('opts_fit') opts_fit=struct(); end %for psg_rayfit
if ~exist('opts_stats') opts_stats=struct(); end %for psg_raystats
if ~exist('opts_lljit') opts_lljit=struct(); end %for psg_lljit_crit
if ~exist('opts_qpred') opts_qpred=struct(); end %for psg_qformpred
if ~exist('data_fullname') data_fullname=[]; end
if ~exist('setup_fullname') setup_fullname=[]; end
%
njit_types=2; %two types of critical jitter, see psj_lljit_crit
tag_coords='_coords_';
tag_choices='_choices_';
%
opts_lljit=filldefault(opts_lljit,'ndraws',100);
opts_lljit.ndraws=getinp('number of draws for critical jitter calculation','d',[1 10^4],opts_lljit.ndraws);
%
opts_stats=filldefault(opts_stats,'nsurrs',100);
opts_stats.nsurrs=getinp('number of surrogates for final confidence limits (0 to omit)','d',[0 10^4],opts_stats.nsurrs);
%
if ~exist('pval') pval=0.05; end
pval=getinp('p-value for confidence limits','f',[0 1],pval); 
%
jit_crit_choice=getinp('choice of critical jitter or 0 to combine','d',[0 njit_types],0);
%
opts_read=filldefault(opts_read,'if_log',1);
[sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred);
nsets=length(ds);
ray_labels=cell(nsets,1);
ray_pair_labels=cell(nsets,1);
mult_labels=cell(nsets,1);
%
angles=cell(nsets,2);
mults=cell(nsets,2);
angles_stats=cell(nsets,2);
mults_stats=cell(nsets,2);
rayfit=cell(nsets,2);
%
for iset=1:nsets
    nrays=rayss{iset}.nrays;
    disp(' ');
    disp(sprintf('set %1.0f: %s',iset,sets{iset}.label))
    nstims=size(ds{iset}{1},1);
    ndims=length(ds{iset});
    dim_list=sets{iset}.dim_list;
    model_dim_max=max(dim_list);
    %
    %read choice file, needed to find error bars
    %!!!! need to have option without error bars, which skips this, and
    %!!!! sets confidence lims and nsurrs to zero
    %
    data_fullname=opts_read_used{iset}.data_fullname;
    disp(sprintf('processing file: %s',strrep(data_fullname,'/','\')));
    if_havechoice=1;
    if contains(data_fullname,'_coords_')
        choices_fullname=strrep(data_fullname,tag_coords,tag_choices);
        disp(sprintf('   choice file:  %s',strrep(choices_fullname,'/','\')));
        if exist(choices_fullname,'file')
            c=load(choices_fullname);
            disp(sprintf('    nstims: %3.0f dims: %3.0f, cols in responses: %3.0f',nstims,ndims,size(c.responses,2)));
        else
            disp(sprintf('choice file not found: %s',strrep(choices_fullname,'/','\')));
            if_havechoice=0;
        end
    else
        if_havechoice=0;
    end
    %
    %compute critical jitters
    %
    jit_rms_list=zeros(model_dim_max,1);
    jit_crits=zeros(model_dim_max,njit_types); %type 1 and type 2 jits
    if if_havechoice
        for idimptr=1:length(dim_list) %compute ray fits and angles
            idim=dim_list(idimptr);
            jit_crits(idim,:)=psg_lljit_crit(pval,ds{iset}{idim},sas{iset}.typenames,c.responses,c.stim_list,opts_lljit);
        end       
        disp('critical jitters by type')
        for itype=1:njit_types
            jits_nan=find(isnan(jit_crits(:,itype)));
            if ~isempty(jits_nan) %if any of the jits are NaN, replace by maximum of non-nan jits
                jits_nonan=setdiff(dim_list,jits_nan);
                if ~isempty(jits_nonan)
                    jit_max=max(jit_crits(jits_nonan,itype));
                    jit_crits(jits_nan,itype)=jit_max;
                end
                disp(cat(2,' critical jitters were NaN for model dimensions:',sprintf(' %2.0f ',jits_nan))); ...
            end
            disp(sprintf(' %8.6f ',jit_crits(:,itype)))
        end
        if jit_crit_choice>0
            jit_rms_list=jit_crits(:,jit_crit_choice);
        else
            jit_rms_list=sqrt(sum(jit_crits.^2,2));
        end
        disp('rms jitter used')
        disp(sprintf(' %8.6f ',jit_rms_list));
    else
        disp('critical jitters not calculated, and confidence limits will be skipped');
    end
    %
    disp(sprintf('stimulus coordinates group along %2.0f rays',nrays));
    %
    ray_counts=full(sparse(rayss{iset}.whichray(rayss{iset}.whichray>0),ones(sum(rayss{iset}.whichray>0),1),1,nrays,1));
    for iray=1:nrays
        disp(sprintf('ray %2.0f: %2.0f points; endpoint: %s',iray,ray_counts(iray),sprintf('%5.2f',rayss{iset}.endpt(iray,:))));
    end
    %find stimulus label at end of ray, in positive direction when posssible
    ray_labels{iset}=cell(1,nrays);
    for iray=1:nrays
        mults_ray=rayss{iset}.mult(rayss{iset}.whichray==iray);
        maxend=intersect(find(abs(rayss{iset}.mult)==max(abs(mults_ray))),find(rayss{iset}.whichray==iray));
        maxend=maxend(find(rayss{iset}.mult(maxend)==max(rayss{iset}.mult(maxend)))); %choose positive direction if possible
        ray_labels{iset}{iray}=strrep(strrep(sas{iset}.spec_labels{maxend},' ',''),'=',''); %strip = and space
        disp(sprintf('ray %2.0f label: %s',iray,ray_labels{iset}{iray})); 
        mult_labels{iset}{iray}=...
            sprintf('%12s%5s    ',ray_labels{iset}{iray},'[-]',ray_labels{iset}{iray},'[+]',ray_labels{iset}{iray},'[-+]');
    end
    ray_pair_labels{iset}=cell(1,nrays*(nrays-1)/2);
    ilab=0;
    for iray=1:nrays-1
        for jray=iray+1:nrays
            ilab=ilab+1;
            ray_pair_labels{iset}{ilab}=cat(2,ray_labels{iset}{iray},':',ray_labels{iset}{jray});
        end
    end
    %
    %for each dimension model, find best-fitting signed and unsigned rays, including the origin
    %
    for idimptr=1:length(dim_list) %compute ray fits and angles
        idim=dim_list(idimptr);
        for if_bid=0:1 %uni- and bi-directional
            opts_stats_use=setfield(opts_stats,'if_bid',if_bid);
            if if_havechoice==0
                opts_stats_use=setfield(opts_stats_use,'nsurrs',0);
            end
            jit_use=jit_rms_list(idim);
            if isnan(jit_use) %will only happen if jitters are all NaN for all dimensions
                jit_use=0;
            end
            [angles{iset,1+if_bid}{idim},mults{iset,1+if_bid}{idim},angles_stats{iset,1+if_bid}{idim},mults_stats{iset,1+if_bid}{idim},rayfit{iset,1+if_bid}{idim}]=...
                psg_raystats(ds{iset}{idim},sas{iset},rayss{iset},jit_rms_list(idim),opts_stats_use);
        end %if_bid
    end %idim_ptr
    %
    %nicely formatted output to console
    %to be done: show gains, make a table
    %
    disp(' ');
    disp('gains along each ray, and for bidirectional fit to each axis')
    disp(cat(2,' dim    ',mult_labels{iset}{:}));
    for idimptr=1:length(dim_list) %display ray fits and angles
        idim=dim_list(idimptr);
        stat_names=fieldnames(angles_stats{iset,1}{idim});
        nstats=length(stat_names);
        for iv=0:nstats
            v=cell(1,2);
            if (iv==0)
                t=sprintf('%3.0f   ',idim);
                for ib=1:2
                    v{ib}=mults{iset,ib}{idim};
                end
            else
                stat_name=stat_names{iv};
                t=sprintf('%6s',stat_name);
                for ib=1:2
                    v{ib}=mults_stats{iset,ib}{idim}.(stat_name);
                end
            end
            for iray=1:nrays
                t=cat(2,t,sprintf(' %15.4f',v{1}.dist_gain(iray,1)),'     '); %gain on negative ray
                t=cat(2,t,sprintf(' %15.4f',v{1}.dist_gain(iray,2)),'     '); %gain on positive ray
                t=cat(2,t,sprintf(' %15.4f',v{2}.dist_gain(iray,1)),'     '); %bidirectional gain
            end
            disp(t);
        end %iv 0=data, >=1: stats
        disp(' ');
    end %idim_ptr
    disp(' ');
    disp('cosines of angles between pos and neg rays, and between bidirectional fits to axes')
    disp(cat(2,' dim    ',sprintf('%20s ',ray_labels{iset}{:}),sprintf('%30s ',ray_pair_labels{iset}{:}))); %header
    for idimptr=1:length(dim_list) %display ray fits and angles
        idim=dim_list(idimptr);
        stat_names=fieldnames(angles_stats{iset,1}{idim});
        nstats=length(stat_names);
        for iv=0:nstats
            v=cell(1,2);
            if (iv==0)
                t=sprintf('%3.0f   ',idim);
                for ib=1:2
                    v{ib}=angles{iset,ib}{idim};
                end
            else
                stat_name=stat_names{iv};
                t=sprintf('%6s',stat_name);
                for ib=1:2
                    v{ib}=angles_stats{iset,ib}{idim}.(stat_name);
                end
            end
            for iray=1:nrays
                t=cat(2,t,sprintf(' %15.4f',v{1}.cosangs(iray,iray,1,2)),'     '); %cosine of angle between pos and neg direction on each axis
            end
            for iray=1:nrays-1
                for jray=iray+1:nrays
                    t=cat(2,t,sprintf(' %25.4f',v{2}.cosangs(iray,jray)),'     '); %cosine of angle between bidirectional fits of two axes
                end
            end
            disp(t);
        end %iv 0=data, >=1: stats
        disp(' ');
    end %idim_ptr
    %create tables for mults and angles
    %All angles in same table but flag opposite angles at origin 
    %include p-value and number of surrogates and draws
    %one row fo each entry, include as metadata the set number, the subject, the config file, (eg bgca3pt), the dimension of the model, the paradigm
    % axes indicated by coordinate(s), up to 2, and length of axis on each coordinate, and with a letter label like p m z (bidirectional)
    % Remove labels  like bm0000
    % > 
    % >> Create tables  for gain bid angle and axis angles with set id subj 
    % >> id orig axis set Value Mean Conf Lim’s sem
    % > 
    % > Use same subj id convention as mtc tables Include dimension
    % > 
    % > Also make for qfm model (later)— can be customized or universal, and  designated in subj id
    % > 
    % mtc_mgm_ramp_tables.mat  
    % 
    % load mtc_mgm_ramp_tables
    % whos
    %   Name           Size               Bytes  Class    Attributes
    % 
    %   t_mdl      19860x16            19288745  table              
    %   t_psy       1498x16             1451657  table              
    % 
    % t_psy
    % t_psy =
    %   1498×16 table
    %     psy_model    subj_model_ID    expt_grp       cgroup1         cgroup2       expt_name       expt_uid          expt_type       plot_deg    thresh_mags_adj    thresh_mags_eblo_adj    thresh_mags_ebhi_adj    ray_angle     bexpon_mags    bexpon_mags_eblo    bexpon_mags_ebhi
    %     _________    _____________    _________    ____________    ____________    __________    _____________    _______________    ________    _______________    ____________________    ____________________    __________    ___________    ________________    ________________
    %      {'psy'}        {'jwb'}       {'expt2'}    {'AB_1_1'  }    {'AC_1_2'  }    {'YDM'   }    {'YDM'      }    {'mixed'      }        0             0.195                0.18                   0.209            3.5084e-15       2.375            2.276               2.462      
    %      {'psy'}        {'jwb'}       {'expt2'}    {'AB_1_1'  }    {'AC_1_2'  }    {'YDM'   }    {'YDM'      }    {'mixed'      }       30              0.22               0.204                   0.234                    30       2.375            2.276               2.462      
    %      {'psy'}        {'jwb'}       {'expt2'}    {'AB_1_1'  }    {'AC_1_2'  }    {'YDM'   }    {'YDM'      }    {'mixed'      }       60             0.277               0.262                   0.293                    60       2.375            2.276               2.462      
    %      {'psy'}        {'jwb'}       {'expt2'}    {'AB_1_1'  }    {'AC_1_2'  }    {'YDM'   }    {'YDM'      }    {'mixed'      }       90             0.325               0.307                   0.344                    90       2.375            2.276               2.462      
    %      {'psy'}        {'jwb'}       {'expt2'}    {'AB_1_1'  }    {'AC_1_2'  }    {'YDM'   }    {'YDM'      }    {'mixed'      }      120             0.309               0.291                   0.327                   120       2.375            2.276               2.462      
    %      {'psy'}        {'jwb'}       {'expt2'}    {'AB_1_1'  }    {'AC_1_2'  }    {'YDM'   }    {'YDM'      }    {'mixed'      }      150              0.26               0.245                   0.276                   150       2.375            2.276               2.462      
    %      {'psy'}        {'jwb'}       {'expt2'}    {'AB_1_1'  }    {'AC_1_2'  }    {'YDM'   }    {'YDM'      }    {'mixed'      }      180             0.211               0.196                   0.226                   180       2.375            2.276               2.462      
    %      {'psy'}        {'jwb'}       {'expt2'}    {'AB_1_1'  }    {'AC_1_2'  }    {'YDM'   }    {'YDM'      }    {'mixed'      }      210             0.214               0.201                   0.228                   210       2.375            2.276               2.462      
    % 

end %iset

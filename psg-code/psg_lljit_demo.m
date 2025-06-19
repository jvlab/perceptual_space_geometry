%psg_lljit_demo: significance testing for perceptual space geometry coords
%
% Note that this reads choice files directly, and bypasses psg_read_choicedata, so:
% *  the choice data files have tags in order of stim_list (while the
%    coordinate files are in order of typenames)
% *  the conversion from 'N_Repeats(D(ref, s1) > D(ref, s2))' to < is not carried out
% 
%  See also: PSG_GET_COORDSETS, PSG_LLJIT, PSG_LLJIT_CRIT.
%
if ~exist('opts_read') opts_read=struct();end %for psg_get_coordsets
if ~exist('opts_lljit') opts_lljit=struct;end
%
if ~exist('dim_sel') dim_sel=3; end %selected dimension to focus on 
if ~exist('pval_sel') pval_sel=0.05; end %p-value to focus on
%
if ~exist('jit_types') jit_types={'gaussian','shell'}; end
njtypes=length(jit_types);
%
opts_read.input_type=1;
opts_read.if_log=1;
opts_read.if_warn=0;
opts_read.if_data_only=1;
opts_read.nfiles_max=Inf;
%
opts_lljit=filldefault(opts_lljit,'if_frozen',1);
tag_coords='_coords_';
tag_choices='_choices_';
opts_lljit.if_frozen=...
    getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],opts_lljit.if_frozen);
%
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(opts_read);
nsets=length(sets);
%variables calculated by psg_lljit
lljit=cell(njtypes,nsets);
opts_lljit_used=cell(njtypes,nsets);
%variables calculated by psg_lljit_crit or locally
jit_crits=cell(njtypes,nsets);
lljit_crit=cell(njtypes,nsets);
opts_lljit_crit_used=cell(njtypes,nsets);
%
for iset=1:nsets
    disp(' ');
    nstims=size(ds{iset}{1},1);
    ndims=length(ds{iset});
    data_fullname=opts_read_used{iset}.data_fullname;
    disp(sprintf('processing file: %s',strrep(data_fullname,'/','\')));
    if_ok=1;
    if contains(data_fullname,'_coords_')
        choices_fullname=strrep(data_fullname,tag_coords,tag_choices);
        disp(sprintf('   choice file:  %s',strrep(choices_fullname,'/','\')));
        if exist(choices_fullname,'file')
            c=load(choices_fullname);
            disp(sprintf('    nstims: %3.0f dims: %3.0f, cols in responses: %3.0f',nstims,ndims,size(c.responses,2)));
        else
            disp('choice file not found');
            if_ok=0;
        end
    else
        if_ok=0;
    end
    if if_ok==1
        for ijt=1:njtypes
            opts_use=setfield(opts_lljit,'jit_type',jit_types{ijt});
            %coordinates are in the order of sas{iset}.typenames but choices are in the order of stim_list
            [lljit{ijt,iset},opts_lljit_used{ijt,iset}]=psg_lljit(ds{iset},sas{iset}.typenames,c.responses,c.stim_list,opts_use);
            disp(sprintf('psg_lljit done with jitter type %s',jit_types{ijt}));
            jit_crit_lljit=lljit{ijt,iset}.jit_crits(find(opts_lljit_used{ijt,iset}.pvals==pval_sel),dim_sel); %critical jitter type 1 from psg_lljit
            %
            [jit_crits{ijt,iset},lljit_crit{ijt,iset},opts_lljit_crit_used{ijt,iset}]=...
                psg_lljit_crit(pval_sel,ds{iset}{dim_sel},sas{iset}.typenames,c.responses,c.stim_list,opts_use);
            disp(sprintf('psg_lljit_crit done with jitter type %s',jit_types{ijt}));
            disp(sprintf('critical jitters for p-value of %9.6f on dim %3.0f',pval_sel,dim_sel));
            disp(sprintf('critical jitter 1 from psg_lljit: %9.6f, from psg_lljit_crit %9.6f diff %12.8f',...
                jit_crit_lljit,jit_crits{ijt,iset}(1),jit_crit_lljit-jit_crits{ijt,iset}(1)));
            disp(sprintf('critical jitters %9.6f %9.6f  sqrt(sum(squares)) %9.6f',jit_crits{ijt,iset},sqrt(sum(jit_crits{ijt,iset}.^2))));
        end
        if njtypes>1
            disp(sprintf('difference between two jitter types: %9.6f %9.6f',jit_crits{1,iset}-jit_crits{2,iset}));
        end
    end
end
%

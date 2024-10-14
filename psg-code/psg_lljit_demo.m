%psg_lljit_demo: significance testing for perceptual space geometry coords
%
%  See also: PSG_GET_COORDSETS, PSG_LLJIT, PSG_LLJIT_CRIT.
%
if ~exist('opts_read') opts_read=struct();end %for psg_get_coordsets
if ~exist('opts_lljit') opts_lljit=struct;end
%
if ~exist('dim_sel') dim_sel=3; end %selected dimension to focus on 
if ~exist('pval_sel') pval_sel=0.05; end %p-value to focus on
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
lljit=cell(1,nsets);
opts_lljit_used=cell(1,nsets);
%variables calculated by psg_lljit_crit or locally
jit_crit_type1=cell(1,nsets);
jit_crit_type2=cell(1,nsets);
lljit_crit=cell(1,nsets);
opts_lljit_crit_used=cell(1,nsets);
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
        %coordinates are in the order of sas{iset}.typenames but choices are in the order of stim_list
        [lljit{iset},opts_lljit_used{iset}]=psg_lljit(ds{iset},sas{iset}.typenames,c.responses,c.stim_list,opts_lljit);
        disp('psg_lljit done')
        jit_crit_lljit=lljit{iset}.jit_crits(find(opts_lljit_used{iset}.pvals==pval_sel),dim_sel); %critical jitter type 1 from psg_lljit
        %
        [jit_crits{iset},lljit_crit{iset},opts_lljit_crit_used{iset}]=...
            psg_lljit_crit(pval_sel,ds{iset}{dim_sel},sas{iset}.typenames,c.responses,c.stim_list,opts_lljit);
        disp('psg_lljit_crit done');
        disp(sprintf('critical jitters for p-value of %7.4f on dim %3.0f',pval_sel,dim_sel));
        disp(sprintf('critical jitter 1 from psg_lljit: %7.4f, from psg_lljit_crit %7.4f diff %12.8f',...
            jit_crit_lljit,jit_crits{iset}(1),jit_crit_lljit-jit_crits{iset}(1)));
        disp(sprintf('critical jitters %7.4f %7.4f  sqrt(sum(squares)) %7.4f',jit_crits{iset},sqrt(sum(jit_crits{iset}.^2))));
    end
end
%

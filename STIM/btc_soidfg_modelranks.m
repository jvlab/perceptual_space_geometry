%btc_soidfg_modelranks
%
%compute ranks of phenomenological models, and of projections of these models into the ideal-observer subspace
%
% See ..\gr18\figgnd_modeling_notes.docx.
% The projection into th ideal-observer subspace is given by the mapping Q to Q-M'QM, 
%  where M is [aI (1-a)I; aI (1-a)I];
%   
% See also:
% RANK, ORTH, BTC_SOIDFG_DEFINE, BTC_SOIDFG_DEMO, BTC_SOIDFG_MODEL, BTC_SOIFH_MODEL_PROJ.
%
btc_dict=btc_define;
btc_n=length(btc_dict.codel);
%
if ~exist('a_list') a_list=[0.5 0.25]; end %area fraction list
na=length(a_list);
Mlist=zeros(2*btc_n,2*btc_n,na);
for ia=1:na
    a=a_list(ia);
    Mlist(:,:,ia)=[a*eye(btc_n) (1-a)*eye(btc_n);a*eye(btc_n) (1-a)*eye(btc_n)]; %for projectiong with area fraction a
end
%
if ~exist('opts_fit') opts_fit=[]; end
opts_fit=btc_soidfg_define(opts_fit);
if ~exist('aug_opts') aug_opts=[]; end
%
% choose coordinates, model, and model spatial symmetry
%
rank_check_bad=0;
if ~exist('setups')
    coord_list_lib={'g','bc','gbc','bcde','gbcde','gbcdea','gbcdetuvw','gbcdetuvwa'};
    if ~exist('coord_list_ptrs') coord_list_ptrs=[1 2 3 4 6 8]; end
    for ilib=1:length(coord_list_lib)
        disp(sprintf('coordinate list %1.0f->%s',ilib,coord_list_lib{ilib}))
    end
    coord_list_ptrs=getinp('choice(s)','d',[1 length(coord_list_lib)],coord_list_ptrs);
    %
    model_type_lib=opts_fit.model_type_avail;
    if ~exist('model_type_ptrs') model_type_ptrs=[1:length(model_type_lib)];end
    for ilib=1:length(model_type_lib)
        disp(sprintf('model type %1.0f->%s',ilib,model_type_lib{ilib}));
    end
    model_type_ptrs=getinp('choice(s)','d',[1 length(model_type_lib)],model_type_ptrs);
    %
    sym_type_lib=opts_fit.sym_type_avail;
    if ~exist('sym_type_ptrs') sym_type_ptrs=[strmatch('none',sym_type_lib,'exact'),strmatch('hflip',sym_type_lib,'exact'),strmatch('full',sym_type_lib,'exact')]; end
    for ilib=1:length(sym_type_lib)
        disp(sprintf('spatial symmetry type %1.0f->%s',ilib,sym_type_lib{ilib}));
    end
    sym_type_ptrs=getinp('choice(s)','d',[1 length(sym_type_lib)],sym_type_ptrs);
end %setups exist?
%create setups and models
setups=cell(length(coord_list_ptrs),length(sym_type_ptrs),length(model_type_ptrs));
models=cell(length(coord_list_ptrs),length(sym_type_ptrs),length(model_type_ptrs));
%create setups and models
for coord_list_ptr=1:length(coord_list_ptrs)
    for sym_type_ptr=1:length(sym_type_ptrs)
        for model_type_ptr=1:length(model_type_ptrs)
            setups{coord_list_ptr,sym_type_ptr,model_type_ptr}.coords=coord_list_lib{coord_list_ptrs(coord_list_ptr)};
            setups{coord_list_ptr,sym_type_ptr,model_type_ptr}.model_type=model_type_lib{model_type_ptrs(model_type_ptr)};
            setups{coord_list_ptr,sym_type_ptr,model_type_ptr}.sym_type=sym_type_lib{sym_type_ptrs(sym_type_ptr)};
            models{coord_list_ptr,sym_type_ptr,model_type_ptr}=btc_soidfg_model(setups{coord_list_ptr,sym_type_ptr,model_type_ptr},btc_dict);
        end
    end
end %
%
nregressors=zeros(length(coord_list_ptrs),length(sym_type_ptrs),length(model_type_ptrs));
ranks=zeros(length(coord_list_ptrs),length(sym_type_ptrs),length(model_type_ptrs),na+1); %dim 4: raw rank, and ranks after projection by various values of area fraction
comb_ranks=zeros(length(coord_list_ptrs),length(sym_type_ptrs),length(model_type_ptrs),length(model_type_ptrs),na+1); %models combined in pairs: dim 5: raw rank, and ranks after projection
for coord_list_ptr=1:length(coord_list_ptrs)
    coord_list=coord_list_ptrs(coord_list_ptr);
    for sym_type_ptr=1:length(sym_type_ptrs)
        sym_type=sym_type_ptrs(sym_type_ptr);
        disp(' ')
        disp(sprintf(' coord list:  %10s, symmetry type: %15s',coord_list_lib{coord_list},sym_type_lib{sym_type}));
        qs=cell(1,length(model_type_ptrs));
        qas=cell(1,length(model_type_ptrs));
        for model_type_ptr=1:length(model_type_ptrs)
            model_type=model_type_ptrs(model_type_ptr);
            q=models{coord_list_ptr,sym_type_ptr,model_type_ptr}.qforms;
            nregressor=models{coord_list_ptr,sym_type_ptr,model_type_ptr}.nparams;
            %compute projections
            qa=zeros([2*btc_n,2*btc_n,nregressor,na]);
            for ir=1:size(q,3)
                for ia=1:na
                    qa(:,:,ir,ia)=q(:,:,ir)-Mlist(:,:,ia)'*q(:,:,ir)*Mlist(:,:,ia);
                end %ia
            end %ir
            q_reshaped=reshape(q,4*btc_n^2,nregressor);
            qa_reshaped=reshape(qa,4*btc_n^2,nregressor,na);
            qs{model_type_ptr}=q_reshaped;
            qas{model_type_ptr}=qa_reshaped;
            nregressors(coord_list_ptr,sym_type_ptr,model_type_ptr)=nregressor;
            rank_f=rank(q_reshaped);
            rank_a=zeros(1,na);
            rank_a_proj=zeros(1,na);
            model_proj=cell(1,na);
            model_sort=cell(1,na);
            r_proj=cell(1,na);
            for ia=1:na
                rank_a(ia)=rank(qa_reshaped(:,:,ia));
                %check rank against rank of model projection
                [model_proj{ia},model_sort{ia},M_proj,r_proj{ia}]=btc_soidfg_model_proj(models{coord_list_ptr,sym_type_ptr,model_type_ptr},a_list(ia),btc_dict);
                rank_a_proj(ia)=model_proj{ia}.nparams;
            end
            ranks(coord_list_ptr,sym_type_ptr,model_type_ptr,:)=reshape([rank_f,rank_a],[1 1 1 1+na]);
            disp(sprintf('model type  %20s:  regressors: %5.0f, full phenom model rank %5.0f (%3.0f loss)',...
                model_type_lib{model_type},nregressor,rank_f,nregressor-rank_f));
            for ia=1:na
                disp(sprintf('   proj to IO model with area frac  %5.2f: IO model rank %3.0f (%3.0f loss)',a_list(ia),rank_a(ia),rank_f-rank_a(ia)));
                if rank_a_proj(ia)~=rank_a(ia)
                    disp(sprintf(' inconsistency: projected rank from btc_soidfg_model_proj is %4.0f',rank_a_proj(ia)));
                    disp(model_proj{ia})
                    rank_check_bad=rank_check_bad+1;
                end
            end
            %analyze how models combine in pairs
            m1_ptr=model_type_ptr;
            m1_type=model_type;
            q1=models{coord_list_ptr,sym_type_ptr,m1_ptr}.qforms;
            for m2=1:m1_ptr-1
                m2_ptr=m2;
                m2_type=model_type_ptrs(m2_ptr);
                q2=models{coord_list_ptr,sym_type_ptr,m2_ptr}.qforms;
                %combine the parameters
                q_comb=cat(3,q1,q2);
                names_comb=[models{coord_list_ptr,sym_type_ptr,m1_ptr}.param_name,models{coord_list_ptr,sym_type_ptr,m2_ptr}.param_name];
                [name_unique,ptrs_unique]=unique(names_comb);
                q_unique=zeros(2*btc_n,2*btc_n,length(ptrs_unique),1+na); %dim 4=1 for unique regressors, dim 4> 1 for projections into IO models
                q_unique(:,:,:,1)=q_comb(:,:,ptrs_unique);
                %project into IO models
                for ia=1:na
                    for ir=1:length(ptrs_unique)
                        q_unique(:,:,ir,1+ia)=q_unique(:,:,ir,1)-Mlist(:,:,ia)'*q_unique(:,:,ir,1)*Mlist(:,:,ia);
                    end
                end
                %q_orth=orth(reshape(q_unique(:,:,:,1),4*btc_n^2,size(q_unique,3)));
                %q_rank=size(q_orth,2)
                q_rank=rank(reshape(q_unique(:,:,:,1),4*btc_n^2,size(q_unique,3)));
                disp(sprintf('model type %20s combined with %20s: params %3.0f + %3.0f, %3.0f unique names (%3.0f loss), span size %3.0f (%3.0f loss)',...
                	model_type_lib{m1_type},model_type_lib{m2_type},size(q1,3),size(q2,3),size(q_unique,3),size(q_comb,3)-size(q_unique,3),...
                    q_rank,size(q_comb,3)-q_rank));
                comb_ranks(coord_list_ptr,sym_type_ptr,m1_ptr,m2_ptr,1)=q_rank;
                comb_ranks(coord_list_ptr,sym_type_ptr,m2_ptr,m1_ptr,1)=q_rank;
                for ia=1:na
                    %q_orth=orth(reshape(q_unique(:,:,:,1+ia),4*btc_n^2,size(q_unique,3)));
                    %q_rank=size(q_orth,2)
                    q_rank=rank(reshape(q_unique(:,:,:,1+ia),4*btc_n^2,size(q_unique,3)));
                    disp(sprintf(' proj to IO model with area frac %5.2f: model rank %3.0f (%3.0f loss)',...
                        a_list(ia),q_rank,comb_ranks(coord_list_ptr,sym_type_ptr,m1_ptr,m2_ptr,1)-q_rank));
                    comb_ranks(coord_list_ptr,sym_type_ptr,m1_ptr,m2_ptr,1+ia)=q_rank;
                    comb_ranks(coord_list_ptr,sym_type_ptr,m2_ptr,m1_ptr,1+ia)=q_rank;
                end %ia
            end %m_ptr3
        end %model_type_ptr
    end %coord_list_ptr
end %sym_type_ptr
disp(' ');
disp('summary, single models');
for model_type_ptr=1:length(model_type_ptrs)
    model_type=model_type_ptrs(model_type_ptr);
    disp(sprintf('model type  %20s:',model_type_lib{model_type}));
    disp(sprintf('maximum loss of rank for individual models (regressors - full rank)          : %3.0f',max(max(max(nregressors(:,:,model_type_ptr)-ranks(:,:,model_type_ptr,1))))));
    for ia=1:na
        disp(sprintf('maximum loss of rank between phenomenological model and IO model with a=%5.2f: %3.0f',a_list(ia),max(max(max(ranks(:,:,model_type_ptr,1)-ranks(:,:,model_type_ptr,ia+1))))));
    end
end
disp(' ');
disp('summary, model pairs');
for m1_ptr=1:length(model_type_ptrs)
    m1_type=model_type_ptrs(m1_ptr);
    for m2_ptr=1:m1_ptr-1
        m2_type=model_type_ptrs(m2_ptr);
        disp(sprintf('model types  %20s and %20s:',model_type_lib{m1_type},model_type_lib{m2_type}));
        for ia=1:na
            disp(sprintf('maximum loss of rank between phenomenological model and IO model with a=%5.2f: %3.0f',a_list(ia),...
                max(max(max(comb_ranks(:,:,m1_ptr,m2_ptr,1)-comb_ranks(:,:,m1_ptr,m2_ptr,ia+1))))));
        end
    end
end
disp(sprintf('comparison with btc_soidfg_model_proj rank: %4.0f inconsistencies.',rank_check_bad));
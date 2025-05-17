function [model_proj,model_sort,M,r]=btc_soidfg_model_proj(model,a,btc_dict)
% [model_proj,model_sort,M,r]=btc_soidfg_model_proj(model,a,btc_dict) projects a model into the ideal-observer subspace
% defined by Q->Q-M'QM, where M is a % See ..\gr18\figgnd_modeling_notes.docx.
%
% The projection into th ideal-observer subspace is given by the mapping Q to Q-M'QM, 
%  where M is [aI (1-a)I; aI (1-a)I];
%
% input:
%  model: figure-ground model structure, typcally created by btc_soidfg_model
%  a: a number in [0 1], corresponding to the fraction of the stimulus area that is the figure
%  btc_dict: structure returned by btc_define; can be omitted
%
% output:
%  model_proj: structure defining the model after projection.  fields listed below may be changed, as well as nparams.
%    This has a field "area_frac", set equal to a,
%    params that have been (non-trivially) modified by projection have a * added to their names
%    params that do not add to the rank are deleted
%  model_sort: model after sorting the parameters by parameter type
%    the following fields may be reordered but otherwise unchanged:
%     qforms  
%     param_name
%     param_type
%     param_name_equiv
%  M: the above matrix M
%  r: a structure of intermediate results
%     r.warn: warning messages (unexpected changes in a projection or loss of rank)
%     r.sort_order:  origin of the parameters in model_sort: model_sort.qforms=model.qforms(:,:,r.sort_order)
%     r.params_orig.(type): pointers to each type of parameter in original model
%     r.proj_order:  origin of the parameters in model_proj, model_proj.qforms=projection of model.qforms(:,:,r.proj_order).  
%        original qforms slices that do not increase rank are omitted
%     r.proj_change:  amount that the projection changed, for projections included in model_proj
%
% model parameters are sorted so that they can be added into model_proj in
% a logical fashion, giving priority to the fmg (figure - ground) terms (which should be un-altered by the projection)
%
%  See also:  BTC_SOIDFG_DEFINE, BTC_DEFINE, BTC_SOIDFG_MODEL_TEST, BTC_SOIDFG_MODEL, BTC_SOIDFG_MODELRANKS, RANK, ORTH.
%
type_sorted={'self_fmg','cross_fmg','self_fig','cross_fig','self_gnd','cross_gnd','self_fig_plus_self_gnd','cross_fig_plus_cross_gnd'}; %,then [any other self],[any other cross],[any other]}
%every entry of type_sorted must be in type_unchanged or type_keeprank
type_unchanged={'self_fmg','cross_fmg'}; %parameter types that should not change under projection
type_keeprank={'self_fig','cross_fig','self_gnd','cross_gnd','self_fig_plus_self_gnd','cross_fig_plus_cross_gnd'}; %these parameters change under projection but should maintain rank
type_misc={}; %other types encountered
proj_tol=10^-5; %tolerance for checking that a projection has not changed something
%
if (nargin<=2)
    btc_dict=btc_define([]);
end
codel=btc_dict.codel;
btc_n=length(codel); %10
%
r=[];
m=rmfield(model,{'qforms','param_name','param_type','param_name_equiv'});
model_proj=m;
model_sort=m;
r=[];
M=[a*eye(btc_n) (1-a)*eye(btc_n);a*eye(btc_n) (1-a)*eye(btc_n)]; %for projectiong with area fraction a
%
% sort order according to parameter type
%
params_used=zeros(1,m.nparams);
type_index=0;
sort_order=[];
while type_index<length(type_sorted) & any(params_used==0)
    type_index=type_index+1;
    param_indices=strmatch(type_sorted{type_index},model.param_type,'exact');
    if ~isempty(param_indices)
        sort_order=[sort_order;param_indices];
        params_used(param_indices)=1;
        r.params_orig.(type_sorted{type_index})=param_indices(:)';
    end
end
%self params not included in the above
param_indices=intersect(find(params_used==0),strmatch('self',model.param_type));
if ~isempty(param_indices)
    type_misc{1,end+1}='self_misc';
    sort_order=[sort_order;param_indices];
    params_used(param_indices)=1;
    r.params_orig.(type_misc{end})=param_indices(:)';
end
%cross params not included in the above
param_indices=intersect(find(params_used==0),strmatch('cross',model.param_type));
if ~isempty(param_indices)
    type_misc{1,end+1}='cross_misc';
    sort_order=[sort_order;param_indices];
    params_used(param_indices)=1;
    r.params_orig.(type_misc{end})=param_indices(:)';
end
%any remaining params
param_indices=find(params_used==0);
if ~isempty(param_indices)
    type_misc{1,end+1}='misc';
    sort_order=[sort_order;param_indices];
    r.params_orig.(type_misc{end})=param_indices(:)';
end
%
model_sort.qforms=model.qforms(:,:,sort_order);
model_sort.param_name=model.param_name(sort_order);
model_sort.param_type=model.param_type(sort_order);
model_sort.param_name_equiv=model.param_name_equiv(sort_order);
r.sort_order=sort_order(:)';
%
model_proj.area_frac=a;
%first install the fmg parameters, these should not change under projection
new_ind=0;
r.warn=[];
r.proj_order=[];
model_proj.qforms=[];
for icat=1:3 %3 categories: expect no change, expect change but increase in rank, and misc
    switch icat
        case 1
            type_cat=type_unchanged;
            if_expect_change=0;
        case 2
            type_cat=type_keeprank;
            if_expect_change=1;
        case 3
            type_cat=type_misc;
            if_expect_change=1;
    end
    for itype=1:length(type_cat)
        param_type=type_cat{itype};
        if isfield(r.params_orig,param_type)
            for ip=1:length(r.params_orig.(param_type))
                old_ind=r.params_orig.(param_type)(ip);
                q=model.qforms(:,:,old_ind);
                q_proj=q-M'*q*M; %project it
                q_chg=max(abs(q(:)-q_proj(:))); %change under projection
                if (q_chg>=proj_tol)
                    proj_mark='*';
                    if if_expect_change==0
                        warn_msg=sprintf('parameter %25s changes by %12.5f under projection, expect no change',model.param_name{old_ind},q_chg);
                        r.warn=strvcat(r.warn,warn_msg);
                    end
                else
                    proj_mark='';
                end
                %check rank
                q_rank=cat(2,reshape(model_proj.qforms,[4*btc_n^2,new_ind]),q_proj(:));
                if rank(q_rank)<=new_ind
                    %rank has not increased
                    warn_msg=sprintf('parameter %25s does not increase the rank, and is not added',model.param_name{old_ind});
                    r.warn=strvcat(r.warn,warn_msg);
                else
                    new_ind=new_ind+1;
                    model_proj.qforms(:,:,new_ind)=q_proj;
                    model_proj.param_name{new_ind}=cat(2,model.param_name{old_ind},proj_mark);
                    model_proj.param_type{new_ind}=model.param_type{old_ind};
                    neq=length(model.param_name_equiv{old_ind});
                    for ieq=1:neq
                        model_proj.param_name_equiv{new_ind}{ieq}=cat(2,model.param_name_equiv{old_ind}{ieq},proj_mark);
                    end
                    r.proj_change(1,new_ind)=q_chg; 
                    r.proj_order(1,new_ind)=old_ind;
                end
            end
        end %isfield
    end %itype
end %icat
%install the parameters that are expected to change after projection
model_proj.nparams=new_ind;
return

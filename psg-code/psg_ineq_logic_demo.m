%psg_ineq_logic_demo: demonstrate psg_ineq_logic, psg_ineq_edgecount
%
%   See also: PSG_INEQ_LOGIC, PSG_TENTLIKE_DEMO, PSG_UMI_TRIPLIKE_DEMO, PSG_INEQ_EDGECOUNT.
%
% 31Aug23: minor formatting changes in output
%
types=psg_ineq_logic();
type_names=fieldnames(types);
p=struct;
for itype=1:length(type_names)
    type_name=type_names{itype};
    nc_list=types.(type_name).nc;
    if isempty(nc_list)
        nc_list=[1 3 6];
    end
    for nc_ptr=1:length(nc_list)
        nc=nc_list(nc_ptr);
        disp(sprintf(' inequality logic for type %32s, nc=%3.0f',type_name,nc));
        p.(type_name){nc}=psg_ineq_logic(nc,type_name,setfield([],'if_log',1));
    end
end
%test for symmetries
disp('testing flip symmetry (changing sign of comparison)');
for itype=1:length(type_names)
    type_name=type_names{itype};
    nc_list=types.(type_name).nc;
    if isempty(nc_list)
        nc_list=[1 3 6];
    end
    for nc_ptr=1:length(nc_list)
        nc=nc_list(nc_ptr);
        p_orig=p.(type_name){nc};
        p_flipped=p_orig;
        for id=1:nc
            p_flipped=flip(p_flipped,id);
        end
        ndiff=sum(p_flipped(:)~=p_orig(:));
        if (ndiff==0) & types.(type_name).flip_invar==1
            okstring='OK, expected invariant';
        elseif (ndiff>0) & types.(type_name).flip_invar==0
            okstring='OK, expected change';
        else
            okstring='not OK';
        end
        disp(sprintf(' type %32s, nc=%3.0f number changed: %7.0f (%6s)',type_name,nc,ndiff,okstring));
    end
end
%test for symmetries
disp('testing cycle symmetry (cycling the stimuli)');
for itype=1:length(type_names)
    type_name=type_names{itype};
    nc_list=types.(type_name).nc;
    if isempty(nc_list)
        nc_list=[1 3 6];
    end
    for nc_ptr=1:length(nc_list)
        nc=nc_list(nc_ptr);
        p_orig=p.(type_name){nc};
        if (nc>1)
            p_cycled=permute(p_orig,[2:nc 1]);
        else
            p_cycled=p_orig;
        end
        ndiff=sum(p_cycled(:)~=p_orig(:));
        if (ndiff==0) & types.(type_name).cycle_invar==1
            okstring='OK, expected invariant';
        elseif (ndiff>0) & types.(type_name).cycle_invar==0
            okstring='OK, expected change';
        else
            okstring='not OK';
        end
        disp(sprintf(' type %32s, nc=%3.0f number changed: %7.0f (%6s)',type_name,nc,ndiff,okstring));
    end
end

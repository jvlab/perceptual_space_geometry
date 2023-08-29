%psg_logic_v_demo: demonstrate psg_ineq_logic, with attention to "V" notation
% in ori_char_sym manuscript
%
%   See also: PSG_INEQ_LOGIC, PSG_LOGIC_DEMO.
%
V=struct;
V_edge_counts=struct;
V_binary=struct;
V_all=struct;
%
%table of correspondences
%
corresp.sym.type='exclude_sym';
corresp.sym.nc=3;
corresp.umi.type='exclude_umi_trans';
corresp.umi.nc=3;
corresp.symtent.type='exclude_trans_tent';
corresp.symtent.nc=6;
corresp.addtree.type='exclude_addtree_trans';
corresp.addtree.nc=6;
%
inclusions={{'sym','umi'},{'umi','sym'},{'symtent','addtree'},{'addtree','symtent'}};
%
vtypes=fieldnames(corresp);
for iv=1:length(vtypes)
    vtype=vtypes{iv};
    nc=corresp.(vtype).nc;
    combs=2.^(nc-[0:nc])*factorial(nc)./factorial([0:nc])./factorial(nc-[0:nc]); %total number of regions with a given number of edges
    V_all.(vtype)=combs;
    [parts,ou]=psg_ineq_logic(corresp.(vtype).nc,corresp.(vtype).type);
    V.(vtype)=1-parts;
    V_edge_counts.(vtype)=combs-ou.edge_counts;
    V_binary.(vtype)=1-ou.partitions_nz;
end
for iv=1:length(vtypes)
    vtype=vtypes{iv};
    nc=corresp.(vtype).nc;
    disp(sprintf('V %-10s, taken from 1 - %s ',vtype,corresp.(vtype).type));
    for id=0:nc
        disp(sprintf(' %5.0f included regions with %4.0f ties (out of %5.0f)',...
            V_edge_counts.(vtype)(id+1),id,V_all.(vtype)(id+1)));
    end
end
% check inclusion
for iinc=1:length(inclusions)
    disp(sprintf('testing that V%-10s includes V%-10s',inclusions{iinc}{1},inclusions{iinc}{2}));
    Vnum=V.(inclusions{iinc}{1});
    Vden=V.(inclusions{iinc}{2});
    disp(sprintf('full:    %4.0f exceptions',sum(Vden(:).*(Vden(:)-Vnum(:)))));
    Vnum=V_binary.(inclusions{iinc}{1});
    Vden=V_binary.(inclusions{iinc}{2});
    disp(sprintf('binary:  %4.0f exceptions',sum(Vden(:).*(Vden(:)-Vnum(:)))));

end

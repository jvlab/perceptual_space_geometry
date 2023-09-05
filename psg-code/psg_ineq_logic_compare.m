%psg_ineq_logic_compare
% compare (and test) a priori log likelihoods and inequality structures based on tetrahedra,
% looking at addtree with transitivity/symmetry conditions derived
% from tents (i.e., 3-cycles), vs. 4-cycles, vs. 3-cycles and 4-cycles
% 
% also compare with a priori ratios structures based on tents for addtree and sym
% also compare with a priori ratios based on triplets for sym
%
% See also  PSG_INEQ_LOGIC, PSG_INEQ_LOGIC_TETRA.
%
ptetra=struct;
qtetra=struct;
ptent=struct;
qtent=struct;
ptriplet=struct;
qtriplet=struct;
%
apriori_adt=struct;
apriori_sym=struct;
%
nc_tetra=12;
nc_tent=6;
nc_triplet=3;
%
qtot_tetra=2^nc_tetra;
qtot_tent=2^nc_tent;
qtot_triplet=2^nc_triplet;
%
%calculations for tetrahedra
%
[partitions,opts_used]=psg_ineq_logic(nc_tetra,'exclude_trans_tetra');ptetra.sym.tent=partitions;qtetra.sym.tent=opts_used.partitions_nz; %tetrahedra excluded based on tents
[partitions,opts_used]=psg_ineq_logic(nc_tetra,'exclude_trans_tetra_cyc4');ptetra.sym.cyc4=partitions;qtetra.sym.cyc4=opts_used.partitions_nz; %tetrahedra excluded based on 4-cycles
[partitions,opts_used]=psg_ineq_logic(nc_tetra,'exclude_trans_tetra_all');ptetra.sym.all=partitions;qtetra.sym.all=opts_used.partitions_nz; %tetrahedra excluded based on tents and 4-cycles
%
[partitions,opts_used]=psg_ineq_logic(nc_tetra,'exclude_addtree_trans_tetra');ptetra.adt.tent=partitions;qtetra.adt.tent=opts_used.partitions_nz; %exclude tetras on addtree and trans for tents
%not a standard optoin in psg_ineq_logic
[p1,opts_used]=psg_ineq_logic(nc_tetra,'exclude_addtree_trans_tetra');q1=opts_used.partitions_nz; %exclude tetras on addtree
ptetra.adt.cyc4=min(1,p1+ptetra.sym.cyc4);qtetra.adt.cyc4=min(1,q1+qtetra.sym.cyc4); %exclude tetras based on addtree and trans for 4-cycle
[partitions,opts_used]=psg_ineq_logic(nc_tetra,'exclude_addtree_trans_tetra_all');ptetra.adt.all=partitions;qtetra.adt.all=opts_used.partitions_nz; %exclude tetras on addtree and trans for tents and 4-cycles
%
%calculations for tent
%
[partitions,opts_used]=psg_ineq_logic(nc_tent,'exclude_trans_tent');ptent.sym.tent=partitions;qtent.sym.tent=opts_used.partitions_nz; %excluded based on tents
[partitions,opts_used]=psg_ineq_logic(nc_tent,'exclude_addtree_trans');ptent.adt.tent=partitions;qtent.adt.tent=opts_used.partitions_nz; %excluded based on tents
%
%calculations for triplets
%
[partitions,opts_used]=psg_ineq_logic(nc_triplet,'exclude_trans');ptriplet.sym.triplet=partitions;qtriplet.sym.triplet=opts_used.partitions_nz; %excluded based on triplets
%
condits=fieldnames(ptetra);
for icondit=1:length(condits)
    condit=condits{icondit};
    disp(sprintf('for condition %s, considering tetrahedra',condit))
    disp(sprintf('        masked via tents: %7.0f  via cyc4: %7.0f   via all  %7.0f   via tents but not cyc4 %7.0f   via cyc4 but not tents %7.0f  via both %7.0f',...
        sum(qtetra.(condit).tent(:)),sum(qtetra.(condit).cyc4(:)),sum(qtetra.(condit).all(:)),...
        sum(max(0,qtetra.(condit).tent(:)-qtetra.(condit).cyc4(:))),...
        sum(max(0,qtetra.(condit).cyc4(:)-qtetra.(condit).tent(:))),...
        sum(qtetra.(condit).cyc4(:)&qtetra.(condit).tent(:))));
    disp(sprintf('      unmasked via tents: %7.0f  via cyc4: %7.0f   via all  %7.0f',...
        qtot_tetra-sum(qtetra.(condit).tent(:)),qtot_tetra-sum(qtetra.(condit).cyc4(:)),qtot_tetra-sum(qtetra.(condit).all(:))));
    disp(sprintf('for condition %s, considering tents',condit)) 
    disp(sprintf('        masked via tents: %7.0f',sum(qtent.(condit).tent(:))));
    disp(sprintf('      unmasked via tents: %7.0f',qtot_tent-sum(qtent.(condit).tent(:))));
    disp(' ');
end
%compare addtree with sym, vs sym alone
apriori_adt.tetra_tent=(qtot_tetra-sum(qtetra.adt.tent(:)))/(qtot_tetra-sum(qtetra.sym.tent(:)));
apriori_adt.tetra_cyc4=(qtot_tetra-sum(qtetra.adt.cyc4(:)))/(qtot_tetra-sum(qtetra.sym.cyc4(:)));
apriori_adt.tetra_all=(qtot_tetra-sum(qtetra.adt.all(:)))/(qtot_tetra-sum(qtetra.sym.all(:)));
apriori_adt.tent=(qtot_tent-sum(qtent.adt.tent(:)))/(qtot_tent-sum(qtent.sym.tent(:)));
%
disp(' for addtree and symmetry/transitivity vs. symmetry/transitivity alone')
disp(sprintf('apriori ratios via tents: %7.4f  via cyc4: %7.4f   via all  %7.4f; considering tents only: %7.4f',...
    apriori_adt.tetra_tent,apriori_adt.tetra_cyc4,apriori_adt.tetra_all,apriori_adt.tent));
disp(sprintf('apriori llrs   via tents: %7.4f  via cyc4: %7.4f   via all  %7.4f; considering tents only: %7.4f',...
    log(apriori_adt.tetra_tent),log(apriori_adt.tetra_cyc4),log(apriori_adt.tetra_all),log(apriori_adt.tent)));
disp(' ')
%compare sym, vs no constraints
apriori_sym.tetra_tent=(qtot_tetra-sum(qtetra.sym.tent(:)))/qtot_tetra;
apriori_sym.tetra_cyc4=(qtot_tetra-sum(qtetra.sym.cyc4(:)))/qtot_tetra;
apriori_sym.tetra_all=(qtot_tetra-sum(qtetra.sym.all(:)))/qtot_tetra;
apriori_sym.tent=(qtot_tent-sum(qtent.sym.tent(:)))/qtot_tent;
apriori_sym.triplet=(qtot_triplet-sum(qtriplet.sym.triplet(:)))/qtot_triplet;
disp(' for symmetry/transitivity vs. no constraints')
disp(sprintf('apriori ratios via tents: %7.4f  via cyc4: %7.4f   via all  %7.4f; considering tents only: %7.4f; considering triplets only: %7.4f',...
    apriori_sym.tetra_tent,apriori_sym.tetra_cyc4,apriori_sym.tetra_all,apriori_sym.tent,apriori_sym.triplet));
disp(sprintf('apriori llrs   via tents: %7.4f  via cyc4: %7.4f   via all  %7.4f; considering tents only: %7.4f; considering triplets only: %7.4f',...
    log(apriori_sym.tetra_tent),log(apriori_sym.tetra_cyc4),log(apriori_sym.tetra_all),log(apriori_sym.tent),log(apriori_sym.triplet)));


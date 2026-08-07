%rigidity_constraints_demo: demonstrate rigidity constraints for an odorant panel designnpts
%
% novlp2 refers to overlaps between each panel and the next, so for the
% 2-panel case, there are overlaps on both sides if cyclic, only on one side if not
%
%  See also: RIGIDITY_CONSTRAINTS.
%
if ~exist('nset_list') nset_list=[2:6];end  %number of panels; no elements less than 2
if ~exist('npanel_list') npanel_list=17; end %number in each panel 
if ~exist('novlp2_list') novlp2_list=[0:7]; end %number in pairwise overlaps
if ~exist('ncommon_list') ncommon_list=[0:5];end %number in common to all in addition to pairwise
if ~exist('ifcyc_list') ifcyc_list=[0 1]; end % if pairwise overlaps are cyclic
%
results.colnames={...
    'tot_odorants_per_panel',...
    'odorants_in_2_adjacent_panels',...
    'additional_odorants_common_to_all_panels',...
    'number_of_panels',...
    'if_cyclic_overlaps',...
    'distance_count_within_panels',...
    'number_of_distances_with_repeated_measures',...
    'total_distances_if_all_measured_in_one_set',...
    'max_allowed_dimension_based_on_number_of_distances_measured',...
    'max_allowed_dimension_based_on_points_not_duplicated_between_sets',...
    'max_allowed_dimension_overall'};
%
disp('  tot odorants  odorants in   addl odorants    number       if       total         dists       repeated   total dists    max dim by   max dim by    max dim');
disp('   per panel   2 adj panels   in all panels  of panels    cyclic    odorants   within panels   measures   if combined   constraints    free pts     overall ');
%
row_ptr=0;
for ipanel=1:length(npanel_list)
    npanel=npanel_list(ipanel);
    for iovlp2=1:length(novlp2_list)
        novlp2=novlp2_list(iovlp2);
        for icommon=1:length(ncommon_list)
            ncommon=ncommon_list(icommon);
            for iset=1:length(nset_list)
                nsets=nset_list(iset);
                for ifcyc_ptr=1:length(ifcyc_list)
                    ifcyc=ifcyc_list(ifcyc_ptr);
                    nstep=npanel-ncommon-novlp2; %increment in first odor index with each panel
                    if (nstep>=0)
                        ncyc=(npanel-ncommon)*nsets-novlp2*(nsets+ifcyc-1);
                        ntot=ncyc+ncommon;
                        ovlps=zeros(ntot,nsets);
                        ovlps(ncyc+[1:ncommon],:)=1;
                        for kset=1:nsets
                            ovlps(1+mod((kset-1)*nstep+[0:npanel-1-ncommon],ncyc),kset)=1;
                        end
            %            disp(sprintf(' odorants per panel: %3.0f, odors in 2 panels: %3.0f, odors in all panels: %3.0f, panels: %2.0f',npanel,novlp2,ncommon,nsets))
                        [params,opts_used]=rigidity_constraints(ovlps,struct());
                        npairs=params.npairs;
                        nmax=params.nmax;
                        counts=params.counts;
                        nrepeated=sum(counts(:))/2-npairs;
                        row_ptr=row_ptr+1;
                        results.data(row_ptr,:)=[npanel novlp2 ncommon nsets ifcyc ntot npairs nrepeated nmax params.dmax_constraints params.dmax_free params.dmax];
                        disp(sprintf('%10.0f   ',results.data(row_ptr,:)));
                    end %nsole>0
                end %ifcyc
            end %iset
        end %icommon
    end %iovlp2
end %ipanel

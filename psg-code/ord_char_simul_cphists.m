%ord_char_simul_cphists: show choice probability histograms from rank order anahysis of toy models
%
% run after load results from the output workspace of ord_char_simul_demo.
%
% See also: ORD_CHAR_SIMUL_DEMO, ORD_CHAR_SIMUL_PLOT.
%
if ~exist('if_poisson')
    if_poisson=0;
end
%show underlying and empiric choice probabilities
for itrial_ptr=0:ntrials_list
    figure; %choice probabilities
    set(gcf,'Position',[100 100 1200 800]);
    set(gcf,'Numbertitle','off');
    if (itrial_ptr>0)
        name_string=sprintf('empiric choice probabilities, %4.2f trials, if_poisson=%1.0f',trials_list(itrial_ptr),if_poisson);
    else
        name_string='choice probabilities';
    end
    set(gcf,'Name',name_string);
    for irule=1:nrules
        button=rules(irule,1);
        sigma=rules(irule,2);
        for itype=1:ndist_types
            dist_type=dist_types{itype};
            subplot(ndist_types,nrules,irule+(itype-1)*nrules);
            if (itrial_ptr>0)
                ntrials=trials_list(itrial_ptr);
                cp_list=results{itype,irule}.C{itrial_ptr}/ntrials;
                bincenters_cp=[0:ntrials]/ntrials;
            else
                cp_list=results{itype,irule}.cp_table(:,4);
                bincenters_cp=[0.5:nbins_cp-0.5]*(1/(nbins_cp-1));
            end
            hist([cp_list;1-cp_list],bincenters_cp);
            set(gca,'XLim',[0 1]);
            set(gca,'XTick',[0 0.5 1]);
            if (itype==1)
                title(results{itype,irule}.rule_string_nice)
            end
            if (irule==1)
                ylabel(results{itype,irule}.dist_type,'Interpreter','none');
            end
        end %itype
    end %irule
    axes('Position',[0.01,0.02,0.01,0.01]); %for text
    text(0,0,cat(2,name_string,'; ',lab_string),'Interpreter','none','FontSize',10);
    axis off;
end
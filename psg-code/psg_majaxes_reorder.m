function [zplot,xtick_labels,match_list]=psg_majaxes_reorder(z,typenames_plot,typenames_std)
% [zplot,xtick_labels,match_list]=psg_majaxes_reorder(z,typenames_plot,typenames_std)
% selects and reorder a matrix for plotting in psg_majaxes
%
%   See also: PSG_MAJAXES.
% 
nstims_plot=length(typenames_plot);
zplot=zeros(nstims_plot,size(z,2));
xtick_labels=cell(nstims_plot,1);
match_list=zeros(nstims_plot,1);
for istim=1:nstims_plot
    imatch=min(strmatch(typenames_plot{istim},typenames_std,'exact'));
    if ~isempty(imatch)
        match_list(istim)=imatch;
        zplot(istim,:)=z(imatch,:);
        xtick_labels{istim}=typenames_std{imatch};
    else
        xtick_labels{istim}=cat(2,'[',typenames_plot{istim},']');
    end
end
return

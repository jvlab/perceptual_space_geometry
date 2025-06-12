%psg_raystats_dbplot_axes: module to adjust axes in psg_raystats_dbplot
%
%key variables:
%ylims(irow,icol,ipage,:)=get(gca,'YLim');
%subhandles{irow,icol,ipage}=gca;
%
% See also: PSG_RAYSTATS_DBPLOT_AXES.
%
axis_opt=-1;
axis_opts={'revert to original',...
    'set all axes to extreme of values in use',...
    'set all axes to specified values'};
subh_dims=[size(ylims,1),size(ylims,2),size(ylims,3)];
ylim_global(1)=min(reshape(ylims(:,:,:,1),prod(subh_dims),1),[],1,'omitnan');
ylim_global(2)=max(reshape(ylims(:,:,:,2),prod(subh_dims),1),[],1,'omitnan');
ylim_spec=ylim_global;
while axis_opt~=0
    disp('options for ordinate');
    for iadj=1:length(axis_opts)
        disp(sprintf('%2.0f->%s',iadj,axis_opts{iadj}));
    end
    disp(' 0->done');
    axis_opt=getinp('choice','d',[0 length(axis_opts)]);
    if axis_opt==3
        ylim_spec=getinp('axis range','f',[-Inf Inf],ylim_spec);
    end
    %
    for irow=1:subh_dims(1)
        for icol=1:subh_dims(2)
            for ipage=1:subh_dims(3)
                if ~isempty(subhandles{irow,icol,ipage})
                    subhandle=subhandles{irow,icol,ipage};
                    switch axis_opt
                        case 1
                            set(subhandle,'YLim',ylims(irow,icol,ipage,:));
                        case 2
                            set(subhandle,'YLim',ylim_global);
                        case 3
                            set(subhandle,'YLim',ylim_spec);
                    end %axis_opt
                end %isempty
            end %ipage
        end %icol
    end %irow
end %while

%psg_dist_heatmaps_demo: show heamaps of distances
%
% dists(:,:,id) is a stack of distances, for each dimension requested,
% stimuli are in same order as in data file
%
% allows for a reordering for display, but this does not affect dists
% 
% display_orders, if present, is a structure, in which each field is theif_us
% name of a custom display order, and the value is a cell array of the 'typenames' in desired order of appearance
% (not all need to be present)
%
%  See also: PSG_GET_COORDSETS, PSG_READ_COORDDATA, COOTODSQ, NICESUBP.
%
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('font_size') font_size=8; end
%
nsets=1;
opts_read=filldefault(opts_read,'if_log',1);
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(opts_read,[],[],nsets);
%
nstims=size(sas{1}.typenames);
disp('stimulus typenames, in file order:')
disp(sas{1}.typenames);
%
if ~exist('display_orders')
    display_orders=struct;
    order_names=cell(0);
    nords=0;
else
    order_names=fieldnames(display_orders);
    nords=length(order_names);
end
disp('-1->order alphabetically');
disp(' 0->order as in file');
for iord=1:nords
    disp(sprintf(' %1.0f->order %s',iord,order_names{iord}));
end
order_choice=getinp('choice','d',[-1 nords],0);
switch order_choice
    case 0
        ip=[1:nstims];
        order_string='file order';
    case -1
        [sorted,ip]=sort(sas{1}.typenames);
        order_string='alphabetical order';
    otherwise
        order_string=cat(2,'order: ',order_names{order_choice});
        ip=[];
        multmatch=0;
        nomatch=0;
        for istim=1:length(display_orders.(order_names{order_choice}))
            ptr=strmatch(display_orders.(order_names{order_choice}){istim},sas{1}.typenames,'exact');
            if length(ptr)>1
                multmatch=multmatch+1;
            elseif ~isempty(ptr)
                ip=[ip,ptr];
            else
                nomatch=nomatch+1;
            end
        end
        disp(sprintf('%2.0f stimuli found, %2.0f not found, %2.0f multiple matches',length(ip),nomatch,multmatch));
end
%
name_string=cat(2,'distances: ',sets{1}.label,' (',order_string,')');
%
nd_list=zeros(1,length(ds{1})); %take into accouint that all dimensions may not be present
for id=1:length(ds{1})
    nd_list(id)=size(ds{1}{id},2);
end
dims_to_show=getinp('one or more dimensions to show','d',[min(nd_list(nd_list>0)) max(nd_list)]);
dims_to_show=dims_to_show(ismember(dims_to_show,nd_list)); %preserve oritinal order
%calculate all distances so they can be scaled together
dists=zeros(sas{1}.nstims,sas{1}.nstims,length(dims_to_show));
for id=1:length(dims_to_show)
    coords=ds{1}{find(nd_list==dims_to_show(id))};
    dists(:,:,id)=sqrt(cootodsq(coords)); 
end
%
[nr,nc]=nicesubp(length(dims_to_show),0.7);
figure;
set(gcf,'Position',[100 100 1200 800]);
set(gcf,'NUmberTitle','off');
set(gcf,'Name',name_string);
for id=1:length(dims_to_show)
    subplot(nr,nc,id)
    imagesc(dists(ip,ip,id),[0 max(dists(:))]);
    title(sprintf('ndims: %2.0f',dims_to_show(id)));
    set(gca,'FontSize',font_size);
    set(gca,'XTick',[1:sas{1}.nstims]);
    set(gca,'XTickLabel',sas{1}.typenames(ip));
    set(gca,'YTick',[1:sas{1}.nstims]);
    set(gca,'YTickLabel',sas{1}.typenames(ip));
    axis square
    colorbar;
end
%
axes('Position',[0.01,0.01,0.01,0.01]); %for text
text(0,0,name_string,'Interpreter','none','FontSize',8);
axis off;

%psg_plotcoords_quads
%
%special=purpose script to visualize results of "quads" experiment
%
%   See also:  PSG_READ_COORDDATA, PSG_PLOTCOORDS.
%
if_connect=getinp('1 to connect bc anchors, 2 to connect all anchors, 3 to connect de anchors (all with mixes)','d',[0 3]);
%
if ~exist('fn_data') fn_data='./psg_data/quads4pt_coords_MC_sess01_10.mat';end
if ~exist('fn_setup') fn_setup='./psg_data/quads4pt9.mat'; end
if ~exist('dim_plot') dim_plot=3; end %choose 2,3, or 4
if ~exist('dim_plot_list') dim_plot_list=[1:dim_plot]; end
%
opts_read=struct;
opts_read.if_log=1;
opts_read.if_auto=1;
%
%local defaults for plots
%
opts_plot=struct;
opts_plot.if_use_rays=0;
opts_plot.label_sets=1;
opts_plot.label_list='typenames';
opts_plot.marker_size=10;
%
%
[d,sa,opts_read_used]=psg_read_coorddata(fn_data,fn_setup,opts_read);
sa.typenames=strrep(sa.typenames,'0400','');
sa.typenames=strrep(sa.typenames,'0200','');
disp('modified typenames');
disp(sa.typenames');
%
figure;
set(gcf,'Position',[100 100 1200 800]);
subplot(1,1,1);
axis_handle=gca;
%
%plot subsets of points, each with its own color
%
color_table.anchors='b';
color_table.anchors_p='b';
color_table.anchors_m='b';
color_table.mixes='r';
color_table.random='g';
if ~exist('subsets')
    subsets=cell(0);
    subsets{1}.members=[17:20];
    subsets{1}.opts.tag_text='anchors_p';
    subsets{1}.opts.marker_noray='+';
    %
    subsets{2}.members=[21:24];
    subsets{2}.opts.tag_text='anchors_m';
    subsets{2}.opts.marker_noray='*';
    %
    subsets{3}.members=[1:16];
    subsets{3}.opts.tag_text='mixes';
    %
    subsets{4}.members=[25];
    subsets{4}.opts.tag_text='random';
    subsets{4}.opts.marker_noray='o';
end
opts_plot_subsets_used=[];
for isub=1:length(subsets)
    subsets{isub}.opts.color_norays=color_table.(subsets{isub}.opts.tag_text); %get color from color_table and tag_text
    opts_plot_subsets=subsets{isub}.opts;
    opts_plot_subsets.if_legend=0; %no legend
    op_fields=fieldnames(opts_plot);
    for k=1:length(op_fields)
        opts_plot_subsets=filldefault(opts_plot_subsets,op_fields{k},opts_plot.(op_fields{k}));
    end
    if isub>1
        opts_plot_subsets.if_just_data=1;
        opts_plot_subsets.tet_show=[];
    end
    sa_subset=sa;
    sa_subset.typenames=sa.typenames(subsets{isub}.members);
    opts_plot_subsets_used{isub}=psg_plotcoords(d{dim_plot}(subsets{isub}.members,:),dim_plot_list,sa_subset,[],opts_plot_subsets);
end
%
%collect tags for legend
%
hc=get(gca,'Children');
hl=cell(0);
ht=[];
for isub=1:length(subsets)
    for k=1:length(hc)
        if 1==min(strfind(hc(k).Tag,subsets{isub}.opts.tag_text))
            hl=[hl;hc(k)];
            ht=strvcat(ht,subsets{isub}.opts.tag_text);
        end
    end
end
if if_connect>0
    %
    %connect some points
    %
    opts_plot_connect=opts_plot;
    opts_plot_connect.if_legend=0;
    opts_plot_connect.tag_text='';
    opts_plot_connect.if_just_data=1;
    opts_plot_connect.connect_only=1;
    opts_plot_connect.color_norays_connect_mode=3; %half each color
    %
    %connect each b or c or d or e anchor with any mix that contains it
    %
    switch if_connect
        case 1
            lets={'b','c'};
            opts_plot_connect.connect_line_type='--';
        case 2
            lets={'b','c','d','e'};
        case 3
            lets={'d','e'};
            opts_plot_connect.connect_line_type=':';
    end
    signs={'p','m'};
    for ilets=1:length(lets)
        if if_connect==2
            if ilets<=2
                opts_plot_connect.connect_line_type='--';
            else
                opts_plot_connect.connect_line_type=':';
            end
        end
        for isigns=1:length(signs)
            lab=cat(2,lets{ilets},signs{isigns});
            select_anchor=strmatch(lab,sa.typenames,'exact');
            select_mix=[];
            for k=1:length(sa.typenames)
                if contains(sa.typenames{k},lab)
                    select_mix=[select_mix,k];
                end
            end
            select_mix=setdiff(select_mix,select_anchor);
            %set up datasets for this correspondence
            dc=zeros(1,dim_plot,1+length(select_mix));
            dc(1,:,1)=d{dim_plot}(select_anchor,:);
            opts_plot_connect.connect_list=[ones(length(select_mix),1),1+[1:length(select_mix)]'];
            opts_plot_connect.color_connect_sets_norays{1}=color_table.anchors;
            for ids=1:length(select_mix)
                dc(1,:,1+ids)=d{dim_plot}(select_mix(ids),:);
                opts_plot_connect.color_connect_sets_norays{ids+1}=color_table.mixes;
            end
            opts_plot_connect_used=psg_plotcoords(dc,dim_plot_list,sa_subset,[],opts_plot_connect);
        end
    end
end %if_connect
%
legend(hl,ht,'Location','NorthEast','Interpreter','none');
%
axes('Position',[0.01,0.04,0.01,0.01]); %for text
text(0,0,fn_data,'Interpreter','none','FontSize',8);
axis off;

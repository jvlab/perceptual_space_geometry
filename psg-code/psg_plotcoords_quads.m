%psg_plotcoords_quads
%
%special=purpose script to visualize results of "quads" experiment
%
%   See also:  PSG_READ_COORDDATA, PSG_PLOTCOORDS.
%
if ~exist('fn_data') fn_data='./psg_data/quads4pt_coords_MC_sess01_01.mat';end
if ~exist('fn_setup') fn_setup='./psg_data/quads4pt9.mat'; end
if ~exist('dim_plot') dim_plot=3; end %choose 2,3, or 4
opts_plot=struct;
opts_plot.if_use_rays=0;
opts_plot.label_sets=1;
opts_plot.label_list='typenames';
%
opts_read=struct;
opts_read.if_log=1;
opts_read.if_auto=1;
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
anchors=[17:25];
sa_anchors=sa;
sa_anchors.typenames=sa.typenames(anchors);
ou=psg_plotcoords(d{dim_plot}(anchors,:),[1:dim_plot],sa_anchors,[],setfields(opts_plot,{'color_norays','axis_handle'},{'b',axis_handle}));
%
hold on;
%
mixes=[1:16];
sa_mixes=sa;
sa_mixes.typenames=sa.typenames(mixes);
ou=psg_plotcoords(d{dim_plot}(mixes,:),[1:dim_plot],sa_mixes,[],setfields(opts_plot,{'color_norays','axis_handle','if_just_data','tet_show'},{'r',axis_handle,1,[]}));
%
axes('Position',[0.01,0.04,0.01,0.01]); %for text
text(0,0,fn_data,'Interpreter','none','FontSize',8);
axis off;

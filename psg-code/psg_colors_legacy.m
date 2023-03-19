function opts_used=psg_colors_legacy(opts_plot)
%simple utility to get legay colors if asked
%
% opts_plot: options structure for psg_plotcoords
%
% opts_used: opts_plot, with 'color' field replaced if requested
%
%  See also:  PSG_PLOTCOORDS, PSG_VISUALIZE_DEMO, PSG_PROCRUSTES_DEMO, PSG_TYPENAMES2COLORS.
%
colors_legacy.g=[.3 .3 .3];
colors_legacy.b=[1 0 0];
colors_legacy.c=[0 .7 0];
colors_legacy.a=[0 0 1];
if (getinp('1 to use legacy colors for g,b,c,a','d',[0 1],0))
    opts_plot=setfield(opts_plot,'colors',colors_legacy);
end
opts_used=opts_plot;
return

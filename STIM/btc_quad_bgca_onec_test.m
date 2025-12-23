%btc_quad_bgca_onec_test: test calculation of change in image states when
%some white or black checks are modified, as in btc_quad+bgca_onec_make.
%
% See btc_quad_onecomponent_notes.docx for details.
%
%   See also:  BTC_DEFINE, BTC_AUGCOORDS, BTC_MAKEMAPS, BTC_QUAD_BGCA_MAKE.
%
dict=btc_define;
%
if ~exist('mapsize') mapsize=1024; end
%
if ~exist('h') h=0.23; end %fraction of checks to flip, not simple
if ~exist('y') y=0.43; end %nonzero param value, large and not simple
%
opts_component=struct;
opts_component.area=[mapsize,mapsize];
%
%add to a random map
disp(' initial map: all params zero');
g=0;
b=0;
c=0;
d=0;
a=0;
qp=b;
rp=0;
g_pert=h;
b_pert=h^2+b*(1-h).^2;
c_pert=h^2+c*(1-h).^2;
d_pert=h^2+d*(1-h).^2;
t_pert=NaN;
a_pert=h^4+(qp+3*rp/4)*(1-h).^2*h^2+a*(1-h)^4;
stats_pert_ideal=[g_pert b_pert c_pert d_pert d_pert t_pert t_pert t_pert t_pert a_pert]';
s=struct;
aug=btc_augcoords(s,dict);
map_init=btc_makemaps(aug.method{1},opts_component,dict);
%
stats_init_ideal=aug.method{1}.vec(:);
counts=btc_map2counts(map_init);
stats_init_found=btc_corrs2vec(getcorrs_p2x2(counts/sum(counts(:))))';
mask=double(rand(mapsize,mapsize)<h); %|g| of these are white
map_pert=map_init;
mc_flip=find(map_init(:)==0);
map_pert(mc_flip)=mask(mc_flip);
counts_pert=btc_map2counts(map_pert);
stats_pert_found=btc_corrs2vec(getcorrs_p2x2(counts_pert/sum(counts_pert(:))))';
disp('              initial stats          perturbed stats');
disp('     ideal    found    fnd-idl     ideal    found   fnd-idl');
disp([stats_init_ideal stats_init_found stats_init_found-stats_init_ideal stats_pert_ideal stats_pert_found stats_pert_found-stats_pert_ideal]);
%
%add to a map with just b
%
disp(' initial map: b is nonzero');
g=0;
b=y;
c=0;
d=0;
a=0;
qp=b;
rp=0;
g_pert=h;
b_pert=h^2+b*(1-h).^2;
c_pert=h^2+c*(1-h).^2;
d_pert=h^2+d*(1-h).^2;
t_pert=NaN;
a_pert=h^4+(qp+3*rp/4)*(1-h).^2*h^2+a*(1-h)^4;
stats_pert_ideal=[g_pert b_pert c_pert d_pert d_pert t_pert t_pert t_pert t_pert a_pert]';
s=struct;
s.b=b;
s.a=0;
aug=btc_augcoords(s,dict);
map_init=btc_makemaps(aug.method{1},opts_component,dict);
%
stats_init_ideal=aug.method{1}.vec(:);
counts=btc_map2counts(map_init);
stats_init_found=btc_corrs2vec(getcorrs_p2x2(counts/sum(counts(:))))';
mask=double(rand(mapsize,mapsize)<h); %|g| of these are white
map_pert=map_init;
mc_flip=find(map_init(:)==0);
map_pert(mc_flip)=mask(mc_flip);
counts_pert=btc_map2counts(map_pert);
stats_pert_found=btc_corrs2vec(getcorrs_p2x2(counts_pert/sum(counts_pert(:))))';
disp('              initial stats          perturbed stats');
disp('     ideal    found    fnd-idl     ideal    found   fnd-idl');
disp([stats_init_ideal stats_init_found stats_init_found-stats_init_ideal stats_pert_ideal stats_pert_found stats_pert_found-stats_pert_ideal]);
%
%add to a map with just d
%
disp(' initial map: d is nonzero');
g=0;
b=0;
c=0;
d=y;
ep=0;
a=0;
qp=b+c;
rp=d+ep;
g_pert=h;
b_pert=h^2+b*(1-h).^2;
c_pert=h^2+c*(1-h).^2;
d_pert=h^2+d*(1-h).^2;
e_pert=h^2+ep*(1-h).^2;
t_pert=NaN;
a_pert=h^4+(qp+3*rp/4)*(1-h).^2*h^2+a*(1-h)^4;
stats_pert_ideal=[g_pert b_pert c_pert d_pert e_pert t_pert t_pert t_pert t_pert a_pert]';
s=struct;
s.d=d;
s.a=0;
aug=btc_augcoords(s,dict);
map_init=btc_makemaps(aug.method{1},opts_component,dict);
%
stats_init_ideal=aug.method{1}.vec(:);
counts=btc_map2counts(map_init);
stats_init_found=btc_corrs2vec(getcorrs_p2x2(counts/sum(counts(:))))';
mask=double(rand(mapsize,mapsize)<h); %|g| of these are white
map_pert=map_init;
mc_flip=find(map_init(:)==0);
map_pert(mc_flip)=mask(mc_flip);
counts_pert=btc_map2counts(map_pert);
stats_pert_found=btc_corrs2vec(getcorrs_p2x2(counts_pert/sum(counts_pert(:))))';
disp('              initial stats          perturbed stats');
disp('     ideal    found    fnd-idl     ideal    found   fnd-idl');
disp([stats_init_ideal stats_init_found stats_init_found-stats_init_ideal stats_pert_ideal stats_pert_found stats_pert_found-stats_pert_ideal]);
%
%add to a map with just a
%
disp(' initial map: a is nonzero');
g=0;
b=0;
c=0;
d=0;
ep=0;
a=y;
qp=b+c;
rp=d+ep;
g_pert=h;
b_pert=h^2+b*(1-h).^2;
c_pert=h^2+c*(1-h).^2;
d_pert=h^2+d*(1-h).^2;
t_pert=NaN;
a_pert=h^4+(qp+3*rp/4)*(1-h).^2*h^2+a*(1-h)^4;
stats_pert_ideal=[g_pert b_pert c_pert d_pert d_pert t_pert t_pert t_pert t_pert a_pert]';
s=struct;
s.a=a;
aug=btc_augcoords(s,dict);
map_init=btc_makemaps(aug.method{1},opts_component,dict);
%
stats_init_ideal=aug.method{1}.vec(:);
counts=btc_map2counts(map_init);
stats_init_found=btc_corrs2vec(getcorrs_p2x2(counts/sum(counts(:))))';
mask=double(rand(mapsize,mapsize)<h); %|g| of these are white
map_pert=map_init;
mc_flip=find(map_init(:)==0);
map_pert(mc_flip)=mask(mc_flip);
counts_pert=btc_map2counts(map_pert);
stats_pert_found=btc_corrs2vec(getcorrs_p2x2(counts_pert/sum(counts_pert(:))))';
disp('              initial stats          perturbed stats');
disp('     ideal    found    fnd-idl     ideal    found   fnd-idl');
disp([stats_init_ideal stats_init_found stats_init_found-stats_init_ideal stats_pert_ideal stats_pert_found stats_pert_found-stats_pert_ideal]);



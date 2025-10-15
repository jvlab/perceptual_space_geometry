%spokes_setup_create
%script to create some standard spoke setups, and also (after 14Oct25), combinations of 
% spokes setups and setups with stimuli from a resource, e.g., created with
% btc_quad_bcde_make
%
% 08Dec22: add setup 13, two axes each polarity, six levels and cmax_sets
% 05Jun23: add setups with 12 and 24 spokes (14-15), and 5x5 grid (16)
% 06Jun23: add setups with 5x5 grid in a quadrant (17-20); btc_pair_choices defined
% 07Jun23: add combinations to btc_pair_choices
% 29Jan24: add combinations to setup 2, four axes, each polarity
% 20Jan25: add combinations gtva, detv, btuv to setup 2
% 14Oct25: add setup with resource with quad specification for bcde (setup 21)
%
%   See also:  PSG_SPOKES_SETUP, SPOKES_LAYOUT_DEMO, BTC_QUAD_BCDE_MAKE.

if ~exist('cmax_sets')
    cmax_sets=cell(1);
    cmax_sets{1}.desc='standard for thresholds, low';
    cmax_sets{1}.vals=[0.2 0.4 0.4 0.5 0.5 1.0 1.0 1.0 1.0 0.8];
    cmax_sets{2}.desc='standard for thresholds';
    cmax_sets{2}.vals=[0.4 0.6 0.6 0.7 0.7 1.0 1.0 1.0 1.0 1.0];
    cmax_sets{3}.desc='standard for thresholds, high';
    cmax_sets{3}.vals=[0.6 0.9 0.9 0.9 0.9 1.0 1.0 1.0 1.0 1.0];
    cmax_sets{4}.desc='standard for thresholds but beta-diags match beta-cards';
    cmax_sets{4}.vals=[0.4 0.6 0.6 0.6 0.6 1.0 1.0 1.0 1.0 1.0];
end
%
pm={'p','m'};
%
btc_pair_choices{1}={'b','c'};
btc_pair_choices{2}={'d','e'};
btc_pair_choices{3}={'g','b'};
btc_pair_choices{4}={'t','v'};
btc_pair_choices{5}={'u','w'};
btc_pair_choices{6}={'g','t'};
btc_pair_choices{7}={'g','a'};
btc_pair_choices{8}={'b','a'};
btc_pair_choices{9}={'t','a'};
%
isetup=1;
spoke_setups{isetup}.name='two  axes, each polarity, mix with 1';
spoke_setups{isetup}.ndims=2; %choose two dimensions
spoke_setups{isetup}.btc_choices=btc_pair_choices;
spoke_setups{isetup}.nspokes=8;
angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
spoke_setups{isetup}.endpoints=sign(round(2*[cos(angs)',sin(angs')]));
spoke_setups{isetup}.mixing=spoke_setups{isetup}.endpoints; %mixing as in a 2d plane
spoke_setups{isetup}.nclevs=3;
%
isetup=2;
spoke_setups{isetup}.name='four axes, each polarity';
spoke_setups{isetup}.ndims=4; %choose four dimensions
spoke_setups{isetup}.btc_choices{1}={'b','d','c','e'};
spoke_setups{isetup}.btc_choices{2}={'b','g','c','a'};
spoke_setups{isetup}.btc_choices{4}={'t','u','v','w'}; %bug fix 20Jan25
spoke_setups{isetup}.btc_choices{3}={'d','g','e','a'};
spoke_setups{isetup}.btc_choices{5}={'g','t','v','a'};
spoke_setups{isetup}.btc_choices{5}={'g','t','v','a'};
spoke_setups{isetup}.btc_choices{6}={'d','e','t','v'};
spoke_setups{isetup}.btc_choices{7}={'b','t','u','v'};
spoke_setups{isetup}.nspokes=8;
angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
spoke_setups{isetup}.endpoints=[cos(angs)',sin(angs')];
spoke_setups{isetup}.mixing=[eye(spoke_setups{isetup}.ndims);-eye(spoke_setups{isetup}.ndims)]; %no mixing; each spoke is an axis
spoke_setups{isetup}.nclevs=3;
%
isetup=3;
spoke_setups{isetup}.name='six  axes, each polarity';
spoke_setups{isetup}.ndims=6; %choose six dimensions
spoke_setups{isetup}.btc_choices{1}={'g','b','c','d','e','a'};
spoke_setups{isetup}.btc_choices{2}={'b','c','t','u','v','w'};
spoke_setups{isetup}.nspokes=12;
angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
spoke_setups{isetup}.endpoints=[cos(angs)',sin(angs')];
spoke_setups{isetup}.mixing=[eye(spoke_setups{isetup}.ndims);-eye(spoke_setups{isetup}.ndims)]; %no mixing; each spoke is an axis
spoke_setups{isetup}.nclevs=3;
%
isetup=4;
spoke_setups{isetup}.name='nine axes, each polarity';
spoke_setups{isetup}.ndims=9; %choose six dimensions
spoke_setups{isetup}.btc_choices{1}={'g','b','c','d','e','t','u','v','w'};
spoke_setups{isetup}.btc_choices{2}={'b','c','d','e','t','u','v','w','a'};
spoke_setups{isetup}.nspokes=18;
angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
spoke_setups{isetup}.endpoints=[cos(angs)',sin(angs')];
spoke_setups{isetup}.mixing=[eye(spoke_setups{isetup}.ndims);-eye(spoke_setups{isetup}.ndims)]; %no mixing; each spoke is an axis
spoke_setups{isetup}.nclevs=2;
%
isetup=5;
spoke_setups{isetup}.name='two axes, each polarity, same-sign mixture, mix with 1/sqrt2';
spoke_setups{isetup}.ndims=2; %choose two dimensions
spoke_setups{isetup}.btc_choices{1}={'b','c'};
spoke_setups{isetup}.btc_choices{2}={'d','e'};
spoke_setups{isetup}.btc_choices{3}={'g','b'};
%    spoke_setups{isetup}.btc_choices{4}={'g','a'}; %out of domain with gmax=0.2, amax=0.8
%    spoke_setups{isetup}.btc_choices{5}={'b','a'}; %out of domain with bmax=0.4, amax=0.8
spoke_setups{isetup}.nspokes=6;
angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
spoke_setups{isetup}.mixing=[1 0;1/sqrt(2) 1/sqrt(2);0 1;-1 0;-1/sqrt(2) -1/sqrt(2);0 -1]; %mix one pair of axes in only same-sign directions
spoke_setups{isetup}.endpoints=sign(spoke_setups{isetup}.mixing);
spoke_setups{isetup}.nclevs=3;
%
isetup=6;
spoke_setups{isetup}=spoke_setups{isetup-1};
spoke_setups{isetup}.name='two axes, each polarity, opposite-sign mixture, mix with 1/sqrt2';
spoke_setups{isetup}.mixing=[1 0;0 1;-1/sqrt(2) 1/sqrt(2);-1 0;0 -1;1/sqrt(2) -1/sqrt(2)]; %mix one pair of axes in only opposite-sign directions
spoke_setups{isetup}.endpoints=sign(spoke_setups{isetup}.mixing);
%
isetup=7;
spoke_setups{isetup}=spoke_setups{isetup-2};
spoke_setups{isetup}.name='two axes, each polarity, same-sign mixture, mix with 1/2';
spoke_setups{isetup}.btc_choices{4}={'g','a'};
spoke_setups{isetup}.btc_choices{5}={'b','a'};
spoke_setups{isetup}.mixing=[1 0;1/2 1/2;0 1;-1 0;-1/2 -1/2;0 -1]; %mix one pair of axes in only same-sign directions
%
isetup=8;
spoke_setups{isetup}=spoke_setups{isetup-1};
spoke_setups{isetup}.name='two axes, each polarity, opposite-sign mixture, mix with 1/2';
spoke_setups{isetup}.mixing=[1 0;0 1;-1/2 1/2;-1 0;0 -1;1/2 -1/2]; %mix one pair of axes in only opposite-sign directions
spoke_setups{isetup}.endpoints=sign(spoke_setups{isetup}.mixing);
%
isetup=9;
spoke_setups{isetup}=spoke_setups{1};
spoke_setups{isetup}.name='two  axes, each polarity, mix with 1/sqrt(2)';
angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
spoke_setups{isetup}.mixing=[cos(angs)',sin(angs')];
%
isetup=10;
spoke_setups{isetup}=spoke_setups{isetup-1};
spoke_setups{isetup}.name='two  axes, each polarity, mix with 1/2';
spoke_setups{isetup}.btc_choices{6}={'g','a'};
spoke_setups{isetup}.btc_choices{7}={'b','a'};
spoke_setups{isetup}.mixing=round(2*spoke_setups{isetup}.mixing)/2; %change 1/sqrt(2) to 1/2
%
isetup=11;
spoke_setups{isetup}.name='two  axes, each polarity, two same-sign mixes';
spoke_setups{isetup}.ndims=2; %choose two dimensions
spoke_setups{isetup}.btc_choices{1}={'b','c'};
spoke_setups{isetup}.btc_choices{2}={'g','b'};
spoke_setups{isetup}.btc_choices{3}={'g','c'};
spoke_setups{isetup}.nspokes=8;
angs=2.*pi*[0:3 6:9]/12;
spoke_setups{isetup}.endpoints=[cos(angs)',sin(angs')];
spoke_setups{isetup}.mixing=[cos(angs)',sin(angs')];
spoke_setups{isetup}.nclevs=3;
%
isetup=12;
spoke_setups{isetup}=spoke_setups{isetup-1};
spoke_setups{isetup}.name='two  axes, each polarity, two oppo-sign mixes';
angs=2.*pi*[3:6 9:12]/12;
spoke_setups{isetup}.endpoints=[cos(angs)',sin(angs')];
spoke_setups{isetup}.mixing=[cos(angs)',sin(angs')];
spoke_setups{isetup}.nclevs=3;
%
isetup=13;
spoke_setups{isetup}.name='two  axes, each polarity';
spoke_setups{isetup}.ndims=2; %choose four dimensions
spoke_setups{isetup}.btc_choices{1}={'b','c'};
spoke_setups{isetup}.btc_choices{2}={'d','e'};
spoke_setups{isetup}.btc_choices{3}={'t','v'};
spoke_setups{isetup}.btc_choices{4}={'t','u'};
spoke_setups{isetup}.btc_choices{5}={'g','b'};
spoke_setups{isetup}.btc_choices{6}={'g','c'};
spoke_setups{isetup}.btc_choices{7}={'g','d'};
spoke_setups{isetup}.btc_choices{8}={'g','e'};
spoke_setups{isetup}.btc_choices{9}={'g','a'};
spoke_setups{isetup}.nspokes=4;
angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
spoke_setups{isetup}.endpoints=[cos(angs)',sin(angs')];
spoke_setups{isetup}.mixing=[eye(spoke_setups{isetup}.ndims);-eye(spoke_setups{isetup}.ndims)]; %no mixing; each spoke is an axis
spoke_setups{isetup}.nclevs=6;
%
isetup=14;
spoke_setups{isetup}.name='two  axes, 12 spokes';
spoke_setups{isetup}.ndims=2; %choose two dimensions
spoke_setups{isetup}.btc_choices=btc_pair_choices;
spoke_setups{isetup}.nspokes=12;
angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
spoke_setups{isetup}.endpoints=[cos(angs)',sin(angs')];
spoke_setups{isetup}.mixing=spoke_setups{isetup}.endpoints; %mixing as in a 2d plane
spoke_setups{isetup}.nclevs=2;
%
isetup=15;
spoke_setups{isetup}.name='two  axes, 24 spokes';
spoke_setups{isetup}.ndims=2; %choose two dimensions
spoke_setups{isetup}.btc_choices=btc_pair_choices;
spoke_setups{isetup}.nspokes=24;
angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
spoke_setups{isetup}.endpoints=[cos(angs)',sin(angs')];
spoke_setups{isetup}.mixing=spoke_setups{isetup}.endpoints; %mixing as in a 2d plane
spoke_setups{isetup}.nclevs=1;
%
isetup=16;
spoke_setups{isetup}.name='two  axes, 5x5 grid';
spoke_setups{isetup}.ndims=2; %choose two dimensions
spoke_setups{isetup}.btc_choices=btc_pair_choices;
spoke_setups{isetup}.nspokes=24; %(kludge:  actually 4 cardinal, 4 oblique, and 8 "knight move")
[ex,ey]=meshgrid(-1:0.5:1);
ptr_origin=(1+size(ex(:),1))/2;
not_origin=setdiff(1:size(ex(:),1),ptr_origin);
ex=ex(not_origin); %remove origin
ey=ey(not_origin); %remove origin
spoke_setups{isetup}.endpoints=[ex(:),ey(:)];
spoke_setups{isetup}.mixing=spoke_setups{isetup}.endpoints; %mixing as in a 2d plane
spoke_setups{isetup}.nclevs=1;
%
for id1=1:2
    for id2=1:2
        isetup=isetup+1;
        spoke_setups{isetup}.name=sprintf('  two axes, quadrant, ax 1: %s ax 2: %s',pm{id1},pm{id2});
        spoke_setups{isetup}.btc_choices=btc_pair_choices;
        spoke_setups{isetup}.nspokes=24;
        [ex,ey]=meshgrid(0:0.25:1);
        ex=ex-id1+1;
        ey=ey-id2+1;
        ex=ex(:);
        ey=ey(:);
        ptr_origin=intersect(find(ex==0),find(ey==0));
        not_origin=setdiff(1:size(ex(:),1),ptr_origin);
        ex=ex(not_origin); %remove origin
        ey=ey(not_origin); %remove origin
        spoke_setups{isetup}.endpoints=[ex(:),ey(:)];
        spoke_setups{isetup}.mixing=spoke_setups{isetup}.endpoints; %mixing as in a 2d plane
        spoke_setups{isetup}.nclevs=1;      
    end
end
%
isetup=isetup+1; %isetup=21
spoke_setups{isetup}.name='four axes and 16 quad combinations';
spoke_setups{isetup}.ndims=4; %choose four dimensions
spoke_setups{isetup}.btc_choices={{'b','c','d','e'}};
spoke_setups{isetup}.nspokes=8;
angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
spoke_setups{isetup}.endpoints=[cos(angs)',sin(angs')];
spoke_setups{isetup}.mixing=[eye(spoke_setups{isetup}.ndims);-eye(spoke_setups{isetup}.ndims)]; %no mixing; each spoke is an axis
spoke_setups{isetup}.nclevs=1;
spoke_setups{isetup}.clev_fracvals_def=2/3; %for singletons, 2/3 of the way to the max value of a 4-axis, 3-point experiment
spoke_setups{isetup}.resource_template='../stim/btc_quad_bcde_make_t20_sc4_it5000.mat';
spoke_setups{isetup}.need_resource=1;
%
disp(sprintf('%3.0f setups made',length(spoke_setups)));


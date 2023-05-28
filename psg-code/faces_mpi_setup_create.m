%faces_mpi_setup_create
%script to create some face setups for mpi dataset
%
% modeled after spokes_setup_create.
%
%   See also:  FACES_MPI_PSG_SETUP.

if ~exist('faces_mpi_setups')
    faces_mpi_setups=cell(0);
end
% if ~exist('cmax_sets')
%     cmax_sets=cell(1);
%     cmax_sets{1}.desc='standard for thresholds';
%     cmax_sets{1}.vals=[0.2 0.4 0.4 0.5 0.5 1.0 1.0 1.0 1.0 0.8];
%     cmax_sets{2}.desc='standard for thresholds';
%     cmax_sets{2}.vals=[0.4 0.6 0.6 0.7 0.7 1.0 1.0 1.0 1.0 1.0];
%     cmax_sets{3}.desc='standard for thresholds, high';
%     cmax_sets{3}.vals=[0.6 0.9 0.9 0.9 0.9 1.0 1.0 1.0 1.0 1.0];
%     cmax_sets{4}.desc='standard for thresholds but beta-diags match beta-cards';
%     cmax_sets{4}.vals=[0.4 0.6 0.6 0.6 0.6 1.0 1.0 1.0 1.0 1.0];
% end
% %
% isetup=1;
% spoke_setups{isetup}.name='two  axes, each polarity, mix with 1';
% spoke_setups{isetup}.ndims=2; %choose two dimensions
% spoke_setups{isetup}.btc_choices{1}={'b','c'};
% spoke_setups{isetup}.btc_choices{2}={'d','e'};
% spoke_setups{isetup}.btc_choices{3}={'g','b'};
% spoke_setups{isetup}.btc_choices{4}={'t','v'};
% spoke_setups{isetup}.btc_choices{5}={'u','w'};
% spoke_setups{isetup}.nspokes=8;
% angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
% spoke_setups{isetup}.endpoints=sign(round(2*[cos(angs)',sin(angs')]));
% spoke_setups{isetup}.mixing=spoke_setups{isetup}.endpoints; %mixing as in a 2d plane
% spoke_setups{isetup}.nclevs=3;
% %
% isetup=2;
% spoke_setups{isetup}.name='four axes, each polarity';
% spoke_setups{isetup}.ndims=4; %choose four dimensions
% spoke_setups{isetup}.btc_choices{1}={'b','d','c','e'};
% spoke_setups{isetup}.btc_choices{2}={'b','g','c','a'};
% spoke_setups{isetup}.btc_choices{3}={'t','u','v','w'};
% spoke_setups{isetup}.nspokes=8;
% angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
% spoke_setups{isetup}.endpoints=[cos(angs)',sin(angs')];
% spoke_setups{isetup}.mixing=[eye(spoke_setups{isetup}.ndims);-eye(spoke_setups{isetup}.ndims)]; %no mixing; each spoke is an axis
% spoke_setups{isetup}.nclevs=3;
% %
% isetup=3;
% spoke_setups{isetup}.name='six  axes, each polarity';
% spoke_setups{isetup}.ndims=6; %choose six dimensions
% spoke_setups{isetup}.btc_choices{1}={'g','b','c','d','e','a'};
% spoke_setups{isetup}.btc_choices{2}={'b','c','t','u','v','w'};
% spoke_setups{isetup}.nspokes=12;
% angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
% spoke_setups{isetup}.endpoints=[cos(angs)',sin(angs')];
% spoke_setups{isetup}.mixing=[eye(spoke_setups{isetup}.ndims);-eye(spoke_setups{isetup}.ndims)]; %no mixing; each spoke is an axis
% spoke_setups{isetup}.nclevs=3;
% %
% isetup=4;
% spoke_setups{isetup}.name='nine axes, each polarity';
% spoke_setups{isetup}.ndims=9; %choose six dimensions
% spoke_setups{isetup}.btc_choices{1}={'g','b','c','d','e','t','u','v','w'};
% spoke_setups{isetup}.btc_choices{2}={'b','c','d','e','t','u','v','w','a'};
% spoke_setups{isetup}.nspokes=18;
% angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
% spoke_setups{isetup}.endpoints=[cos(angs)',sin(angs')];
% spoke_setups{isetup}.mixing=[eye(spoke_setups{isetup}.ndims);-eye(spoke_setups{isetup}.ndims)]; %no mixing; each spoke is an axis
% spoke_setups{isetup}.nclevs=2;
% %
% isetup=5;
% spoke_setups{isetup}.name='two axes, each polarity, same-sign mixture, mix with 1/sqrt2';
% spoke_setups{isetup}.ndims=2; %choose two dimensions
% spoke_setups{isetup}.btc_choices{1}={'b','c'};
% spoke_setups{isetup}.btc_choices{2}={'d','e'};
% spoke_setups{isetup}.btc_choices{3}={'g','b'};
% %    spoke_setups{isetup}.btc_choices{4}={'g','a'}; %out of domain with gmax=0.2, amax=0.8
% %    spoke_setups{isetup}.btc_choices{5}={'b','a'}; %out of domain with bmax=0.4, amax=0.8
% spoke_setups{isetup}.nspokes=6;
% angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
% spoke_setups{isetup}.mixing=[1 0;1/sqrt(2) 1/sqrt(2);0 1;-1 0;-1/sqrt(2) -1/sqrt(2);0 -1]; %mix one pair of axes in only same-sign directions
% spoke_setups{isetup}.endpoints=sign(spoke_setups{isetup}.mixing);
% spoke_setups{isetup}.nclevs=3;
% %
% isetup=6;
% spoke_setups{isetup}=spoke_setups{isetup-1};
% spoke_setups{isetup}.name='two axes, each polarity, opposite-sign mixture, mix with 1/sqrt2';
% spoke_setups{isetup}.mixing=[1 0;0 1;-1/sqrt(2) 1/sqrt(2);-1 0;0 -1;1/sqrt(2) -1/sqrt(2)]; %mix one pair of axes in only opposite-sign directions
% spoke_setups{isetup}.endpoints=sign(spoke_setups{isetup}.mixing);
% %
% isetup=7;
% spoke_setups{isetup}=spoke_setups{isetup-2};
% spoke_setups{isetup}.name='two axes, each polarity, same-sign mixture, mix with 1/2';
% spoke_setups{isetup}.btc_choices{4}={'g','a'};
% spoke_setups{isetup}.btc_choices{5}={'b','a'};
% spoke_setups{isetup}.mixing=[1 0;1/2 1/2;0 1;-1 0;-1/2 -1/2;0 -1]; %mix one pair of axes in only same-sign directions
% %
% isetup=8;
% spoke_setups{isetup}=spoke_setups{isetup-1};
% spoke_setups{isetup}.name='two axes, each polarity, opposite-sign mixture, mix with 1/2';
% spoke_setups{isetup}.mixing=[1 0;0 1;-1/2 1/2;-1 0;0 -1;1/2 -1/2]; %mix one pair of axes in only opposite-sign directions
% spoke_setups{isetup}.endpoints=sign(spoke_setups{isetup}.mixing);
% %
% isetup=9;
% spoke_setups{isetup}=spoke_setups{1};
% spoke_setups{isetup}.name='two  axes, each polarity, mix with 1/sqrt(2)';
% angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
% spoke_setups{isetup}.mixing=[cos(angs)',sin(angs')];
% %
% isetup=10;
% spoke_setups{isetup}=spoke_setups{isetup-1};
% spoke_setups{isetup}.name='two  axes, each polarity, mix with 1/2';
% spoke_setups{isetup}.btc_choices{6}={'g','a'};
% spoke_setups{isetup}.btc_choices{7}={'b','a'};
% spoke_setups{isetup}.mixing=round(2*spoke_setups{isetup}.mixing)/2; %change 1/sqrt(2) to 1/2
% %
% isetup=11;
% spoke_setups{isetup}.name='two  axes, each polarity, two same-sign mixes';
% spoke_setups{isetup}.ndims=2; %choose two dimensions
% spoke_setups{isetup}.btc_choices{1}={'b','c'};
% spoke_setups{isetup}.btc_choices{2}={'g','b'};
% spoke_setups{isetup}.btc_choices{3}={'g','c'};
% spoke_setups{isetup}.nspokes=8;
% angs=2.*pi*[0:3 6:9]/12;
% spoke_setups{isetup}.endpoints=[cos(angs)',sin(angs')];
% spoke_setups{isetup}.mixing=[cos(angs)',sin(angs')];
% spoke_setups{isetup}.nclevs=3;
% %
% isetup=12;
% spoke_setups{isetup}=spoke_setups{isetup-1};
% spoke_setups{isetup}.name='two  axes, each polarity, two oppo-sign mixes';
% angs=2.*pi*[3:6 9:12]/12;
% spoke_setups{isetup}.endpoints=[cos(angs)',sin(angs')];
% spoke_setups{isetup}.mixing=[cos(angs)',sin(angs')];
% spoke_setups{isetup}.nclevs=3;
% %
% isetup=13;
% spoke_setups{isetup}.name='two  axes, each polarity';
% spoke_setups{isetup}.ndims=2; %choose four dimensions
% spoke_setups{isetup}.btc_choices{1}={'b','c'};
% spoke_setups{isetup}.btc_choices{2}={'d','e'};
% spoke_setups{isetup}.btc_choices{3}={'t','v'};
% spoke_setups{isetup}.btc_choices{4}={'t','u'};
% spoke_setups{isetup}.btc_choices{5}={'g','b'};
% spoke_setups{isetup}.btc_choices{6}={'g','c'};
% spoke_setups{isetup}.btc_choices{7}={'g','d'};
% spoke_setups{isetup}.btc_choices{8}={'g','e'};
% spoke_setups{isetup}.btc_choices{9}={'g','a'};
% spoke_setups{isetup}.nspokes=4;
% angs=2.*pi*[0:spoke_setups{isetup}.nspokes-1]/spoke_setups{isetup}.nspokes;
% spoke_setups{isetup}.endpoints=[cos(angs)',sin(angs')];
% spoke_setups{isetup}.mixing=[eye(spoke_setups{isetup}.ndims);-eye(spoke_setups{isetup}.ndims)]; %no mixing; each spoke is an axis
% spoke_setups{isetup}.nclevs=6;
% 

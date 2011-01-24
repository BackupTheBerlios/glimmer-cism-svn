
%% make a simple domain w/ a square of non-zero thickness ... for testing 
%% advection schemes

clear all

thck = 10*ones( 101, 101 );
thck( 5:25, 5:25 ) = 200;
% thck( 90:95,90:95) = 150;

usrf = thck;
topg = zeros( size( thck ) );

figure, imagesc( thck ), axis xy
axis equal

uvelhom = 10*ones(5,100,100);
vvelhom = uvelhom;

uvelbc = uvelhom; vvelbc = vvelhom;

kinbcmask = zeros( 100, 100 );
kinbcmask(:,:) = 1;
% kinbcmask(1:80,1:80) = 1;

save advection.mat thck topg usrf uvelhom vvelhom kinbcmask uvelbc vvelbc




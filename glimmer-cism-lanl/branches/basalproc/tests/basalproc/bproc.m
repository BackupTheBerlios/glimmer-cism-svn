
%% file to create .mat input variables for Marion's coupled runs
%%
%% To be read by .py script to make .nc input file

rho_i = 910;     % ice density kg/m^3
rho_w = 1028;    % ocean water density kg/m^3

H0 = 1e3;
T0 = -10;
q0 = -7e-2;
slope = 1e-4;
b0 = 0.5;

r = 22;
c = 45;
l = 11;     % no of vert levels

dew = 5e3; 
dns = 5e3;

x = [0:5e3:5e3*c-1];
y = [0:5e3:5e3*r-1]';

thck = H0 * ones( r, c );
airt = T0 * ones( r, c );
qgeo = q0 * ones( r, c );
acab = b0 * ones( r, c );

% topg = repmat( linspace( 1e3, 1e3-(slope*(c-1)*dew), c), r, 1 );

temp = repmat( linspace( 0, -H0*(rho_i/rho_w)-50, c-3), r, 1 ) + 0;        %% force to be near floatationg at ds end
topg = [ temp, repmat(temp(:,end-3), 1, 3 ) ];

usrf = topg + thck;

% %% put buffer of zero thickness cells around perimeter (for remapping)
thck(1:3,:) = 0; thck(:,1:3) = 0; thck(end-2:end,:) = 0; thck(:,end-2:end) = 0;
usrf(1:3,:) = 0; usrf(:,1:3) = 0; usrf(end-2:end,:) = 0; usrf(:,end-2:end) = 0;

figure(1), imagesc( x/1e3, y/1e3, thck ), axis xy, axis equal, axis tight, colorbar
xlabel( 'x (km)' ), ylabel( 'y (km)' ), title( 'thickness (m)' )

figure(2), imagesc( x/1e3, y/1e3, topg ), axis xy, axis equal, axis tight, colorbar
xlabel( 'x (km)' ), ylabel( 'y (km)' ), title( 'basal topg (m)' )

figure(3), imagesc( x/1e3, y/1e3, usrf ), axis xy, axis equal, axis tight, colorbar
xlabel( 'x (km)' ), ylabel( 'y (km)' ), title( 'upper surface (m)' )


%% load in old till map
% load ~/Home/GLAM/GLIMGLAM/SENS/new_UPB/trunk/GLAM/Tillggl           % Marion's path
load ~/work/modeling/glam-stream-marion-new/trunk/GLAM/Tillggl      % Steve's path

minTauf = Tillggl;
% minTauf = 5e3 * ones( size( minTauf ) );       % for debugging
beta = 5e1*ones(size(minTauf));
% beta(8:14,:) = 3e1;
% beta(10:12,:) = 0.5e1;

figure(4), imagesc( x/1e3, y/1e3, minTauf/1e3 ), axis xy, axis equal, axis tight, colorbar
xlabel( 'x (km)' ), ylabel( 'y (km)' ), title( 'Tau0 (kPa)' )

figure(5), imagesc( x/1e3, y/1e3, airt ), axis xy, axis equal, axis tight, colorbar
xlabel( 'x (km)' ), ylabel( 'y (km)' ), title( 'airt temp (C)' )

figure(6), imagesc( x/1e3, y/1e3, qgeo ), axis xy, axis equal, axis tight, colorbar
xlabel( 'x (km)' ), ylabel( 'y (km)' ), title( 'geo flux (W m^2)' )

%% add a kinbcmask field to specify where 0 flux bcs are
kinbcmask = zeros(size(minTauf));
kinbcmask(1:3,:) = 1; kinbcmask(end-2:end,:) = 1; kinbcmask(:,1:3) = 1;

uvelhom = zeros( l, r-1, c-1 );
for i=1:l
    uvelhom(i,:,:) = zeros(size(minTauf));
end

vvelhom = uvelhom; 

figure(7), imagesc( x/1e3, y/1e3, kinbcmask ), axis xy, axis equal, axis tight, colorbar
xlabel( 'x (km)' ), ylabel( 'y (km)' ), title( 'kinbcmask' )

figure(8), imagesc( x/1e3, y/1e3, acab ), axis xy, axis equal, axis tight, colorbar
xlabel( 'x (km)' ), ylabel( 'y (km)' ), title( 'acab (m/a)' )

tauf=minTauf;

% cd ~/Home/Glimmer2/glimmer-cism-lanl/branches/basalproc/tests/basalproc     % Marion's path
cd /Users/sprice/work/modeling/cism_new/branches/tests/basalproc            % Steve's path

save bproc.mat airt acab qgeo usrf topg thck beta tauf kinbcmask uvelhom vvelhom 

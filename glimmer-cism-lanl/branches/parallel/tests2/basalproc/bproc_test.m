
%% file to create .mat input variables for Marion's coupled runs
%%
%% To be read by .py script to make .nc input file

double_res = 0;     %% flag to double res of grid
% double_res = 1;

rho_i = 910;     % ice density kg/m^3
rho_w = 1028;    % ocean water density kg/m^3

H0 = 1e3;
T0 = -10;
q0 = -7e-2;
slope = 1e-3;
b0 = 0.5;

r = 22; c = 45;
dew = 5e3; dns = 5e3;
x = [0:dew:dew*c-1]; xs = ( x(1:end-1)+x(2:end) )/2;
y = [0:dns:dns*r-1]'; ys = ( y(1:end-1)+y(2:end) )/2;


if( double_res == 1 )
    dew = 1/2*dew; dns = 1/2*dns;    % make these active to double resolution
    r = 2*r; c = 2*c;
    x2 = [0:dew:dew*c-1]; x2s = ( x2(1:end-1)+x2(2:end) )/2;
    y2 = [0:dns:dns*r-1]'; y2s = ( y2(1:end-1)+y2(2:end) )/2;
end

l = 11;     % no of vert levels


thck = H0 * ones( r, c );
artm = T0 * ones( r, c );
bheatflx = q0 * ones( r, c );
acab = b0 * ones( r, c );

topg = repmat( linspace( 1e3, 1e3-(slope*(c-1)*dew), c), r, 1 );      %% for debugging

% temp = repmat( linspace( (-H0*(rho_i/rho_w)-50)/1.25, -H0*(rho_i/rho_w)-50, c-3), r, 1 ) + 0;        %% force to be near floatationg at ds end
% topg = [ temp, repmat(temp(:,end-3), 1, 3 ) ];

usrf = topg + thck;

% %% put buffer of zero thickness cells around perimeter (for remapping)
% thck(1:3,:) = 0; thck(:,1:3) = 0; thck(end-2:end,:) = 0; thck(:,end-2:end) = 0;
% usrf(1:3,:) = 0; usrf(:,1:3) = 0; usrf(end-2:end,:) = 0; usrf(:,end-2:end) = 0;

figure(1), imagesc( x/1e3, y/1e3, thck ), axis xy, axis equal, axis tight, colorbar
xlabel( 'x (km)' ), ylabel( 'y (km)' ), title( 'thickness (m)' )

figure(2), imagesc( x/1e3, y/1e3, topg ), axis xy, axis equal, axis tight, colorbar
xlabel( 'x (km)' ), ylabel( 'y (km)' ), title( 'basal topg (m)' )

figure(3), imagesc( x/1e3, y/1e3, usrf ), axis xy, axis equal, axis tight, colorbar
xlabel( 'x (km)' ), ylabel( 'y (km)' ), title( 'upper surface (m)' )

%% load in old till map
%  load ~/Home/GLAM/GLIMGLAM/SENS/new_UPB/trunk/GLAM/Tillggl           % Marion's path
% load ~/work/modeling/glam-stream-marion-new/trunk/GLAM/Tillggl      % Steve's path

% minTauf = Tillggl;
% ind = find( minTauf < 5e3 ); minTauf(ind) = 5e3;

% minTauf = 10e3 * ones( r-1, c-1 );       % for debugging
minTauf = 1e7 * ones( r-1, c-1 );       
% minTauf(5:end-4,5:end-4) = 5e3;
minTauf(5:end-4,:) = 5e3;
if( double_res == 1 )
    minTauf(10:end-9,:) = 5e3;
end 

if( double_res == 1 )
    minTauf = interp2( xs, ys, minTauf, x2s, y2s, 'nearest' );
    ind = find( isnan( minTauf ) ); minTauf(ind) = max( max( minTauf ) );
end 

beta = 1e6*ones(size(minTauf));
ind = find( minTauf <= 5e3 ); beta(ind) = 1e2;
% beta(8:14,:) = 3e1;
% beta(10:12,:) = 0.5e1;

figure(4), imagesc( x/1e3, y/1e3, minTauf/1e3 ), axis xy, axis equal, axis tight, colorbar
xlabel( 'x (km)' ), ylabel( 'y (km)' ), title( 'Tau0 (kPa)' )

figure(5), imagesc( x/1e3, y/1e3, artm ), axis xy, axis equal, axis tight, colorbar
xlabel( 'x (km)' ), ylabel( 'y (km)' ), title( 'artm temp (C)' )

figure(6), imagesc( x/1e3, y/1e3, bheatflx ), axis xy, axis equal, axis tight, colorbar
xlabel( 'x (km)' ), ylabel( 'y (km)' ), title( 'geo flux (W m^2)' )

% %% add a kinbcmask field to specify where 0 flux bcs are
kinbcmask = zeros(size(minTauf));
% kinbcmask(1:4,:) = 1; kinbcmask(end-3:end,:) = 1; kinbcmask(:,1:4) = 1;

uvelhom = zeros( l, r-1, c-1 );
for i=1:l
    uvelhom(i,:,:) = zeros(size(minTauf));
end

vvelhom = uvelhom; 

figure(7), imagesc( x/1e3, y/1e3, kinbcmask ), axis xy, axis equal, axis tight, colorbar
xlabel( 'x (km)' ), ylabel( 'y (km)' ), title( 'kinbcmask' )

figure(8), imagesc( x/1e3, y/1e3, acab ), axis xy, axis equal, axis tight, colorbar
xlabel( 'x (km)' ), ylabel( 'y (km)' ), title( 'acab (m/a)' )

tauf = minTauf;

%  cd ~/Home/Glimmer2/glimmer-cism-lanl/branches/basalproc/tests/basalproc     % Marion's path
cd /Users/sprice/work/modeling/cism-parallel/tests/basalproc            % Steve's path

save bproc.mat artm acab bheatflx usrf topg thck beta tauf kinbcmask uvelhom vvelhom 

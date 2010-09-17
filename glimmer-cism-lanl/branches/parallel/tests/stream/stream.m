
%% PLOT MODEL STREAM/PLASTIC BED RESULTS AND COMPARE W/ RAYMOND/SCHOOF
%% SOLUTIONS

clear all

cd /Users/sprice/work/modeling/cism-parallel/tests/stream/

rho = 910;
g = -9.81;

n = 3;
m = 2;

r = 20;         % !! r needs to be even # divisible by 2 !!
c = 30;
levels = 5;

L = 1.5e4; 
A = 1e-16;
B = A^(-1/n);

H = 1000;
dsdx = -1.0e-3;

% ice stream width
w = 50e3;
W = w / 2;

yy = linspace( 0, W, (r-4)/2 );

dy = mean( diff( yy ) );
dx = dy;

taud = rho * g * H * dsdx;

% Schoof yield stress distribution
tau0s = taud * abs( yy / L ).^m;

% yield stress (for Raymond solution only)
% tau0r = mean( tau0s );
tau0r = 5e3;

% Raymond solution
ur = 2 * A / (n+1) * ( (taud - tau0r)/H )^n * ( W^(n+1) - yy.^(n+1) );


us = -2*taud^3*L^4/(B^3*H^3) * ( ((yy/L).^4 - (m+1)^(4/m))/4 - 3*( abs(yy/L).^(m+4) ...
    - (m+1)^(1+4/m) )/((m+1)*(m+4)) + 3*( abs(yy/L).^(2*m+4) - (m-1)^(2+4/m) )/((m+1)^2*(2*m+4)) ...
    - ( abs(yy/L).^(3*m+4) - (m+1)^(3+4/m) )/ ( (m+1)^3*(3*m+4)) );

sscale = 2*taud^3*L^4/(B^3*H^3);

ind = find( abs( yy ) >= W ); us(ind) = min( min( us ) );

us = us - min( us );


figure(198), clf
subplot(2,1,1), hold on
plot( yy/1e3, ur - min(ur), 'k-', 'linewidth', 1.0 ), hold on, grid on
xlabel( 'dist across flow (m)'), ylabel( 'velocity (m/a)'), title( 'Raymond solution' )
subplot(2,1,2), hold on
plot( yy/1e3, tau0r*ones(size(yy))/1e3, 'k-', 'linewidth', 1.0 ), hold on, grid on
xlabel( 'dist across flow (m)'), ylabel( 'yield stress (kPa)')

figure(199), clf
subplot(2,1,1), hold on
plot( yy/1e3, us - min(us), 'k-', 'linewidth', 1.0 ), hold on, grid on
xlabel( 'dist across flow (m)'), ylabel( 'velocity (m/a)'), title( 'Schoof solution' )
subplot(2,1,2), hold on
plot( yy/1e3, tau0s/1e3, 'k-', 'linewidth', 1.0 ), hold on, grid on
xlabel( 'dist across flow (m)'), ylabel( 'yield stress (kPa)')


%% put data into .mat file, to be read by .py script w/ config file, to make input .nc files
thck = H * ones( r, c );
topg = repmat( 5e3 + dsdx*dx*[0:c-1], r, 1 );
usrf = topg + thck;


tauf_profile = [ fliplr(tau0s(2:end)), tau0s ]';        % assign Schoof profile
tauf = 1e6 * ones( r-1, c-1 );
tauf(3:end-2,3:end-2) = repmat( tauf_profile, 1, c-5 );     %% no slip at/near up/downstream ends
% tauf(3:end-2,:) = repmat( tauf_profile, 1, c-1 );     %% use periodic bcs for cont. along flow

% tauf_profile = tau0r*ones(r-7,1);        % assign Raymond profile
% tauf = 1e10 * ones( r-1, c-1 );
% % tauf(4:end-3,3:end-2) = repmat( tauf_profile, 1, c-5 );     %% no slip at/near up/downstream ends
% tauf(4:end-3,:) = repmat( tauf_profile, 1, c-1 );     %% use periodic bcs for cont. along flow


% figure(200),clf
% subplot(2,2,1),hold on
% imagesc( thck ), axis xy, axis equal, axis tight, colorbar, title( 'thickness (m)' )
% subplot(2,2,2),hold on
% imagesc( usrf ), axis xy, axis equal, axis tight, colorbar, title( 'surface (m)' )
% subplot(2,2,3),hold on
% imagesc( topg ), axis xy, axis equal, axis tight, colorbar, title( 'bed (m)' )
% subplot(2,2,4),hold on
% imagesc( tauf/1e3 ), axis xy, axis equal, axis tight, colorbar, title( 'yield stress (kPa)' )



% %% optional: add analytic solution at up/downstream ends as kin vel bc
kinbcmask = zeros(size(tauf));
uvelhom = zeros( levels, r-1, c-1 );
vvelhom = zeros( levels, r-1, c-1 );

% comment these out for now kin bc
ur_profile = repmat( [ 0 0 fliplr(ur) ur(2:end) 0 0 ], levels, 1 );     
us_profile = repmat( [ 0 0 fliplr(us) us(2:end) 0 0 ], levels, 1 );
% uvelhom( :, :, end ) = ur_profile; uvelhom( :, :, 1 ) = ur_profile; 
uvelhom( :, :, end ) = us_profile; uvelhom( :, :, 1 ) = us_profile; 
kinbcmask(:,1) = 1; kinbcmask(:,end) = 1;
ind = find( tauf <= 0 ); tauf(ind) = 1e-10;



save stream.mat usrf topg thck tauf uvelhom vvelhom kinbcmask levels 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot output


%% open file
filename = 'stream.out.nc'; 
ncid = netcdf.open( filename, 'nowrite' );

%% get id for variable names (use ncview to find var names?)
id_uvel = netcdf.inqvarid( ncid, 'uvelhom' );
id_vvel = netcdf.inqvarid( ncid, 'vvelhom' );
id_btractx = netcdf.inqvarid( ncid, 'btractx' );
id_btracty = netcdf.inqvarid( ncid, 'btracty' );

uvel = permute( netcdf.getvar(ncid, id_uvel ), [ 2 1 3 ] );
vvel = permute( netcdf.getvar(ncid, id_vvel ), [ 2 1 3 ] );
btractx = netcdf.getvar(ncid, id_btractx )';
btracty = netcdf.getvar(ncid, id_btracty )';

yy2 = [ yy yy(end)+dx yy(end)+2*dx ];
yy2 = [ -fliplr(yy2(2:end)), yy2 ];

figure(198)
subplot(2,1,1), hold on
plot( yy2/1e3, uvel(:,round(c/2),1), 'bo:' )
plot( yy2/1e3, uvel(:,end,1), 'b*' )              %% boundary value
legend( 'analytic', 'model', 'model boundary' )
subplot(2,1,2), hold on
plot( yy2/1e3, -btractx(:,round(c/2),1)/1e3, 'bo:' )


figure(199)
subplot(2,1,1), hold on
plot( yy2/1e3, uvel(:,round(c/2),1), 'bo:' )
plot( yy2/1e3, uvel(:,end,1), 'b*' )              %% boundary value
legend( 'analytic', 'model', 'model boundary' )
subplot(2,1,2), hold on
plot( yy2/1e3, -btractx(:,round(c/2),1)/1e3, 'bo:' )


%% Plot flat shelf profiles

f1 = '/home/gladish/research/gcp_resources/nc_files/last_time.flat_shelf_100_m_sept15.2.out.nc';
f2 = '/home/gladish/research/gcp_resources/nc_files/last_time.flat_shelf_200_m_sept15.2.out.nc';
f3 = '/home/gladish/research/gcp_resources/nc_files/last_time.flat_shelf_500_m_sept15.2.out.nc';

m = 5;
n = 100;

hy = 40000.0 / n;
hx = hy;

kinbcw = 2;
x0 = (n-kinbcw-0.5)*hy;
x1 = (4.0)*hy;
A = 1.0*10^(-16);   % Pa^(-3)* a^(-1)
n_points = 1000;

u0 = -1000.0;

rhoi = 910.0;
rhoo = 1028.0;
g = 9.81;

fs = 18;

hold on

[~,y0,~,y1,thck,~,vvel] = nc_read(f1, -1);
[ycvel,ycthk, vvelc, thkc] = centerline_profile(y0,y1,vvel,thck);
h0 = 100.0;
[y,v] = flat_shelf_ice(A,rhoi,rhoo,g,x0,x1,n_points, u0, h0);
plot(ycvel/1000, vvelc, 'b*');
plot(y/1000, v, 'k');

[~,y0,~,y1,thck,~,vvel] = nc_read(f2, -1);
[ycvel,ycthk, vvelc, thkc] = centerline_profile(y0,y1,vvel,thck);
h0 = 200.0;
[y,v] = flat_shelf_ice(A,rhoi,rhoo,g,x0,x1,n_points, u0, h0);
plot(ycvel/1000, vvelc, 'r*');
plot(y/1000, v, 'k');

[~,y0,~,y1,thck,~,vvel] = nc_read(f3, -1);
[ycvel,ycthk, vvelc, thkc] = centerline_profile(y0,y1,vvel,thck);
h0 = 500.0;
[y,v] = flat_shelf_ice(A,rhoi,rhoo,g,x0,x1,n_points, u0, h0);
plot(ycvel/1000, vvelc, 'g*');
plot(y/1000, v, 'k');

title('Diagnosed Velocity - flat shelf','FontSize',fs);
ylabel('m/year','FontSize',fs);
xlabel('km','FontSize',fs);

legend( 'model vel (100 m thk)','exact vel (100 m thk)', ...
        'model vel (200 m thk)','exact vel (200 m thk)', ...
        'model vel (500 m thk)','exact vel (500 m thk)', ...
       'Location','BestOutside');    

hold off


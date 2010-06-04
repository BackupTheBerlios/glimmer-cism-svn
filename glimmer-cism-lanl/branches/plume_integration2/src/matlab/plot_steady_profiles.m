%% Plot steady profiles

f1 = 'iso1.nc';
f2 = 'iso2.nc';
f4 = 'iso4.nc';

ystart = 40*1000-1000;
yend = 3500;
A = 4.6*10^(-18) * (365.242/365.25)^(-1);
n = 1000;
u0 = -1000;
h0 = 1000;
rhoi = 910.0;
rhoo = 1028.0;
g = 9.81;
k = 0.9875;

rhoo = rhoo*k;

[~,y0,~,y1,thck,~,vvel] = nc_read(f1, 1390);
[ycvel,ycthk, vvelc, thkc] = centerline_profile(y0,y1,vvel,thck);

bmelt = 0.0;
[y,v,h,w] = steady_ice_1(k*A,rhoi,rhoo,g,ystart,yend,n, u0, h0, bmelt);
subplot(2,1,1);
hold on
plot(ycvel/1000, vvelc, 'b*');
plot(y/1000, v, 'b');
title('Along flow velocity');
ylabel('meters/year');
xlabel('km');
subplot(2,1,2);
hold on
plot(ycthk/1000, thkc, 'b*');
plot(y/1000, h, 'b');
title('Along flow thickness');
ylabel('meters');
xlabel('km');

[~,y0,~,y1,thck,~,vvel] = nc_read(f2, 1498);
[ycvel,ycthk, vvelc, thkc] = centerline_profile(y0,y1,vvel,thck);
bmelt = -10.0;
[y,v,h,w] = steady_ice_1(k*A,rhoi,rhoo,g,ystart,yend,n, u0, h0, bmelt);
subplot(2,1,1);
plot(ycvel/1000, vvelc, 'r*');
plot(y/1000, v, 'r');
subplot(2,1,2);
plot(ycthk/1000, thkc, 'r*');
plot(y/1000, h, 'r');

[~,y0,~,y1,thck,~,vvel] = nc_read(f4, 1738);
[ycvel,ycthk, vvelc, thkc] = centerline_profile(y0,y1,vvel,thck);
bmelt = -25.0;
[y,v,h,w] = steady_ice_1(k*A,rhoi,rhoo,g,ystart,yend,n, u0, h0, bmelt);
subplot(2,1,1);
plot(ycvel/1000, vvelc, 'g*');
plot(y/1000, v, 'g');
legend('0.0 m/year', '0.0 m/year Analytic Solution', ...
       '10.0 m/year','10.0 m/year Analytic Solution', ...
       '25.0 m/year','25.0 m/year Analytic Solution', ...
       'Location','NorthWest');
hold off
subplot(2,1,2);
plot(ycthk/1000, thkc, 'g*');
plot(y/1000, h, 'g');
legend('0.0 m/year', '0.0 m/year Analytic Solution', ...
       '10.0 m/year','10.0 m/year Analytic Solution', ...
       '25.0 m/year','25.0 m/year Analytic Solution', ...
       'Location','NorthWest');
hold off

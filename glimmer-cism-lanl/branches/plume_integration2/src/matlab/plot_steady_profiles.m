%% Plot steady profiles

f1 = '/archive/cvg222/gc_output/1d_isoT_fixedA_5/1d_isoT_fixedA.out.5.nc';
f2 = '/archive/cvg222/gc_output/1d_isoT_fixedA_6/1d_isoT_fixedA.out.4.nc';
f4 = '/archive/cvg222/gc_output/1d_isoT_fixedA_7/1d_isoT_fixedA.out.5.nc';

ystart = 40*1000-1000;
yend = 3500;
A = 1.0*10^(-16) * (365.242/365.25)^(-1);
n = 1000;
u0 = -1000;
h0 = 1000;
rhoi = 910.0;
rhoo = 1028.0;
g = 9.81;

k = 0.96;
%k = 1.0;
rhoo = rhoo*k;

fs = 18;

figure;

[~,y0,~,y1,thck,~,vvel] = nc_read(f1, -1);
[ycvel,ycthk, vvelc, thkc] = centerline_profile(y0,y1,vvel,thck);
bmelt = 0.0;
[y,v,h,~] = steady_ice_1(A,rhoi,rhoo,g,ystart,yend,n, u0, h0, bmelt);
subplot(2,1,1);
hold on
plot(ycvel/1000, vvelc, 'b*');
plot(y/1000, v, 'b');
title('Velocity','FontSize',fs);
ylabel('m/year','FontSize',fs);
xlabel('km','FontSize',fs);
subplot(2,1,2);
hold on
plot(ycthk/1000, thkc, 'b*');
plot(y/1000, h, 'b');
title('Thickness','FontSize',fs);
ylabel('meters','FontSize',fs);
xlabel('km','FontSize',fs);

[~,y0,~,y1,thck,~,vvel] = nc_read(f2, -1);
[ycvel,ycthk, vvelc, thkc] = centerline_profile(y0,y1,vvel,thck);
bmelt = -10.0;
[y,v,h,~] = steady_ice_1(A,rhoi,rhoo,g,ystart,yend,n, u0, h0, bmelt);
subplot(2,1,1);
plot(ycvel/1000, vvelc, 'r*');
plot(y/1000, v, 'r');
subplot(2,1,2);
plot(ycthk/1000, thkc, 'r*');
plot(y/1000, h, 'r');

[~,y0,~,y1,thck,~,vvel] = nc_read(f4, -1);
[ycvel,ycthk, vvelc, thkc] = centerline_profile(y0,y1,vvel,thck);
bmelt = -25.0;
[y,v,h,w] = steady_ice_1(A,rhoi,rhoo,g,ystart,yend,n, u0, h0, bmelt);
subplot(2,1,1);
plot(ycvel/1000, vvelc, 'g*');
plot(y/1000, v, 'g');
legend('0.0 m/year', '0.0 m/year Exact Sol.', ...
       '10.0 m/year','10.0 m/year Exact Sol.', ...
       '25.0 m/year','25.0 m/year Exact Sol.', ...
       'Location','BestOutside') %,'FontSize',fs);
hold off
subplot(2,1,2);
plot(ycthk/1000, thkc, 'g*');
plot(y/1000, h, 'g');
legend('0.0 m/year', '0.0 m/year Exact Sol.', ...
       '10.0 m/year','10.0 m/year Exact Sol.', ...
       '25.0 m/year','25.0 m/year Exact Sol.', ...
       'Location','BestOutside') %,'FontSize',fs);
hold off

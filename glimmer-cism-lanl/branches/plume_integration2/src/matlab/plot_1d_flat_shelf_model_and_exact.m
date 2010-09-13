%% Plot flat shelf profiles

f1 = '/scratch/gc_output/l9/sunstudio.out.nc';
%f2 = '/scratch/gc_output/l6/sunstudio.out.nc';

hx = 1000;
hy = hx;
m = 5;
n = 41;
kinbcw = 2;
x0 = (n-kinbcw-0.5)*hy;
x1 = (4-0.5)*hy;
A = 1.0*10^(-16);   % Pa^(-3)* a^(-1)
n_points = 1000;
u0 = -1000;
h0 = 1000;
rhoi = 910.0;
rhoo = 1028.0;
g = 9.81;

fs = 18;

[~,y0,~,y1,thck,~,vvel] = nc_read(f1, -1);
[ycvel,ycthk, vvelc, thkc] = centerline_profile(y0,y1,vvel,thck);
[y,v] = flat_shelf_ice(A,rhoi,rhoo,g,x0,x1,n_points, u0, h0);

%subplot(2,1,1);
hold on
plot(ycvel/1000, vvelc, 'b*');
plot(y/1000, v, 'r');
title('Velocity','FontSize',fs);
ylabel('m/year','FontSize',fs);
xlabel('km','FontSize',fs);


%[~,y0,~,y1,thck,~,vvel] = nc_read(f2, -1);
%[ycvel,ycthk, vvelc, thkc] = centerline_profile(y0,y1,vvel,thck);
%plot(ycvel/1000, vvelc, 'r*');

%subplot(2,1,2);
%hold on
%plot(ycthk/1000, thkc, 'b*');
%plot(y/1000, h, 'b');
%title('Thickness','FontSize',fs);
%ylabel('meters','FontSize',fs);
%xlabel('km','FontSize',fs);

legend('computed v vel (m/year)', 'exact sol. v vel (m/year)', ...
       'Location','BestOutside');    
%       '10.0 m/year','10.0 m/year Exact Sol.', ...
%       '25.0 m/year','25.0 m/year Exact Sol.', ...

hold off


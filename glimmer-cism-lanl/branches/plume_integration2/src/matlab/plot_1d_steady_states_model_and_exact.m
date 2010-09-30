%% Plot steady profiles

f0 = '~/research/gcp_resources/nc_files/one_d_0_bmlt_15.4.out.nc';
f1 = '~/research/gcp_resources/nc_files/one_d_10_bmlt_15.4.out.nc';
f2 = '~/research/gcp_resources/nc_files/one_d_25_bmlt_15.4.out.nc';

%f0 = '~/research/gcp_resources/nc_files/1d_0_bmlt_sept15.3.out.nc';
%f1 = '~/research/gcp_resources/nc_files/1d_10_bmlt_sept15.3.out.nc';
%f2 = '~/research/gcp_resources/nc_files/1d_25_bmlt_sept15.3.out.nc';

f0 = '~/research/gcp_resources/nc_files/one_d_0_bmlt_sept15.40.out.nc';
f1 = '~/research/gcp_resources/nc_files/one_d_10_bmlt_sept15.40.out.nc';
f2 = '~/research/gcp_resources/nc_files/one_d_25_bmlt_sept15.40.out.nc';

m = 5;
n = 80;
%n= 100;
n = 40;

hx = 40000.0 / n;
hy = hx;

kinbcw = 2;
ystart = (n-kinbcw-0.5)*hy;
yend = (4)*hy;

A = 1.0*10^(-16); 

n_points = 1000;
u0 = -1000;
h0 = 1000;

rhoi = 910.0;
rhoo = 1028.0;
g = 9.81;

fs = 18;
lwidth = 1;
lwidth_dot = 1;

plot0 = true;
plot1 = true;
plot2 = true;
plot5 = false;

[~,y0,~,y1,thck,~,vvel] = nc_read(f0, -1);
[ycvel,ycthk, vvelc, thkc] = centerline_profile(y0,y1,vvel,thck);
acab = 0.0;
[y,v,h,~] = steady_ice_1(A,rhoi,rhoo,g,ystart,yend,n_points, u0, h0, acab);

subplot(2,1,1);
hold on
plot(ycvel/1000, vvelc, 'b*','LineWidth',lwidth_dot);
plot((y-0.5*hy)/1000, v, 'b','LineWidth',lwidth);
title('Steady Ice Velocity','FontSize',fs);
ylabel('m/year','FontSize',fs);
xlabel('km','FontSize',fs);

subplot(2,1,2);
hold on
plot(ycthk/1000, thkc, 'b*','LineWidth',lwidth_dot);
plot((y)/1000, h, 'b','LineWidth',lwidth);
title('Steady Ice Thickness','FontSize',fs);
ylabel('meters','FontSize',fs);
xlabel('km','FontSize',fs);

if (plot1) 
[~,y0,~,y1,thck,~,vvel] = nc_read(f1, -1);
[ycvel,ycthk, vvelc, thkc] = centerline_profile(y0,y1,vvel,thck);
acab = -10.0;
[y,v,h,~] = steady_ice_1(A,rhoi,rhoo,g,ystart,yend,n_points, u0, h0, acab);

subplot(2,1,1);
plot(ycvel/1000, vvelc, 'k*','LineWidth',lwidth_dot);
plot((y-0.5*hy)/1000, v, 'k','LineWidth',lwidth);

subplot(2,1,2);
plot(ycthk/1000, thkc, 'k*','LineWidth',lwidth_dot);
plot(y/1000, h, 'k','LineWidth',lwidth);
end

if (plot2) 
[~,y0,~,y1,thck,~,vvel] = nc_read(f2, -1);
[ycvel,ycthk, vvelc, thkc] = centerline_profile(y0,y1,vvel,thck);
acab = -25.0;
[y,v,h,~] = steady_ice_1(A,rhoi,rhoo,g,ystart,yend,n_points, u0, h0, acab);

subplot(2,1,1);
plot(ycvel/1000, vvelc, 'r*','LineWidth',lwidth_dot);
plot((y-0.5*hy)/1000, v, 'r','LineWidth',lwidth);

subplot(2,1,2);
plot(ycthk/1000, thkc, 'r*','LineWidth',lwidth_dot);
plot(y/1000, h, 'r','LineWidth',lwidth);
end

if (plot5)
[~,y0,~,y1,thck,~,vvel] = nc_read(f5, -1);
[ycvel,ycthk, vvelc, thkc] = centerline_profile(y0,y1,vvel,thck);
acab = -5.0;
[y,v,h,w] = steady_ice_1(A,rhoi,rhoo,g,ystart,yend,n_points, u0, h0, acab);

subplot(2,1,1);
plot(ycvel/1000, vvelc, 'g*','LineWidth',lwidth_dot);
plot((y-0.5*hy)/1000, v, 'g','LineWidth',lwidth);

subplot(2,1,2);
plot(ycthk/1000, thkc, 'g*','LineWidth',lwidth_dot);
plot(y/1000, h, 'g','LineWidth',lwidth);

end

if (plot5)  
subplot(2,1,1);
legend('0.0 m/year', '0.0 m/year Exact Sol.', ...
       '5.0 m/year', '5.0 m/year Exact Sol.', ...
       '10.0 m/year','10.0 m/year Exact Sol.', ...
       '25.0 m/year','25.0 m/year Exact Sol.', ...
       'Location','BestOutside') %,'FontSize',fs);
hold off;
subplot(2,1,2);
legend('0.0 m/year', '0.0 m/year Exact Sol.', ...
       '5.0 m/year', '5.0 m/year Exact Sol.', ...
       '10.0 m/year','10.0 m/year Exact Sol.', ...
       '25.0 m/year','25.0 m/year Exact Sol.', ...
       'Location','BestOutside') %,'FontSize',fs);
hold off;
else
    
subplot(2,1,1);
legend('0.0 m/year', '0.0 m/year Exact Sol.', ...
       '10.0 m/year','10.0 m/year Exact Sol.', ...
       '25.0 m/year','25.0 m/year Exact Sol.', ...
       'Location','BestOutside') %,'FontSize',fs);
hold off
subplot(2,1,2);
legend('0.0 m/year', '0.0 m/year Exact Sol.', ...
       '10.0 m/year','10.0 m/year Exact Sol.', ...
       '25.0 m/year','25.0 m/year Exact Sol.', ...
       'Location','BestOutside') %,'FontSize',fs);
hold off

end
   
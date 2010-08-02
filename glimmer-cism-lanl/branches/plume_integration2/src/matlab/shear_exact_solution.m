modeloutput = '~/computation/gc/testruns/shear_stress_test2/sunstudio_8.out.nc';
[~,y0,~,y1,thck,~,vvel] = nc_read(modeloutput, -1);
[yvel,ythk, velc, thkc] = centerline_profile(y0,y1,vvel,thck);

A = 10.0^(-16);
rhoi = 910.0;
rhoo = 1028.0;
g = 9.81;
Vin = -1000.0;
Hin = 600.0;

n = length(y0);
S = double(y0(n-1));
yout = double(y0(4)) ;
yin = S;

ny = 50;
dy = (yin - yout) / (ny-1);

a = 0.0;

tauxy0 = 10.0*10^3;

L = 5*1000.0;

k = (1-rhoi/rhoo)*(rhoi*g/4.0);
q = tauxy0/(2*L);

fs = 20;

f = @(y) Vin*Hin - a*(yin - y);

% u = [v; w]

odefun = @(y,u) [ A*(k*f(y)/u(1) - q*u(2)*u(1)/f(y))^3.0;
                  f(y)/u(1) ];

bcfun = @(u0, u1) [u0(2); u1(1)-Vin];

solinit.x = yout:dy:yin;

solinit.y = zeros(2,ny);
solinit.y(1,:) = Vin;
solinit.y(2,:) = 0:(Hin*S*0.5/(ny-1)):Hin*S*0.5;

sol = bvp4c(odefun, bcfun,solinit);

h = f(sol.x) ./ sol.y(1,:);

h_Pat = @(y) Hin + 2*tauxy0/(rhoi*g*L*(1-rhoi/rhoo))*(y-yin);

figure;

subplot(2,1,1);
hold on
plot(sol.x / 1000.0, sol.y(1,:),'r*');
%p = polyfit(sol.x, sol.y(1,:),1);
%plot(sol.x / 1000.0, p(1)*sol.x + p(2) ,'g*');
plot(yvel / 1000.0, velc, 'k*');

xlabel('y (km)','FontSize',fs);
ylabel('velocity','FontSize',fs);
set(gca,'FontSize',fs);
title('Centerline velocities (tauxy0 = 10 kPa)','FontSize',fs);
legend('quasi-analytic', ...
       'calculated','Location','BestOutside') ;

hold off

subplot(2,1,2);
hold on
plot(sol.x / 1000.0, h, 'b*');
%p = polyfit(sol.x, h, 1);
%plot(sol.x / 1000.0, p(1)*sol.x + p(2),'g*');
plot(ythk / 1000.0, thkc, 'k*');
plot(ythk / 1000.0, h_Pat(ythk), 'g*');
xlabel('y (km)','FontSize',fs);
ylabel('height','FontSize',fs);
title('Centerline thickness (tauxy0 = 10 kPa)','FontSize',fs);
set(gca,'FontSize',fs);
legend('quasi-analytic','calculated','Paterson','Location','BestOutside') ;
hold off

%subplot(3,1,3);
%hold on
%plot(sol.x / 1000.0, sol.y(2,:) ./ h, 'b*');
%ylabel('w / h', 'FontSize',fs);
%xlabel('y(km)', 'FontSize',fs);
%set(gca, 'FontSize',fs);
%hold off


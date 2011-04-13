function [] = plot_velocity_overlay(data)

fs = 9;
figure('Units','centimeters','Position',[0 0 40 20]);
x = data.x;
y = data.y;
su = data.su;
sv = data.sv;
u = data.u;
v = data.v;

ref_speed = 0.1;
stride = 2;
su(1+stride,1+stride) = 0.0;
sv(1+stride,1+stride) = ref_speed;
scale = 2.5;

subplot(1,2,1);

hold on
contourf(x/1000.0,y/1000.0,data.draft',20,'EdgeColor','None');colorbar;
quiver(x(1:stride:end)/1000,y(1:stride:end)/1000,...
      su(1:stride:end,1:stride:end)',...
      sv(1:stride:end,1:stride:end)',...
     scale,'w');
%colormap hsv
hold off
%quiverscale(max(max(su)),max(max(sv)),gca);

%contourf(x/1000.0,y/1000.0,sqrt(su.*su + sv.*sv)',20);colorbar; caxis([0 0.05])
  title(strcat(['plume velocity (ref speed is ',sprintf('%0.2f',ref_speed),' m/s',...
					       ') with draft contours']),...
        'FontSize',fs);

subplot(1,2,2);

speed = sqrt(su.*su + sv.*sv);
hold on
contourf(data.x/1000.0,data.y/1000.0,speed',20,'EdgeColor','None') ;colorbar;
contour(data.x/1000.0,data.y/1000.0,data.draft',20,'w');
caxis([min(min(speed)), max(max(speed))]);
title('plume speed (m/s) with draft contours','FontSize',fs);

end

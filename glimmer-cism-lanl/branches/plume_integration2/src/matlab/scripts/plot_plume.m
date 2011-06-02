function [] = plot_plume(data )

figure(1);
fs = 12;

x = data.x;
y = data.y;
su = data.su;
sv = data.sv;
u = data.u;
v = data.v;
train = data.train;
bmelt = data.bmelt;
temp = data.temp;

stride = 1;
scale = 2.0;

%subplot(2,2,1);
figure(1);
clf;
hold on
contourf(x/1000.0,y/1000.0,-data.draft',40,'EdgeColor','None');colorbar;
caxis([0 800])
quiver(x(1:stride:end)/1000,y(1:stride:end)/1000,...
      su(1:stride:end,1:stride:end)',...
      sv(1:stride:end,1:stride:end)',...
     scale,'k');
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Ocean Velocities with contours of Ice Draft (m)','FontSize',fs);
hold off
print('-depsc','/home/gladish/research/visuals/plume_trial_8_plume_vel');
%subplot(2,2,2);
%contourf(x/1000.0,y/1000.0,sqrt(su.*su + sv.*sv)',20);colorbar; caxis([0 0.05])
%  title(strcat(['plume velocity (ref speed is ',sprintf('%0.2f',ref_speed),' m/s',...
%					       ') with draft contours']),...
%        'FontSize',fs);
figure(1);
clf;
speed = sqrt(su.*su + sv.*sv);
hold on
contourf(data.x/1000.0,data.y/1000.0,data.bmelt',20,'EdgeColor','None') ;colorbar;
%contour(data.x/1000.0,data.y/1000.0,data.draft',30,'w');
caxis([min(min(bmelt)), max(max(bmelt))]);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Basal Melt Rate (m/year)','FontSize',fs);

hold off
print('-depsc','/home/gladish/research/visuals/plume_trial_8_meltrate');


figure(1);
clf;
hold on
contourf(data.x/1000.0,data.y/1000.0,train',...
         30,'EdgeColor','None') ;colorbar;
contour(data.x/1000.0,data.y/1000.0,train',[0 0],'k')
%contour(data.x/1000.0,data.y/1000.0,data.draft',20,'w');
%contour(data.x/1000.0,data.y/1000.0,data.draft',30,'w');
caxis([min(min(train)), max(max(train))]);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Entrainment/Detrainment Rate (m/year)','FontSize',fs);
hold off
print('-depsc','/home/gladish/research/visuals/plume_trial_8_train');

figure(1);
clf;
hold on
contourf(data.x/1000.0,data.y/1000.0,temp',...
         30,'EdgeColor','None') ;colorbar;
%contour(data.x/1000.0,data.y/1000.0,data.draft',20,'w');
%contour(data.x/1000.0,data.y/1000.0,data.draft',30,'w');
caxis([min(min(temp)), max(max(temp))]);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Ocean Mixed Layer Temperature (C)','FontSize',fs);
hold off
print('-depsc','/home/gladish/research/visuals/plume_trial_8_temp');

figure(1);
clf;

[m,n] = size(data.draft);
channel_amp = zeros(m,n);
mean_draft = zeros(m,n);
max_draft = zeros(m,n);

for i=1:n
  channel_amp(:,i) = data.draft(:,i) - min(data.draft(:,i));
  mean_draft(:,i) = sum(data.draft(:,i))/m;
  max_draft(:,i) = min(data.draft(:,i));
end

%subplot(1,1,1);
  hold on
  contourf(x/1000.0,y/1000.0,channel_amp',20,'EdgeColor','None') ;colorbar;
%  contour(x/1000.0,y/1000.0,mean_draft',20,'w')
  [C,h] = contour(x/1000.0,y/1000.0,max_draft',20,'w');
  caxis([min(min(channel_amp)) max(max(channel_amp))]);
  xlabel('Across shelf distance (km)','FontSize',fs);
  ylabel('Along shelf distance (km)','FontSize',fs);
  clabel(C,'FontSize',fs,'Color','k')
  title('Channel depth with deepest-draft contours','FontSize',fs);

  print('-depsc','/home/gladish/research/visuals/plume_trial_8_channels');  
end


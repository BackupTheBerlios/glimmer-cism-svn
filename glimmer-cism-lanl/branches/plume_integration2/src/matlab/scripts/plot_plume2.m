function [] = plot_plume2(data,dice,timeslice,flabel)

figure(1);
fs = 14;

if (timeslice < 0) 
    timeslice = size(data.su, 3);
end 

x = data.x;
y = data.y;
su = squeeze(data.su(:,:,timeslice));
sv = squeeze(data.sv(:,:,timeslice));
u = squeeze(data.u(:,:,timeslice));
v = squeeze(data.v(:,:,timeslice));
train = squeeze(data.train(:,:,timeslice));
bmelt = squeeze(data.bmelt(:,:,timeslice));
temp = squeeze(data.temp(:,:,timeslice));
draft = squeeze(data.draft(:,:,timeslice));
thk_t = squeeze(dice.thk_t(:,:,timeslice));

grad = squeeze(data.grad(:,:,timeslice));
stride = 1;
scale = 3.0;
lw = 1.0;

figure(1);
clf;


hold on
contourf(x/1000.0,y/1000.0,grad',40,'EdgeColor','None');colorbar('FontSize',fs);
%caxis([0 0.2])
quiver(x(1:stride:end)/1000,y(1:stride:end)/1000,...
      su(1:stride:end,1:stride:end)',...
      sv(1:stride:end,1:stride:end)',...
     scale,'w','LineWidth',lw);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Ocean Velocities with contourfs of ice draft gradient','FontSize',fs);
hold off
fname = strcat(['/home/gladish/research/visuals/plume_',flabel,'_grad']);
print('-depsc',fname);

%subplot(2,2,1);
figure(1);
clf;
hold on
contourf(x/1000.0,y/1000.0,-draft',40,'EdgeColor','None');colorbar('FontSize',fs);
caxis([000 800])
%quiver(x(1:stride:end)/1000,y(1:stride:end)/1000,...
%      su(1:stride:end,1:stride:end)',...
%      sv(1:stride:end,1:stride:end)',...
%     scale,'k','LineWidth',lw);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Ocean Velocities with contours of Ice Draft (m)','FontSize',fs);
title('Contours of Ice Draft (m)','FontSize',fs);
hold off
print('-depsc',strcat(['/home/gladish/research/visuals/plume_',flabel,'_plume_vel']));

%subplot(2,2,4);
%contourf(x/1000.0,y/1000.0,sqrt(su.*su + sv.*sv)',20);colorbar('FontSize',fs); caxis([0 0.05])
%  title(strcat(['plume velocity (ref speed is ',sprintf('%0.2f',ref_speed),' m/s',...
%					       ') with draft contours']),...
%        'FontSize',fs);

figure(1);
clf;
speed = sqrt(su.*su + sv.*sv);
hold on
contourf(x/1000.0,y/1000.0,bmelt',20,'EdgeColor','None') ;colorbar('FontSize',fs);
%contour(x/1000.0,y/1000.0,draft',30,'w');
caxis([min(min(bmelt)), max(max(bmelt))]);
caxis([0 40]);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Basal Melt Rate (m/year)','FontSize',fs);
hold off
print('-depsc',strcat(['/home/gladish/research/visuals/plume_',flabel,'_meltrate']));

figure(1);
clf;
hold on
%subplot(2,2,3)
contourf(x/1000.0,y/1000.0,train'/1000.0,...
         30,'EdgeColor','None') ;colorbar('FontSize',fs);
contour(x/1000.0,y/1000.0,train'/1000.0,[0 0],'k')
%contour(x/1000.0,y/1000.0,draft',20,'w');
%contour(x/1000.0,y/1000.0,draft',30,'w');
caxis([min(min(train/1000.0)), max(max(train/1000.0))]);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Entrainment/Detrainment Rate (1000m/year)','FontSize',fs);
hold off
print('-depsc',strcat(['/home/gladish/research/visuals/plume_',flabel,'_train']));

figure(1);
clf;
%subplot(2,2,3);
hold on
contourf(x/1000.0,y/1000.0,temp',...
         30,'EdgeColor','None') ;colorbar('FontSize',fs);
%contour(data.x/1000.0,data.y/1000.0,data.draft',20,'w');
%contour(data.x/1000.0,data.y/1000.0,data.draft',30,'w');
caxis([min(min(temp)), max(max(temp))]);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Ocean Mixed Layer Temperature (C)','FontSize',fs);
hold off
print('-depsc',strcat(['/home/gladish/research/visuals/plume_',flabel,'_temp']));

figure(1);
clf;

[m,n] = size(draft);
channel_amp = zeros(m,n);
mean_draft = zeros(m,n);
max_draft = zeros(m,n);

for i=1:n
  channel_amp(:,i) = draft(:,i) - min(draft(:,i));
  mean_draft(:,i) = sum(draft(:,i))/m;
  max_draft(:,i) = min(draft(:,i));
end

%subplot(2,2,2);
  hold on
  contourf(x/1000.0,y/1000.0,channel_amp',20,'EdgeColor','None') ;colorbar('FontSize',fs);
%  contour(x/1000.0,y/1000.0,mean_draft',20,'w')
  [C,h] = contour(x/1000.0,y/1000.0,max_draft',10,'w');
  caxis([min(min(channel_amp)) max(max(channel_amp))]);
  caxis([0 300]);
  xlabel('Across shelf distance (km)','FontSize',fs);
  ylabel('Along shelf distance (km)','FontSize',fs);
  %clabel(C,'manual','FontSize',fs,'Color','w')
  title('Channel depth (m) with deepest-draft contours','FontSize',fs);

  %print('-depsc',strcat(['/home/gladish/research/visuals/plume_',flabel,'_channels']));


figure(1);
clf;
hold on;
ratio = zeros(m,n);
for i=1:m
    for j=1:n
        if (grad(i,j) ~= 0)
            ratio(i,j) = log(bmelt(i,j))/grad(i,j);
        end
    end
end
contourf(x/1000.0,y/1000.0,ratio',20,'EdgeColor','None');
colorbar('FontSize',fs);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('bmelt / basal ice gradient','FontSize',fs);
print('-depsc',strcat(['/home/gladish/research/visuals/plume_',flabel,'_bmlt_grad']));
hold off

figure(1);
clf;
hold on;
[m,n,k] = size(dice.thk);
j_start = 3;
%dx = dice.x0(2)-dice.x0(1);
%dy = dice.y0(2)-dice.y0(1);
contourf(dice.Xgrid(:,j_start:n)/1000.0, ...
         dice.Ygrid(:,j_start:n)/1000.0, ...
         min(5.0,dice.vel_div(:,j_start:n,timeslice)), ...
         30,'EdgeColor','None');colorbar('FontSize',fs);
xlabel('Across shelf distance (km)','FontSize',fs)
ylabel('Along shelf distance (km)','FontSize',fs);
colorbar('FontSize',fs);
title('Thickness tendency due to ice divergence','FontSize',fs);    
print('-depsc',strcat(['/home/gladish/research/visuals/plume_',flabel,'_ice_div']));
hold off

clf;
hold on;
j_start = 3;
j_stop = 118;
contourf(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
         dice.Ygrid(:,j_start:j_stop)/1000.0, ...
         -dice.flux_div(:,j_start:j_stop,timeslice), ...
         25,'EdgeColor','None');colorbar('FontSize',fs);
xlabel('Across shelf distance (km)','FontSize',fs)
ylabel('Along shelf distance (km)','FontSize',fs);
title('Thickness time tendency due to flux divergence','FontSize',fs);    
contour(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
        dice.Ygrid(:,j_start:j_stop)/1000.0, ...
        dice.thk(:,j_start:j_stop,timeslice), ...
        20,'w');
caxis([min(min(-dice.flux_div(:,j_start:j_stop,timeslice))) ...
       max(max(-dice.flux_div(:,j_start:j_stop,timeslice)))]);
print('-depsc',strcat(['/home/gladish/research/visuals/plume_',flabel,'_flux_div']));
hold off

clf;
hold on;
j_start = 3;
j_stop = 118;
ice_adv = -dice.flux_div(:,j_start:j_stop,timeslice) + ...
           dice.vel_div(:,j_start:j_stop,timeslice);   
contourf(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
         dice.Ygrid(:,j_start:j_stop)/1000.0, ...
         ice_adv, ...
         25,'EdgeColor','None');colorbar('FontSize',fs);
xlabel('Across shelf distance (km)','FontSize',fs)
ylabel('Along shelf distance (km)','FontSize',fs);
title('Thickness time tendency due to ice advection term u\cdot \nabla H','FontSize',fs);    
%contour(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
%        dice.Ygrid(:,j_start:j_stop)/1000.0, ...
%        dice.thk(:,j_start:j_stop,timeslice), ...
%        20,'w');
caxis([0 40]);
print('-depsc',strcat(['/home/gladish/research/visuals/plume_',flabel,'_ice_adv']));
hold off

clf;
hold on;
j_start = 3;
j_stop = 118;
y_adv = - dice.y_adv(:,j_start:j_stop,timeslice);
   
contourf(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
         dice.Ygrid(:,j_start:j_stop)/1000.0, ...
         y_adv, ...
         25,'EdgeColor','None');colorbar('FontSize',fs);
xlabel('Across shelf distance (km)','FontSize',fs)
ylabel('Along shelf distance (km)','FontSize',fs);
title('Thickness time tendency due to ice advection term v_0*H_y','FontSize',fs);    
%contour(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
%        dice.Ygrid(:,j_start:j_stop)/1000.0, ...
%        dice.thk(:,j_start:j_stop,timeslice), ...
%        20,'w');
caxis([0 40]);
print('-depsc',strcat(['/home/gladish/research/visuals/plume_',flabel,'_y_adv']));
hold off

clf;
hold on;
j_start = 3;
j_stop = 118;
ice_def = - dice.flux_div(:,j_start:j_stop,timeslice) + ...
            dice.y_adv(:,j_start:j_stop,timeslice);
   
contourf(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
         dice.Ygrid(:,j_start:j_stop)/1000.0, ...
         ice_def, ...
         25,'EdgeColor','None');colorbar('FontSize',fs);
xlabel('Across shelf distance (km)','FontSize',fs)
ylabel('Along shelf distance (km)','FontSize',fs);
title('Ice deformation influence on H_t','FontSize',fs);    
%contour(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
%        dice.Ygrid(:,j_start:j_stop)/1000.0, ...
%        dice.thk(:,j_start:j_stop,timeslice), ...
%        20,'w');
%caxis([0 40]);
print('-depsc',strcat(['/home/gladish/research/visuals/plume_',flabel,'_ice_def']));
hold off

clf;
hold on;
j_start = 3;
j_stop = 118;
   
contourf(dice.Xgrid(:,:)/1000.0, ...
         dice.Ygrid(:,:)/1000.0, ...
         dice.thk_t(:,:,timeslice), ...
         25,'EdgeColor','None');colorbar('FontSize',fs);
xlabel('Across shelf distance (km)','FontSize',fs)
ylabel('Along shelf distance (km)','FontSize',fs);
title('H_t/H','FontSize',fs);    
caxis([-0.1 0.1]);
print('-depsc',strcat(['/home/gladish/research/visuals/plume_',flabel,'_thk_t']));
hold off
end
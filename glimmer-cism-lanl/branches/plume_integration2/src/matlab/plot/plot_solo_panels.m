
fig1 = figure(1);
clf;

hold on
contourf(x/1000.0,y/1000.0,grad',40,'EdgeColor','None');
colorbar('FontSize',fs);
set(gca,'FontSize',fs2);
quiver(x(1:stride:end)/1000,y(1:stride:end)/1000,...
      su(1:stride:end,1:stride:end)',...
      sv(1:stride:end,1:stride:end)',...
     scale,'w','LineWidth',lw);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Ocean Velocities with contourfs of ice draft gradient','FontSize',fs);
hold off
fname = strcat([fig_dir,'/plume_grad']);
print('-depsc',fname);

figure(1);
clf;
set(fig1,'Position',[1 1 solo_fig_x_size solo_fig_y_size])

hold on
contourf(x/1000.0,y/1000.0,-draft',40,'EdgeColor','None');
colorbar('FontSize',fs);
set(gca,'FontSize',fs2);
%caxis([000 650])
scale=1.5;stride=2;
i0 = 2;
quiver(x(i0:stride:end)/1000,y(1:stride:end)/1000,...
      su(i0:stride:end,1:stride:end)',...
      sv(i0:stride:end,1:stride:end)',...
      scale,'k','LineWidth',lw);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
caxis([0 600]);
title('Contours of Ice Draft (m) with plume velocities','FontSize',fs);
hold off
print('-depsc',strcat([fig_dir,'/plume_draft_vel']));


figure(1);
set(fig1,'Position',[1 1 solo_fig_x_size solo_fig_y_size])
clf;
hold on
contourf(x/1000.0,y/1000.0,-draft',40,'EdgeColor','None');colorbar('FontSize',fs);
set(gca,'FontSize',fs2);
%caxis([000 650])
%quiver(x(1:stride:end)/1000,y(1:stride:end)/1000,...
%      su(1:stride:end,1:stride:end)',...
%      sv(1:stride:end,1:stride:end)',...
%      scale,'k','LineWidth',lw);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
%caxis([0 800]);
title('Contours of Ice Draft (m)','FontSize',fs);
hold off
print('-depsc',strcat([fig_dir,'/plume_draft']));


figure(1);
set(fig1,'Position',[1 1 solo_fig_x_size solo_fig_y_size])
clf;
%subplot(1,2,2);

hold on
contourf(x/1000.0,y/1000.0,bmelt',20,'EdgeColor','None') ;colorbar('FontSize',fs);
set(gca,'FontSize',fs2);
%contour(x/1000.0,y/1000.0,draft',30,'w');
%caxis([min(min(bmelt)), max(max(bmelt))]);
%caxis([0 40]);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Basal Melt Rate (m/year)','FontSize',fs);
hold off
print('-depsc',strcat([fig_dir,'/plume_meltrate']));


figure(1);
clf;
set(fig1,'Position',[1 1 solo_fig_x_size solo_fig_y_size])
hold on
%subplot(2,2,3)
contourf(x/1000.0,y/1000.0,train'/1000.0,...
         30,'EdgeColor','None') ;colorbar('FontSize',fs);
contour(x/1000.0,y/1000.0,train'/1000.0,[0 0],'k')
set(gca,'FontSize',fs2);
caxis([min(min(train/1000.0)), max(max(train/1000.0))]);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Entrainment/Detrainment Rate (1000m/year)','FontSize',fs);
hold off
print('-depsc',strcat([fig_dir,'/plume_train']));


figure(1);
clf;
set(fig1,'Position',[1 1 solo_fig_x_size solo_fig_y_size])
%subplot(2,2,3);
hold on
contourf(x/1000.0,y/1000.0,temp',...
         30,'EdgeColor','None') ;colorbar('FontSize',fs);
%contour(x/1000.0,y/1000.0,draft',20,'w');
%contour(x/1000.0,y/1000.0,draft',30,'w');
caxis([min(min(temp)), max(max(temp))]);
set(gca,'FontSize',fs2);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Plume Temperature (C)','FontSize',fs);
hold off
print('-depsc',strcat([fig_dir,'/plume_temp']));


figure(1);
clf;
set(fig1,'Position',[1 1 solo_fig_x_size solo_fig_y_size])
%subplot(2,2,3);
hold on
contourf(x/1000.0,y/1000.0,salt',...
         30,'EdgeColor','None') ;colorbar('FontSize',fs);
caxis([min(min(salt)), max(max(salt))]);
set(gca,'FontSize',fs2);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Plume Salinity (psu)','FontSize',fs);
hold off
print('-depsc',strcat([fig_dir,'/plume_salt']));


figure(1);
clf;
set(fig1,'Position',[1 1 solo_fig_x_size solo_fig_y_size])
hold on
contourf(x/1000.0,y/1000.0,T_forcing',...
         30,'EdgeColor','None') ;colorbar('FontSize',fs);
caxis([min(min(T_forcing)), max(max(T_forcing))]);
set(gca,'FontSize',fs2);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Thermal Forcing (C)','FontSize',fs);
hold off
print('-depsc',strcat([fig_dir,'/plume_T_forcing']));


figure(1);
clf;
set(fig1,'Position',[1 1 solo_fig_x_size solo_fig_y_size]);
hold on
contourf(x/1000.0,y/1000.0,channel_amp',20,'EdgeColor','None') ;colorbar('FontSize',fs);
set(gca,'FontSize',fs2);
caxis([min(min(channel_amp)) max(max(channel_amp))]);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Channel depth (m)','FontSize',fs);
print('-depsc',strcat([fig_dir,'/plume_channels']));

if (false)
figure(1);
clf;
set(fig1,'Position',[1 1 solo_fig_x_size solo_fig_y_size]);
hold on;
[m,n,k] = size(dice.thk);
contourf(dice.Xgrid(:,j_start:n)/1000.0, ...
         dice.Ygrid(:,j_start:n)/1000.0, ...
         min(5.0,-dice.vel_div(:,j_start:n,timeslice)), ...
         30,'EdgeColor','None');colorbar('FontSize',fs);
set(gca,'FontSize',fs2);
xlabel('Across shelf distance (km)','FontSize',fs)
ylabel('Along shelf distance (km)','FontSize',fs);
colorbar('FontSize',fs);
title('Thickness tendency due to ice divergence','FontSize',fs);    
print('-depsc',strcat([fig_dir,'/plume_ice_div']));
hold off
end

if (false)
figure(1);
clf;
set(fig1,'Position',[1 1 solo_fig_x_size solo_fig_y_size]);
hold on;

contourf(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
         dice.Ygrid(:,j_start:j_stop)/1000.0, ...
         -dice.flux_div(:,j_start:j_stop,timeslice), ...
         25,'EdgeColor','None');colorbar('FontSize',fs);
set(gca,'FontSize',fs2);
xlabel('Across shelf distance (km)','FontSize',fs)
ylabel('Along shelf distance (km)','FontSize',fs);
title('Thickness time tendency due to flux divergence','FontSize',fs);    
contour(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
        dice.Ygrid(:,j_start:j_stop)/1000.0, ...
        dice.thk(:,j_start:j_stop,timeslice), ...
        20,'w');
caxis([min(min(-dice.flux_div(:,j_start:j_stop,timeslice))) ...
       max(max(-dice.flux_div(:,j_start:j_stop,timeslice)))]);
print('-depsc',strcat([fig_dir,'/plume_flux_div']));
hold off
end


figure(1);
clf;
set(fig1,'Position',[1 1 solo_fig_x_size solo_fig_y_size]);
hold on;
contourf(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
         dice.Ygrid(:,j_start:j_stop)/1000.0, ...
         ice_adv, ...
         25,'EdgeColor','None');colorbar('FontSize',fs);
set(gca,'FontSize',fs2);
xlabel('Across shelf distance (km)','FontSize',fs)
ylabel('Along shelf distance (km)','FontSize',fs);
title('Thickness time tendency due to ice advection term u\cdot \nabla H','FontSize',fs);    
%contour(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
%        dice.Ygrid(:,j_start:j_stop)/1000.0, ...
%        dice.thk(:,j_start:j_stop,timeslice), ...
%        20,'w');
caxis([0 80]);
print('-depsc',strcat([fig_dir,'/plume_ice_adv']));
hold off

figure(1);
clf;
set(fig1,'Position',[1 1 solo_fig_x_size solo_fig_y_size]);
hold on;
contourf(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
         dice.Ygrid(:,j_start:j_stop)/1000.0, ...
         y_adv, ...
         25,'EdgeColor','None');colorbar('FontSize',fs);
set(gca,'FontSize',fs2);
xlabel('Across shelf distance (km)','FontSize',fs)
ylabel('Along shelf distance (km)','FontSize',fs);
title('Thickness time tendency due to ice advection term v_0*H_y','FontSize',fs);    
%contour(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
%        dice.Ygrid(:,j_start:j_stop)/1000.0, ...
%        dice.thk(:,j_start:j_stop,timeslice), ...
%        20,'w');
caxis([0 80]);
print('-depsc',strcat([fig_dir,'/plume_y_adv']));
hold off


figure(1)
clf;
set(fig1,'Position',[1 1 solo_fig_x_size solo_fig_y_size]);
hold on;
   
contourf(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
         dice.Ygrid(:,j_start:j_stop)/1000.0, ...
         ice_def, ...
         25,'EdgeColor','None');colorbar('FontSize',fs);
set(gca,'FontSize',fs2);
xlabel('Across shelf distance (km)','FontSize',fs)
ylabel('Along shelf distance (km)','FontSize',fs);
title('Ice deformation influence on H (m/year)','FontSize',fs);    
%contour(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
%        dice.Ygrid(:,j_start:j_stop)/1000.0, ...
%        dice.thk(:,j_start:j_stop,timeslice), ...
%        20,'w');
%caxis([0 40]);
print('-depsc',strcat([fig_dir,'/plume_ice_def']));
hold off

if (false)
figure(1);
set(fig1,'Position',[1 1 solo_fig_x_size solo_fig_y_size]);
clf;
hold on;
contourf(dice.Xgrid(:,:)/1000.0, ...
         dice.Ygrid(:,:)/1000.0, ...
         dice.thk_t(:,:,timeslice), ...
         25,'EdgeColor','None');colorbar('FontSize',fs);
set(gca,'FontSize',fs2);
xlabel('Across shelf distance (km)','FontSize',fs)
ylabel('Along shelf distance (km)','FontSize',fs);
title('H_t/H','FontSize',fs);    
%caxis([-0.1 0.1]);
print('-depsc',strcat([fig_dir,'/plume_thk_t']));
hold off
end

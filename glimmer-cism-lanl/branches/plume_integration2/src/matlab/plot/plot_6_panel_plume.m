fig1 = figure(1);
clf;
set(fig1,'Position',[1 1 six_panel_fig_x_size six_panel_fig_y_size]);

subplot(2,3,3);
hold on
contourf(x/1000.0,y/1000.0,salt',40,'EdgeColor','None');colorbar('FontSize',fs6);
set(gca,'FontSize',fs6);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs6);
ylabel('Along shelf distance (km)','FontSize',fs6);
title('Plume salinity (psu)','FontSize',fs6);
hold off

subplot(2,3,2);
hold on
contourf(x/1000.0,y/1000.0,pdep',40,'EdgeColor','None') ;colorbar('FontSize',fs6);
set(gca,'FontSize',fs6);
xlabel('Across shelf distance (km)','FontSize',fs6);
ylabel('Along shelf distance (km)','FontSize',fs6);
title('Plume thickness (m)','FontSize',fs6);
hold off

subplot(2,3,6);
contourf(x/1000.0,y/1000.0,temp',40,'EdgeColor','None');colorbar('FontSize',fs6);
set(gca,'FontSize',fs6);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs6);
ylabel('Along shelf distance (km)','FontSize',fs6);
title('Plume temperature (C)','FontSize',fs6);
hold off

subplot(2,3,1);
contourf(x/1000.0,y/1000.0,speed'*100.0,40,'EdgeColor','None');colorbar('FontSize',fs6);
set(gca,'FontSize',fs6);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs6);
ylabel('Along shelf distance (km)','FontSize',fs6);
title('Plume speed (cm/s)','FontSize',fs6);

subplot(2,3,5);
hold on
contourf(x/1000.0,y/1000.0,train'/1000.0,40,'EdgeColor','None');
colorbar('FontSize',fs6);
contour(x/1000.0,y/1000.0,train'/1000.0,[0 0],'k')
caxis([-1000 1000]);
set(gca,'FontSize',fs6);
xlabel('Across shelf distance (km)','FontSize',fs6);
ylabel('Along shelf distance (km)','FontSize',fs6);
title('Entrainment/detrainment (km/year)','FontSize',fs6);
hold off

subplot(2,3,4);
hold on;
contourf(x/1000.0,y/1000.0,bmelt',40,'EdgeColor','None') ;colorbar('FontSize',fs6);
set(gca,'FontSize',fs6);
xlabel('Across shelf distance (km)','FontSize',fs6);
ylabel('Along shelf distance (km)','FontSize',fs6);
title('Melt rate (m/year)','FontSize',fs6);
hold off

fname = strcat([fig_dir,'/plume_',flabel,'_6panel_plume']);
print('-depsc',fname);


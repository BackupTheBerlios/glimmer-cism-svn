fig1 = figure(1);
clf;
%set(fig1,'Position',[1 1 six_panel_fig_x_size six_panel_fig_y_size]);

fs_letter = 18;
fs6 = 16;
x_label = -3;
y_label = 42;

subplot(2,3,3);
hold on
contourf(x/1000.0,y/1000.0,salt',40,'EdgeColor','None');colorbar('FontSize',fs6);
set(gca,'FontSize',fs6);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs6);
ylabel('Along shelf distance (km)','FontSize',fs6);
title('Salinity (psu)','FontSize',fs6);
text(x_label,y_label,'c','color','k','FontSize',fs_letter)
hold off

subplot(2,3,2);
hold on
contourf(x/1000.0,y/1000.0,pdep',40,'EdgeColor','None') ;colorbar('FontSize',fs6);
set(gca,'FontSize',fs6);
xlabel('Across shelf distance (km)','FontSize',fs6);
ylabel('Along shelf distance (km)','FontSize',fs6);
title('Thickness (m)','FontSize',fs6);
text(x_label,y_label,'b','color','k','FontSize',fs_letter)
hold off

subplot(2,3,6);
contourf(x/1000.0,y/1000.0,temp',40,'EdgeColor','None');colorbar('FontSize',fs6);
set(gca,'FontSize',fs6);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs6);
ylabel('Along shelf distance (km)','FontSize',fs6);
title('Temperature (\circC)','FontSize',fs6);
text(x_label,y_label,'f','color','k','FontSize',fs_letter)
hold off

subplot(2,3,1);
contourf(x/1000.0,y/1000.0,speed'*100.0,40,'EdgeColor','None');colorbar('FontSize',fs6);
set(gca,'FontSize',fs6);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs6);
ylabel('Along shelf distance (km)','FontSize',fs6);
title('Speed (cm s^{-1})','FontSize',fs6);
%text(50,50,'a','FontSize',fs);
text(x_label,y_label,'a','color','k','FontSize',fs_letter)

subplot(2,3,5);
hold on
contourf(x/1000.0,y/1000.0,train'/1000.0,40,'EdgeColor','None');
colorbar('FontSize',fs6);
contour(x/1000.0,y/1000.0,train'/1000.0,[0 0],'k')
%caxis([-700 50]);
set(gca,'FontSize',fs6);
xlabel('Across shelf distance (km)','FontSize',fs6);
ylabel('Along shelf distance (km)','FontSize',fs6);
title('Entrain/detrain (km a^{-1})','FontSize',fs6);
%title('$\dot{e}-\dot{d}$ (km/a)','interpreter','latex','FontSize',fs6);
text(x_label,y_label,'e','color','k','FontSize',fs_letter)
hold off

subplot(2,3,4);
hold on;
contourf(x/1000.0,y/1000.0,bmelt',40,'EdgeColor','None') ;colorbar('FontSize',fs6);
set(gca,'FontSize',fs6);
xlabel('Across shelf distance (km)','FontSize',fs6);
ylabel('Along shelf distance (km)','FontSize',fs6);
title('Melt rate (m a^{-1})','FontSize',fs6);
text(x_label,y_label,'d','color','k','FontSize',fs_letter)
hold off

fname = strcat([fig_dir,'/plume_6panel_plume']);
print('-depsc',fname);


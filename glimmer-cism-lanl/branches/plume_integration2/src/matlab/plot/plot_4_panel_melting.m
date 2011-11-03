fig1 = figure(1);
clf;
set(fig1,'Position',[1 1 four_panel_fig_x_size four_panel_fig_y_size]);
subplot(2,2,1);
hold on
fs3 = 14;
contourf(x/1000.0,y/1000.0,bmelt',40,'EdgeColor','None');colorbar('FontSize',fs3);
set(gca,'FontSize',fs3);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs3);
ylabel('Along shelf distance (km)','FontSize',fs3);
title('Basal melt rate (m/year)','FontSize',fs3);
hold off

subplot(2,2,2);
hold on
contourf(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
         dice.Ygrid(:,j_start:j_stop)/1000.0, ...
         ice_def,40,'EdgeColor','None') ;colorbar('FontSize',fs3);
set(gca,'FontSize',fs3);
xlabel('Across shelf distance (km)','FontSize',fs3);
title('Influence of ice deformation (m/year)','FontSize',fs3);
hold off

subplot(2,2,3);
contourf(x/1000.0,y/1000.0,-draft',40,'EdgeColor','None');colorbar('FontSize',fs3);
set(gca,'FontSize',fs3);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs3);
ylabel('Along shelf distance (km)','FontSize',fs3);
title('Contours of Ice Draft (m)','FontSize',fs3);
hold off

subplot(2,2,4);
%contourf(x/1000.0,y/1000.0,-draft',40,'EdgeColor','None');colorbar('FontSize',fs3);
%set(gca,'FontSize',fs3);
%colormap jet;
%xlabel('Across shelf distance (km)','FontSize',fs3);
%ylabel('Along shelf distance (km)','FontSize',fs3);
%title('Contours of Ice Draft (m)','FontSize',fs3);
%hold off

contourf(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
         dice.Ygrid(:,j_start:j_stop)/1000.0, ...
         y_adv, ...
         25,'EdgeColor','None');colorbar('FontSize',fs3);

%contourf(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
%         dice.Ygrid(:,j_start:j_stop)/1000.0, ...
%         ice_adv, ...
%         25,'EdgeColor','None');colorbar('FontSize',fs3);
set(gca,'FontSize',fs2);
xlabel('Across shelf distance (km)','FontSize',fs3)
ylabel('Along shelf distance (km)','FontSize',fs3);
title('Thickness time tendency due to ice advection term v_0 H_y','FontSize',fs3);    

%subplot(2,2,3);
%hold on
%contourf(x/1000.0,y/1000.0,T_forcing',40,'EdgeColor','None') ;colorbar('FontSize',fs33);
%set(gca,'FontSize',fs3);
%xlabel('Across shelf distance (km)','FontSize',fs3);
%title('Thermal forcing (^\circ C)','FontSize',fs3);
%hold off

fname = strcat([fig_dir,'/plume_4panel_melt']);
print('-depsc',fname);


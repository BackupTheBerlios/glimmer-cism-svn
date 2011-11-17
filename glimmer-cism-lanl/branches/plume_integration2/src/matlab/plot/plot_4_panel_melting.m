fig1 = figure(1);
clf;
set(fig1,'Position',[1 1 four_panel_fig_x_size four_panel_fig_y_size]);

subplot(2,2,1);
hold on
fs3 = 14;
fs = 14;
fs_label = 18;

contourf(x/1000.0,y/1000.0,bmelt',40,'EdgeColor','None');colorbar('FontSize',fs3);
set(gca,'FontSize',fs3);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs3);
ylabel('Along shelf distance (km)','FontSize',fs3);
title('Basal melt rate (m/year)','FontSize',fs3);
text(-2,42,'a','color','k','FontSize',fs_label);

hold off

subplot(2,2,2);
hold on
contourf(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
         dice.Ygrid(:,j_start:j_stop)/1000.0, ...
         ice_def,40,'EdgeColor','None') ;colorbar('FontSize',fs3);
set(gca,'FontSize',fs3);
xlabel('Across shelf distance (km)','FontSize',fs3);
ylabel('Along shelf distance (km)','FontSize',fs3);
title('Influence of ice deformation (m/year)','FontSize',fs3);
text(-2,31,'b','color','k','FontSize',fs_label);

hold off

subplot(2,2,3);
contourf(x/1000.0,y/1000.0,-draft',40,'EdgeColor','None');colorbar('FontSize',fs3);
set(gca,'FontSize',fs3);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs3);
ylabel('Along shelf distance (km)','FontSize',fs3);
title('Contours of Ice Draft (m)','FontSize',fs3);
text(-2,42,'c','color','k','FontSize',fs_label);
hold off

subplot(2,2,4);
hold on
contourf(dice.Xgrid(:,j_start:j_stop)/1000.0, ...
         dice.Ygrid(:,j_start:j_stop)/1000.0, ...
         y_adv, ...
         25,'EdgeColor','None');colorbar('FontSize',fs3);

set(gca,'FontSize',fs2);
xlabel('Across shelf distance (km)','FontSize',fs3)
ylabel('Along shelf distance (km)','FontSize',fs3);
title('Slab-like advection v_0 H_y','FontSize',fs3);    
text(-2,31,'d','color','k','FontSize',fs_label);
hold off


fname = strcat([fig_dir,'/plume_4panel_melt']);
print('-depsc',fname);


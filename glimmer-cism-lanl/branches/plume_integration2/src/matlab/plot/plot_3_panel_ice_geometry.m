three_panel_fig_x_size = 1500;
three_panel_fig_y_size = 600;

fig1 = figure(1);
clf;

fs3 = 11;

set(fig1,'Position',[1 1 three_panel_fig_x_size three_panel_fig_y_size]);
subplot(1,3,3);
hold on
contourf(x/1000.0,y/1000.0,channel_amp',40,'EdgeColor','None');colorbar('FontSize',fs3);
set(gca,'FontSize',fs3);
colormap jet;
%[C,h] = contour(x/1000.0,y/1000.0,max_draft',10,'w');
caxis([min(min(channel_amp)) max(max(channel_amp))]);
xlabel('Across shelf distance (km)','FontSize',fs3);

%clabel(C,'manual','FontSize',fs3,'Color','w')
%title('Channel depth (m) with deepest-draft contours','FontSize',fs3);
title('Channel depth (m)','FontSize',fs3);
hold off

subplot(1,3,2);
hold on
contourf(x/1000.0,y/1000.0,bmelt',20,'EdgeColor','None') ;colorbar('FontSize',fs3);
set(gca,'FontSize',fs3);
xlabel('Across shelf distance (km)','FontSize',fs3);
title('Basal Melt Rate (m/year)','FontSize',fs3);
hold off

subplot(1,3,1);
hold on
contourf(x/1000.0,y/1000.0,-draft',40,'EdgeColor','None') ;colorbar('FontSize',fs3);
set(gca,'FontSize',fs3);
xlabel('Across shelf distance (km)','FontSize',fs3);
ylabel('Along shelf distance (km)','FontSize',fs3);
title('Ice draft (m)','FontSize',fs3);
hold off

fname = strcat([fig_dir,'/plume_3panel']);
print('-depsc',fname);


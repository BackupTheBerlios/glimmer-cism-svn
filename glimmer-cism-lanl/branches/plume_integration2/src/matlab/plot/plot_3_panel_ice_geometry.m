three_panel_fig_x_size = 1400;
three_panel_fig_y_size = 600;

fig1 = figure(1);
clf;

fs3 = 16;
fs_label = 18;

set(fig1,'Position',[1 1 three_panel_fig_x_size three_panel_fig_y_size]);
subplot(1,3,3);
hold on
contourf(x/1000.0,y/1000.0,channel_amp',40,'EdgeColor','None');colorbar('FontSize',fs3);
set(gca,'FontSize',fs3);
colormap jet;
%%[C,h] = contour(x/1000.0,y/1000.0,max_draft',10,'w');
caxis([min(min(channel_amp)) max(max(channel_amp))]);
xlabel('Across shelf distance (km)','FontSize',fs3);
text(-2,41,'c','color','k','FontSize',fs_label);

%clabel(C,'manual','FontSize',fs3,'Color','w')
%title('Channel depth (m) with deepest-draft contours','FontSize',fs3);
title('Channel depth (m)','FontSize',fs3);
hold off


subplot(1,3,2);
%subplot(1,2,2);
hold on
contourf(x/1000.0,y/1000.0,bmelt',20,'EdgeColor','None') ;colorbar('FontSize',fs3);
set(gca,'FontSize',fs3);
xlabel('Across shelf distance (km)','FontSize',fs3);
title('Basal Melt Rate (m/a)','FontSize',fs3);
text(-2,41,'b','color','k','FontSize',fs_label);
hold off

subplot(1,3,1);
%subplot(1,2,1);
hold on
contourf(x/1000.0,y/1000.0,-draft',40,'EdgeColor','None') ;colorbar('FontSize',fs3);

scale = 1.35;
scale = 1.5;
arrowcolor = 'k';
stride = 3;
arrow_start = 2;
lw = 1.25;
quiver(x(arrow_start:stride:end)/1000,y(arrow_start:stride:end)/1000,...
      su(arrow_start:stride:end,arrow_start:stride:end)',...
      sv(arrow_start:stride:end,arrow_start:stride:end)',...
      scale,arrowcolor,'LineWidth',lw,'ShowArrowHead','on');
rectangle('Position',[15 36.5 4 2.5],'Facecolor','w');
h = quiver(18.5,37,0,1,scale,arrowcolor,'LineWidth',lw);
adjust_quiver_arrowhead_size(h,4);
text(15.5,37.5,'1 m/s');
set(gca,'FontSize',fs3);
xlabel('Across shelf distance (km)','FontSize',fs3);
ylabel('Along shelf distance (km)','FontSize',fs3);
title('Ice draft (m)','FontSize',fs3);
text(-2,41,'a','color','k','FontSize',fs_label);

hold off

fname = strcat([fig_dir,'/plume_3panel']);
print('-depsc',fname);


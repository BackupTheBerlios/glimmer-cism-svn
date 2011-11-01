figure(1);
clf;
set(fig1,'Position',[1 1 solo_fig_x_size solo_fig_y_size]);
hold on
contourf(x/1000.0,y/1000.0,-draft',40,'EdgeColor','None');colorbar('FontSize',fs);
set(gca,'FontSize',fs2);
quiver(dice.xstag(1:stride2:end)/1000, ...
       dice.ystag(1:stride2:end)/1000,...
       ice_u(1:stride2:end,1:stride2:end)',...
       ice_v(1:stride2:end,1:stride2:end)'-0.0,...
       scale,'k','LineWidth',lw);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
caxis([0 600]);
title('Contours of Ice Draft (m) with ice velocities','FontSize',fs);
hold off
print('-depsc',strcat([fig_dir,'/plume_',flabel,'_draft_ice_vel']));



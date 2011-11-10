figure(1);
clf;
fs3 = 14;
subplot(1,2,1);
hold on
contourf(x/1000.0,y/1000.0,-draft',40,'EdgeColor','None');colorbar('FontSize',fs3);
set(gca,'FontSize',fs3);
%caxis([000 650])
quiver(x(1:stride:end)/1000,y(1:stride:end)/1000,...
      su(1:stride:end,1:stride:end)',...
      sv(1:stride:end,1:stride:end)',...
      scale,'k','LineWidth',lw);
colormap jet;
xlabel('Across shelf distance (km)','FontSize',fs3);
ylabel('Along shelf distance (km)','FontSize',fs3);
caxis([0 600]);
title('Contours of Ice Draft (m) with plume velocities','FontSize',fs3);
hold off

subplot(1,2,2);
%plot(y/1000.0,mean(bmelt,1),'b-','LineWidth',3.0);
xlabel('Along shelf distance (km)','FontSize',fs3);
ylabel('melt rate (m/year)','FontSize',fs3);
%title('Cross-shelf averaged melt rate (m/year)','FontSize',fs3);

[ax1,h1,h2] = plotyy(dice.ygrid/1000.0,mean(dice.thk(:,:,end),1), ...
                     dice.ygrid/1000.0,mean(dice.bmlt(:,:,end),1));
set(ax1,'ycolor','k')
  lw=2.0;
set(h1,'Linewidth',lw,'LineStyle','-','color','k');
set(h2,'Linewidth',lw,'LineStyle','-','color','b');

xmin=-2.0;
xmax=40.0;
y2max = 60.5;
y1max = 605;
set(ax1(1),'xlim',[xmin xmax]);
set(ax1(2),'xlim',[xmin xmax]);
set(ax1(1),'ylim',[0 y1max]);
set(ax1(2),'ylim',[0 y2max]);
set(ax1(1),'FontSize',fs3);
set(ax1(2),'FontSize',fs3);
set(ax1(1),'YTick',0:100:y1max);
set(ax1(2),'YTick',0:10:y2max);
xlabel('Along shelf distance (km)','FontSize',fs3);
set(get(ax1(1),'Ylabel'),'String','ice thickness (m)','FontSize',fs3)
set(get(ax1(2),'Ylabel'),'String','melt rate (m/year)','FontSize',fs3)
  title('Cross-shelf averaged ice thickness and basal melt rates','FontSize',fs3);


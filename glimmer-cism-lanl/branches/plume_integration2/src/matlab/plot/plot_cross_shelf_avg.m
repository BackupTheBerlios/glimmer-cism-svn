figure(1);
clf;
fs3 = 14;


[ax1,h1,h2] = plotyy(dice.ygrid/1000.0,mean(dice.thk(:,:,end),1), ...
                     dice.ygrid/1000.0,mean(dice.bmlt(:,:,end),1));
set(ax1(1),'ycolor','k')
set(ax1(2),'ycolor','b')

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

function plot_amb_water_column(ts,ss,zs,fig_dir)

% plot the ambient water column
figure(1);
clf;
set(figure(1),'Position',[1 1 600 600]);
fs = 14;
lw = 2.0;

[ax1,h1,h2] = plotyy(ts,zs,ss,zs); 

set(ax1(1),'xcolor','b')
set(ax1(2),'xcolor','r')
set(ax1(1),'ycolor','k')
set(ax1(2),'ycolor','k')


set(h1,'Linewidth',lw,'LineStyle','-','color','b');
set(h2,'Linewidth',lw,'LineStyle','-','color','r');

set(ax1(1),'xlim',[-2.2 0.25]);
set(ax1(2),'xlim',[34.65 34.77]);
set(ax1(1),'XAxisLocation','Bottom');
set(ax1(2),'XAxisLocation','Top');

set(ax1(1),'ylim',[0 max(zs)]);
set(ax1(2),'ylim',[0 max(zs)]);

set(ax1(1),'FontSize',fs);
set(ax1(2),'FontSize',fs);

set(ax1(1),'YTick',0:100:max(zs));
set(ax1(2),'YTick',0:100:max(zs));
set(ax1(1),'XTick',-2:0.5:0.5);
set(ax1(2),'XTick',34.66:0.02:34.8);

set(get(ax1(1),'Ylabel'),'String','Depth (m)','FontSize',fs)
set(get(ax1(2),'Ylabel'),'String','Depth (m)','FontSize',fs)
set(ax1(1),'Ydir','Reverse');
set(ax1(2),'Ydir','Reverse');

set(get(ax1(1),'Xlabel'),'String','Temperature (\circC)','FontSize',fs);
set(get(ax1(2),'Xlabel'),'String','Salinity (psu)','FontSize',fs);

print ('-depsc',strcat([fig_dir,'/plume_amb_water_column']));

end

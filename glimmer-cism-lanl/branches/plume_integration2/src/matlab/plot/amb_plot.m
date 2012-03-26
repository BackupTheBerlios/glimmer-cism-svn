load piecewise_amb;

figure(1);
clf;
tstride = 2;
sstride = 4;

tdepth = piecewise_amb(1:tstride:end,1);
sdepth = piecewise_amb(1:sstride:end,1);
temp = piecewise_amb(1:tstride:end,2);
salt = piecewise_amb(1:sstride:end,3);

saltclr = 'k';
tempclr = 'k';
saltstyle = '.';
tempstyle = '-';

fs = 13;

ax1 = axes('XAxisLocation','top',...
           'YAxisLocation','left',...
           'Color', 'w',...
           'XColor',tempclr,'YColor','k',...
           'FontSize',fs);
       
h1 = line(temp,-tdepth,'Color',tempclr,'Linestyle',tempstyle,'linewidth',1.0);
h11 =line([-1 -0.1],[0 -800],'Color','r');
%ax1 = gca;
%set(ax1,'XAxisLocation','top')
%set(ax1,'XColor',tempclr,'YColor','k','FontSize',fs);
set(get(ax1,'XLabel'),'String', 'Temperature ( ^\circ C)','FontSize',fs);
set(get(ax1,'Ylabel'),'String', 'Depth (m)','color','k','FontSize',fs);
%ylim([-800 0]);
set(ax1,'ylim',[-825 25]);
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','bottom',...
           'YAxisLocation','right',...
           'Color', 'none',...
           'XColor',saltclr,'YColor','k',...
           'FontSize',fs);
       
h2 = line(salt,-sdepth,'Color',saltclr,'Parent',ax2,'linestyle',saltstyle);
set(get(ax2,'XLabel'),'String', 'Salinity (psu)','FontSize',fs);  
set(ax2,'ylim',[-825 25]);
set(ax2,'xlim',[34.62 34.76]);
%print('-depsc','amb_water.eps');

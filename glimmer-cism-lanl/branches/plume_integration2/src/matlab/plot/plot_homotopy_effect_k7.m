


recalculate = false;
doFigure1 = false;
doFigure2 = true;
doFigure3 = false;

if (recalculate)
%clear all
jobs = getenv('GC_JOBS');

j1name = 'oct30_perturb_usq.9_thk_perturb_k_7.0_amp_25.0';

f1_ice = strcat([jobs,'/',j1name,'/',j1name,'.out.nc']);
f1_ocean = strcat([jobs,'/',j1name,'/plume.',j1name,'.out.nc']);

f2_ice = strcat([jobs,'/',j1name,'/','usq.9.no_homotopy.out.nc']);
f2_ocean = strcat([jobs,'/',j1name,'/plume.usq.9.no_homotopy.out.nc']);

f1_ice
dice1 = nc_ice_read(f1_ice,0,1,675);
f2_ice
dice2 = nc_ice_read(f2_ice,0,1,19);

f1_ocean
dplume1 = nc_plume_read(f1_ocean,0,1,675);
f2_ocean
dplume2 = nc_plume_read(f2_ocean,0,1,19);

m_control = 0.789;
[m,m_applied,in,out,acab,unsteady] = mass_balance(dice1,dplume1,-1);
m1 = m/in;
  m1_a = m_applied/in;
[m,m_applied,in,out,acab,unsteady] = mass_balance(dice2,dplume2,-1);
m2 = m/in;
  m2_a = m_applied/in;

end

 if (doFigure1)
fig1 = figure(1);
set(fig1,'Position',[1 1 1200 800]);
fs = 16;
fs_label=18;
clf;

subset = @(M) M(:,1:120);

subplot(2,2,1);
hold on
set(gca,'FontSize',fs);
xlabel('across shelf distance (km)','FontSize',fs);
ylabel('along shelf distance (km)','FontSize',fs);
contourf(subset(dice1.Xgrid/1000),subset(dice1.Ygrid/1000), ...
         subset(dice1.bmlt(:,:,6)),30,'EdgeColor','None');
colorbar('FontSize',fs);
caxis([0 80]);

subplot(2,2,2);
hold on
set(gca,'FontSize',fs);
xlabel('across shelf distance (km)','FontSize',fs);
ylabel('along shelf distance (km)','FontSize',fs)
hold on
contourf(subset(dice1.Xgrid/1000),subset(dice1.Ygrid/1000),...
         subset(dice2.bmlt(:,:,end)),30,'EdgeColor','None');
colorbar('FontSize',fs);
caxis([0 80]);

subplot(2,2,3);
hold on
set(gca,'FontSize',fs);
xlabel('across shelf distance (km)','FontSize',fs);
ylabel('along shelf distance (km)','FontSize',fs)
hold on
contourf(subset(dice1.Xgrid/1000),subset(dice1.Ygrid/1000), ...
         subset(dice1.lsurf(:,:,6)),30,'EdgeColor','None');
colorbar('FontSize',fs);
caxis([-550 0]);

subplot(2,2,4);
hold on
set(gca,'FontSize',fs);
xlabel('across shelf distance (km)','FontSize',fs);
ylabel('along shelf distance (km)','FontSize',fs)
contourf(subset(dice2.Xgrid/1000),subset(dice2.Ygrid/1000), ...
         subset(dice2.lsurf(:,:,end)),30,'EdgeColor','None');
colorbar('FontSize',fs);
caxis([-550 0]);
end

if (false)
[melt,applied_melt,influx,outflux,total_acab,dvdt] = mass_balance(dice1,dplume1,2:length(dplume1.time));
[melt2,applied_melt2,influx2,outflux2,total_acab2,dvdt2] = mass_balance(dice2,dplume2,2:length(dplume2.time));
end

if (doFigure2)

fig2 =figure(2);
clf;

set(fig2,'Position',[1 1 800 600]);
t = dplume1.time(2:end)-dplume1.time(1);
t2 = dplume2.time(2:end)-dplume2.time(1);

lw=2.0;
fs = 16;
ms = 15;

hl3 = line(t, melt/influx(1));

hl2 = line(t, applied_melt/influx(1));
hl2b = line(t2, applied_melt2/influx(1));

ax1 = gca;

c1 = 'b';
c2 = 'r';

set(ax1,'YColor','b');

ax2 = axes('Position',get(ax1,'Position'),...
	   'XAxisLocation','bottom',...
	   'YAxisLocation','right',...
	   'Color','None',...
	   'XColor','k','YColor','k');


set(ax1,'FontSize',fs);
set(ax2,'FontSize',fs);

hl1 = line(t,dvdt/influx(1));
hl1b = line(t2,dvdt2/influx(1));


set(ax1,'XColor','k','Ycolor',c1);
set(ax2,'XColor','k','Ycolor','k');

set(hl1,'MarkerSize',ms,'LineStyle','-','Color','k');
set(hl1b,'MarkerSize',ms,'LineStyle','.','Color','k');

set(hl2,'MarkerSize',ms,'LineStyle','-','Color',c1);
set(hl2b,'MarkerSize',ms,'LineStyle','.','Color',c1);

set(hl3,'MarkerSize',ms,'LineStyle','-','Color',c2);



set(hl1,'LineWidth',lw);
set(hl1b,'LineWidth',lw);
set(hl2,'LineWidth',lw);
set(hl2b,'LineWidth',lw);
set(hl3,'LineWidth',lw);

set(ax2,'xlim',[-5 170]);
set(ax1,'xlim',[-5 170]);
set(ax2,'ylim',[-15 50]/100);
set(ax1,'ylim',[25 100]/100);
%set(ax2,'YTick',(-5:5:30)/100);
%set(ax1,'YTick',(60:5:100)/100);

xlabel('time (years)','FontSize',fs);
set(get(ax2,'Ylabel'),'String','Volume rate of change (% of influx)','FontSize',fs)
set(get(ax1,'Ylabel'),'String','Total basal melting (% of influx)','FontSize',fs)

%title('Cross-shelf averaged ice thickness and basal melt rates','FontSize',fs);

set(gca,'FontSize',16);

end

if (doFigure3)
figure(3);
clf;
hold on
plot(t,applied_melt);
plot(t2,applied_melt2);
end

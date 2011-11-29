jobs = getenv('GC_JOBS');
jname = 'central_paper';
f_ice = strcat([jobs,'/',jname,'/',jname,'.out.nc']);
f_ocean = strcat([jobs,'/',jname,'/plume.',jname,'.out.nc']);


dice = nc_ice_read(f_ice,-1,1,-1);

fig1 = figure(1);
set(fig1,'Position',[1 1 800 600]);
fs = 16;

subplot(2,1,2);
contourf(dice.Ygrid'/1000,dice.Xgrid'/1000,-dice.lsurf(:,:,end)',20,'EdgeColor','None');
colorbar;caxis([0 550]);
ylabel('Across shelf distance (km)','FontSize',fs)
xlabel('Along shelf distance (km)','FontSize',fs)
set(gca,'FontSize',fs);

subplot(2,1,1);
plot(dice.Ygrid(1,3:end)/1000,dice.bmlt(10,3:end,end),'k.');
xlabel('Along shelf distance (km)','FontSize',fs)
ylabel('Basal melt rate (m/a)','FontSize',fs);
xlim([-1 41]);
set(gca,'FontSize',fs);

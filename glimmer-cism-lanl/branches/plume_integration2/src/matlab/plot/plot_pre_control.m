jobs = getenv('GC_JOBS');
jname = 'central_paper';
f_ice = strcat([jobs,'/',jname,'/',jname,'.out.nc']);
f_ocean = strcat([jobs,'/',jname,'/plume.',jname,'.out.nc']);

dice = nc_ice_read(f_ice,-1,1,-1,1);

fig1 = figure(1);
set(fig1,'Position',[1 1 800 600]);
fs = 16;

subplot(1,2,2);

contourf(dice.Xgrid'/1000,dice.Ygrid'/1000,-dice.lsurf(:,:,end)',20,'EdgeColor','None');
colorbar;caxis([0 550]);
xlabel('Across shelf distance (km)','FontSize',fs)
ylabel('Along shelf distance (km)','FontSize',fs)
text(17,36,'b','FontSize',18);
set(gca,'FontSize',fs);

subplot(1,2,1);
plot(dice.Ygrid(1,3:end)/1000,dice.bmlt(10,3:end,end),'k-','LineWidth',2);
xlabel('Along shelf distance (km)','FontSize',fs)
ylabel('Basal melt rate (m a^{-1})','FontSize',fs);
text(35,23,'a','FontSize',18);
xlim([-1 41]);
set(gca,'FontSize',fs);

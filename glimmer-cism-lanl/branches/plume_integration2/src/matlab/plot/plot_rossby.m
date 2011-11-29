jobs = getenv('GC_JOBS');
jname = 'oct25_high_min_visc_smooth_5000.0_k12amp_25.0_restart_4';

f_ocean = strcat([jobs,'/',jname,'/plume.',jname,'.out.nc']);

f_amb = strcat([jobs,'/',jname,'/ambout']);
[amb_t,amb_s] = nc_read_amb(f_amb);

istart = 0;
istride = 1;
iend = 36;


dplume = nc_plume_read(f_ocean,istart,stride,iend);
dplume_avg = nc_plume_avg(dplume);

res = plume_momentum_bal(dplume_avg,amb_t,amb_s);

x = dplume_avg.x/1000;
y = dplume_avg.y/1000;

f = 2*pi*sin(80*pi/180) / (3600*24.0);

rossby_rad = sqrt(res.g_prime .* dplume_avg.pdep) / (1000*f);
bmlt = dplume_avg.bmelt;

fig2 = figure(2);
set(fig2,'Position',[1 1 1000 600]);

fs = 16;

subplot(1,2,1);
contourf(x,y,rossby_rad',30,'EdgeColor','None');
xlabel('Across shelf (km)','FontSize',fs);
ylabel('Along shelf (km)','FontSize',fs);
set(gca,'FontSize',fs);
colorbar('FontSize',fs);

subplot(1,2,2);
contourf(x,y,bmlt',30,'EdgeColor','None')
xlabel('Across shelf (km)','FontSize',fs);
ylabel('Along shelf (km)','FontSize',fs);
colorbar('FontSize',fs);
set(gca,'FontSize',fs);

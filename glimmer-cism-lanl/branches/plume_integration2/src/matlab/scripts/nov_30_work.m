jobs = getenv('GC_JOBS');
jname = 'oct25_high_min_visc_smooth_5000.0_k12amp_25.0_restart_4';
jfast = 'no_tangle_oct30_perturb_usq.2_v_inflow_1250';
jvfast = 'no_tangle_oct30_perturb_usq.2_v_inflow_1500';
jslow = 'no_tangle_oct30_perturb_usq.2_v_inflow_750';
jmslow = 'no_tangle_oct30_perturb_usq_v_inflow_900';

iname = @(j) strcat([jobs,'/',j,'/',j,'.out.nc']);
pname = @(j) strcat([jobs,'/',j,'/plume.',j,'.out.nc']);

f_ice = iname(jname);
f_ice_fast = iname(jfast);
f_ice_vfast = iname(jvfast);
f_ice_slow = iname(jslow);
f_ice_mslow = iname(jmslow);

f_ocean = pname(jname);
f_ocean_fast = pname(jfast);
f_ocean_vfast = pname(jvfast);
f_ocean_slow = pname(jslow);
f_ocean_mslow = pname(jmslow);

dice =      nc_ice_read(f_ice,      -1,1,-1,2);
dice_fast = nc_ice_read(f_ice_fast, -1,1,-1,2);
dice_vfast =nc_ice_read(f_ice_vfast,-1,1,-1,2);
dice_slow = nc_ice_read(f_ice_slow, -1,1,-1,2);
dice_mslow =nc_ice_read(f_ice_mslow,-1,1,-1,2);

dplume = nc_plume_read(f_ocean,-1,1,-1);
dplume_fast = nc_plume_read(f_ocean_fast,-1,1,-1);
dplume_vfast = nc_plume_read(f_ocean_vfast,-1,1,-1);
dplume_slow = nc_plume_read(f_ocean_slow,-1,1,-1);
dplume_mslow = nc_plume_read(f_ocean_mslow,-1,1,-1);

[m750,am,in,out,acab,dvdt] = mass_balance(dice_slow,dplume_slow,-1);
m750
[m900,am,in,out,acab,dvdt] = mass_balance(dice_mslow,dplume_mslow,-1);
m900
[m1000,am,in,out,acab,dvdt] = mass_balance(dice,dplume,-1);
m1000
[m1250,am,in,out,acab,dvdt] = mass_balance(dice_fast,dplume_fast,-1);
m1250
[m1500,am,in,out,acab,dvdt] = mass_balance(dice_vfast,dplume_vfast,-1);
m1500

fs = 16;

fig1 = figure(1);
clf;
set(fig1,'Position',[1 1 800 600]);

hl1 = line([750,900,1000,1250,1500],[m750,m900,m1000,m1250,m1500]*1.0e-9);
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'),...
	   'XAxisLocation','bottom',...
           'YAxislocation','right',...
            'Color','None',...
	   'XColor','k','YColor','b');

set(ax1,'FontSize',fs);
set(ax2,'FontSize',fs);

hl2 = line([750,900,1000,1250,1500], ...
        [m750,m900,m1000,m1250,m1500] ./  ...
	   ([750,900,1000,1250,1500]*(12e9/1000)));

set(ax1,'XColor','k','Ycolor','b');
set(ax2,'XColor','k','Ycolor','r');

lw = 2.0;
ms = 10.0;
set(hl1,'MarkerSize',ms,'Marker','*','LineStyle','-','LineWidth',lw,'Color','b');
set(hl2,'MarkerSize',ms,'Marker','*','LineStyle','-','LineWidth',lw,'Color','r');

set(ax2,'ylim',[0.6 1]);
set(ax2,'ytick',0:0.1:1);
set(ax1,'ylim',[5 15]);
set(ax1,'ytick',5:2:15);
set(ax1,'xlim',[700 1550]);
set(ax2,'xlim',[700 1550]);

set(get(ax1,'ylabel'),'String','Total melt rate (km^3/a)');
set(get(ax1,'ylabel'),'FontSize',fs);

set(get(ax2,'ylabel'),'String','Total melt rate (percent of influx)');
set(get(ax2,'ylabel'),'FontSize',fs);
set(get(ax1,'xlabel'),'String','Inflow velocity of ice (m/a)','FontSize',fs);

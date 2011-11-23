clear all

jobs = getenv('GC_JOBS');

j1name = 'oct25_high_min_visc_smooth_5000.0_k12amp_25.0_restart_4';

j2name = 'oct30_perturb_usq.8_top_-1.55';
j3name = 'oct30_perturb_usq.8_bottom_0.4';
j4name = 'no_tangle_oct30_perturb_usq_bottom_-0.2';
j5name = 'no_tangle_oct30_perturb_usq_bottom_0.0';


f1_ice = strcat([jobs,'/',j1name,'/',j1name,'.out.nc']);
f1_ocean = strcat([jobs,'/',j1name,'/plume.',j1name,'.out.nc']);

f2_ice = strcat([jobs,'/',j2name,'/',j2name,'.out.nc']);
f2_ocean = strcat([jobs,'/',j2name,'/plume.',j2name,'.out.nc']);

f3_ice = strcat([jobs,'/',j3name,'/',j3name,'.out.nc']);
f3_ocean = strcat([jobs,'/',j3name,'/plume.',j3name,'.out.nc']);

f4_ice = strcat([jobs,'/',j4name,'/',j4name,'.out.nc']);
f4_ocean = strcat([jobs,'/',j4name,'/plume.',j4name,'.out.nc']);

f5_ice = strcat([jobs,'/',j5name,'/',j5name,'.out.nc']);
f5_ocean = strcat([jobs,'/',j5name,'/plume.',j5name,'.out.nc']);

dice1 = nc_ice_read(f1_ice,-1,1,-1);
dice2 = nc_ice_read(f2_ice,-1,1,-1);
dice3 = nc_ice_read(f3_ice,-1,1,-1);
dice4 = nc_ice_read(f4_ice,-1,1,-1);
dice5 = nc_ice_read(f5_ice,-1,1,-1);

dplume1 = nc_plume_read(f1_ocean,-1,1,-1);
dplume2 = nc_plume_read(f2_ocean,-1,1,-1);
dplume3 = nc_plume_read(f3_ocean,-1,1,-1);
dplume4 = nc_plume_read(f4_ocean,-1,1,-1);
dplume5 = nc_plume_read(f5_ocean,-1,1,-1);

lsurf_control     = dice1.lsurf(:,:,end);
lsurf2 = dice2.lsurf(:,:,end);
lsurf3 = dice3.lsurf(:,:,end);
lsurf4 = dice4.lsurf(:,:,end);
lsurf5 = dice5.lsurf(:,:,end);

[m,in,out,acab,unsteady] = mass_balance(dice1,dplume1);
m_control = m/in
[m,in,out,acab,unsteady] = mass_balance(dice2,dplume2);
m2 = m/in - m_control
[m,in,out,acab,unsteady] = mass_balance(dice3,dplume3);
m3 = m/in - m_control
[m,in,out,acab,unsteady] = mass_balance(dice4,dplume4);
m4 = m/in - m_control
[m,in,out,acab,unsteady] = mass_balance(dice5,dplume5);
m5 = m/in - m_control


fig1 = figure(1);
set(fig1,'Position',[1 1 1200 800]);
fs = 16;
fs_label=18;
clf;

y0 = -500;
N = length(flatten_field(lsurf_control));
f = 0.5;
subset_fun = random_subset_fun(f,N);
subset_fun = @(v) v;
ytext = 'perturbed ice draft (m)';
xtext = 'original ice draft (m)';


subplot(2,2,1);
hold on
plot([y0 0],[y0 0],'k-','LineWidth',2.5);
plot(subset_fun(flatten_field(lsurf_control)), ...
     subset_fun(flatten_field(lsurf3)),'r.');
title(strcat(['total melt change = ',sprintf('%+2.1f',100*m3),' percent of influx']),'FontSize',fs);

set(gca,'FontSize',fs);
plot([-425 -425],[y0 0],'k--');
text(-415,-75,'0.3 ^\circ C bottom warming\newline upper limit','FontSize',fs);
text(-575,50,'a.','FontSize',fs_label);

xlabel(xtext,'FontSize',fs);
ylabel(ytext,'FontSize',fs)
xlim([y0 0]);
ylim([y0 0]);


subplot(2,2,2);
hold on
plot([y0 0],[y0 0],'k-','LineWidth',2.5);
plot(subset_fun(flatten_field(lsurf_control)), ...
     subset_fun(flatten_field(lsurf2)),'r.');

plot([-150 -150],[y0 0],'k--');
text(-250,-350,'0.3 ^\circ C upper warming \newline lower limit','FontSize',fs);
text(-575,50,'b.','FontSize',fs_label);
set(gca,'FontSize',fs);
xlabel(xtext,'FontSize',fs);
ylabel(ytext,'FontSize',fs)
xlim([y0 0]);
ylim([y0 0]);
title(strcat(['total melt change = ',sprintf('%+2.1f',100*m2),' percent of influx']),'FontSize',fs);

subplot(2,2,4);
hold on
plot([y0 0],[y0 0],'k-','LineWidth',2.5);
plot(subset_fun(flatten_field(lsurf_control)), ...
     subset_fun(flatten_field(lsurf4)),'b.');
set(gca,'FontSize',fs);
plot([-425 -425],[y0 0],'k--');
text(-475,-100,'0.3 ^\circ C bottom cooling\newline upper limit','FontSize',fs);
text(-575,50,'d.','FontSize',fs_label);
xlabel(xtext,'FontSize',fs);
ylabel(ytext,'FontSize',fs)
xlim([y0 0]);
ylim([y0 0]);

title(strcat(['total melt change = ',sprintf('%+2.1f',100*m4),' percent of influx']),'FontSize',fs);

subplot(2,2,3);
hold on
plot([y0 0],[y0 0],'k-','LineWidth',2.5);
plot(subset_fun(flatten_field(lsurf_control)), ...
     subset_fun(flatten_field(lsurf5)),'b.');
set(gca,'FontSize',fs);

plot([-425 -425],[y0 0],'k--');
text(-415,-100,'0.1 ^\circ C bottom cooling \newline upper limit','FontSize',fs);
text(-575,50,'c.','FontSize',fs_label);
xlabel(xtext,'FontSize',fs);
ylabel(ytext,'FontSize',fs)
xlim([y0 0]);
ylim([y0 0]);
title(strcat(['total melt change = ',sprintf('%+2.1f',100*m5),' percent of influx']),'FontSize',fs);



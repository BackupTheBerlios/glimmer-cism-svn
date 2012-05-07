clear all

jobs = getenv('GC_JOBS');

control_name = 'oct25_high_min_visc_smooth_5000.0_k12amp_25.0_restart_4';

j_sgd_1p0 = 'no_tangle_oct30_perturb_usq_sgd_flux_1.0';
j_sgd_1p25 = 'oct30_perturb_usq.8_sgd_flux_1.25';
j_sgd_0p5 = 'oct30_perturb_usq.8_sgd_flux_0.5';
j_sgd_1p5 = 'oct30_perturb_usq.8_sgd_flux_1.5';
j_sgd_2p0 = 'no_tangle_oct30_perturb_usq_sgd_flux_2.0';

f_control_ice = strcat([jobs,'/',control_name,'/',control_name,'.out.nc']);
f_control_ocean = strcat([jobs,'/',control_name,'/plume.',control_name,'.out.nc']);

f_0p5_ice = strcat([jobs,'/',j_sgd_0p5,'/',j_sgd_0p5,'.out.nc']);
f_0p5_ocean = strcat([jobs,'/',j_sgd_0p5,'/plume.',j_sgd_0p5,'.out.nc']);
f_1p0_ice = strcat([jobs,'/',j_sgd_1p0,'/',j_sgd_1p0,'.out.nc']);
f_1p0_ocean = strcat([jobs,'/',j_sgd_1p0,'/plume.',j_sgd_1p0,'.out.nc']);
f_1p25_ice = strcat([jobs,'/',j_sgd_1p25,'/',j_sgd_1p25,'.out.nc']);
f_1p25_ocean = strcat([jobs,'/',j_sgd_1p25,'/plume.',j_sgd_1p25,'.out.nc']);
f_1p5_ice = strcat([jobs,'/',j_sgd_1p5,'/',j_sgd_1p5,'.out.nc']);
f_1p5_ocean = strcat([jobs,'/',j_sgd_1p5,'/plume.',j_sgd_1p5,'.out.nc']);
f_2p0_ice = strcat([jobs,'/',j_sgd_2p0,'/',j_sgd_2p0,'.out.nc']);
f_2p0_ocean = strcat([jobs,'/',j_sgd_2p0,'/plume.',j_sgd_2p0,'.out.nc']);

minslices = 2;
dice_control = nc_ice_read(f_control_ice,-1,1,-1,minslices);
dplume_control = nc_plume_read(f_control_ocean,-1,1,-1);

dice_0p5 = nc_ice_read(f_0p5_ice,-1,1,-1,minslices);
dice_1p0 = nc_ice_read(f_1p0_ice,-1,1,-1,minslices);
dice_1p25 = nc_ice_read(f_1p25_ice,-1,1,-1,minslices);
dice_1p5 = nc_ice_read(f_1p5_ice,-1,1,-1,minslices);
dice_2p0 = nc_ice_read(f_2p0_ice,-1,1,-1,minslices);


dplume_0p5 = nc_plume_read(f_0p5_ocean,-1,1,-1);
dplume_1p0 = nc_plume_read(f_1p0_ocean,-1,1,-1);
dplume_1p25 = nc_plume_read(f_1p25_ocean,-1,1,-1);
dplume_1p5 = nc_plume_read(f_1p5_ocean,-1,1,-1);
dplume_2p0 = nc_plume_read(f_2p0_ocean,-1,1,-1);


lsurf_control     = dice_control.lsurf(:,:,end);
lsurf_0p5 = dice_0p5.lsurf(:,:,end);
lsurf_1p0 = dice_1p0.lsurf(:,:,end);
lsurf_1p25 = dice_1p25.lsurf(:,:,end);
lsurf_1p5 = dice_1p5.lsurf(:,:,end);
lsurf_2p0 = dice_2p0.lsurf(:,:,end);

[m,in,out,acab,unsteady] = mass_balance(dice_control,dplume_control,-1);
m_control = m/in
[m,in,out,acab,unsteady] = mass_balance(dice_0p5,dplume_0p5,-1);
m_0p5 = m/in - m_control
[m,in,out,acab,unsteady] = mass_balance(dice_1p0,dplume_1p0,-1);
m_1p0 = m/in - m_control;
[m,in,out,acab,unsteady] = mass_balance(dice_1p25,dplume_1p25,-1);
m_1p25 = m/in - m_control;
[m,in,out,acab,unsteady] = mass_balance(dice_1p5,dplume_1p5,-1);
m_1p5 = m/in - m_control;
[m,in,out,acab,unsteady] = mass_balance(dice_2p0,dplume_2p0,-1);
m_2p0 = m/in - m_control;

fig1 = figure(1);
set(fig1,'Position',[1 1 1200 800]);
fs = 16;
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

plot(subset_fun(flatten_field(lsurf_control)), ...
     subset_fun(flatten_field(lsurf_0p5)),'b.');
plot([y0 0],[y0 0],'k-','LineWidth',2.5);
title(strcat(['total melt change = ',sprintf('%+2.1f',100*m_0p5),' percent of influx']),'FontSize',fs);
text(-450,-50,'a) 0.5 km^3/a ','FontSize',fs);
set(gca,'FontSize',fs);
%plot([-425 -425],[y0 0],'k--');
%text(-415,-75,'0.3 ^\circ C bottom warming\newline upper limit','FontSize',fs);
xlabel(xtext,'FontSize',fs);
ylabel(ytext,'FontSize',fs)
xlim([y0 0]);
ylim([y0 0]);


subplot(2,2,2);
hold on

plot(subset_fun(flatten_field(lsurf_control)), ...
     subset_fun(flatten_field(lsurf_1p0)),'b.');
plot([y0 0],[y0 0],'k-','LineWidth',2.5);
%plot([-150 -150],[y0 0],'k--');
%text(-250,-350,'0.3 ^\circ C upper warming \newline lower limit','FontSize',fs);
text(-450,-50,'b) 1.0 km^3/a ','FontSize',fs);
set(gca,'FontSize',fs);
xlabel(xtext,'FontSize',fs);
ylabel(ytext,'FontSize',fs)
xlim([y0 0]);
ylim([y0 0]);
title(strcat(['total melt change = ',sprintf('%+2.1f',100*m_1p0),' percent of influx']),'FontSize',fs);

subplot(2,2,4);
hold on

plot(subset_fun(flatten_field(lsurf_control)), ...
     subset_fun(flatten_field(lsurf_2p0)),'b.');
plot([y0 0],[y0 0],'k-','LineWidth',2.5);
set(gca,'FontSize',fs);
text(-450,-50,'d) 2.0 km^3/a ','FontSize',fs);
xlabel(xtext,'FontSize',fs);
ylabel(ytext,'FontSize',fs)
xlim([y0 0]);
ylim([y0 0]);

title(strcat(['total melt change = ',sprintf('%+2.1f',100*m_2p0),' percent of influx']),'FontSize',fs);

subplot(2,2,3);
hold on

plot(subset_fun(flatten_field(lsurf_control)), ...
     subset_fun(flatten_field(lsurf_1p5)),'b.');
plot([y0 0],[y0 0],'k-','LineWidth',2.5);
set(gca,'FontSize',fs);

text(-450,-50,'c) 1.5 km^3/a ','FontSize',fs);

xlabel(xtext,'FontSize',fs);
ylabel(ytext,'FontSize',fs)
xlim([y0 0]);
ylim([y0 0]);
title(strcat(['total melt change = ',sprintf('%+2.1f',100*m_1p5),' percent of influx']),'FontSize',fs);



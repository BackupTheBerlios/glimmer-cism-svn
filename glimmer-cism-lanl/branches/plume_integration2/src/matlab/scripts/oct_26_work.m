
fs = 14;

jobs = getenv('GC_JOBS');
jname = 'oct25_high_min_visc_smooth_5000.0_k12amp_25.0_restart_1';

f_ice = strcat([jobs,'/',jname,'/',jname,'.out.nc']);
f_ocean = strcat([jobs,'/',jname,'/plume.',jname,'.out.nc']);
f_amb = strcat([jobs,'/',jname,'/ambout']);

istart = 625;
iend = 700;
stride = 1;
%dice = nc_ice_read(f_ice,  istart,stride,iend);

istart = 250;
iend = 280;
stride = 1;
%dplume = nc_plume_read(f_ocean,istart,stride,iend);

%dplume_avg = nc_plume_avg(dplume);
%dice_avg = nc_ice_avg(dice);

%[flat_ocean_transient,flat_ice_transient] = flatten_gc( dplume, dice, 1:length(dice.time) );
%[flat_ocean,flat_ice] = flatten_gc( dplume_avg, dice_avg,  ...
%				    length(dice_avg.time):length(dice_avg.time) );

%save('/home/cvg222/matlab_figs/oct_27_work_dat.mat');
load '/home/cvg222/matlab_figs/oct_27_work_dat.mat'

[amb_t,amb_s] = nc_read_amb(f_amb);
zs = 0:1:800;
ts = amb_t(zs);
ss = amb_s(zs);

% plot the ambient water column
%figure(10);
%clf;
%plot(ts,zs,'k'); set(gca,'ydir','reverse');
%xlabel('ambient water temperature (^\circ C)','fontsize',fs);
%ylabel('depth (m)','fontsize',fs);
%set(gca,'fontsize',fs);

% plot keel crossing
%times = [1];
%plot_keel_crossing(dice_avg,dplume_avg,times);
%plot_draft_sections(dice_avg,2+[1,5,21,41,81,141]);
%plot_plume2(dplume_avg,dice_avg,-1,'oct_27');
%plot_spectral_evolution(dice_avg);
%gc_scatter(flat_ocean,flat_ice,flat_ocean,flat_ice,amb_t);
%plot_geostrophic(dplume_avg,1,amb_t,amb_s)
  plume_momentum_bal(dplume_avg,amb_t,amb_s)


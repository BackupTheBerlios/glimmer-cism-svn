%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% script to carry out some work for plume paper   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

jobs = getenv('GC_JOBS');
jname = 'oct25_high_min_visc_smooth_5000.0_k12amp_25.0_restart_2';
jname = 'central_paper_no_channels_restart_3';

f_ice = strcat([jobs,'/',jname,'/',jname,'.out.nc']);
f_ocean = strcat([jobs,'/',jname,'/plume.',jname,'.out.nc']);

istart = 440;
iend = 441;
stride = 1;
dice = nc_ice_read(f_ice,  istart,stride,iend);

%istart = 330;
%iend = 350;
%stride = 1;
dplume = nc_plume_read(f_ocean,istart,stride,iend);

f_amb = strcat([jobs,'/',jname,'/ambout']);
[amb_t,amb_s] = nc_read_amb(f_amb);
zs = 0:1:800;
ts = amb_t(zs);
ss = amb_s(zs);

dplume_avg = nc_plume_avg(dplume);
dice_avg = nc_ice_avg(dice);

%[flat_ocean_transient,flat_ice_transient] = flatten_gc( dplume, dice, 1:length(dice.time) );
[flat_ocean,flat_ice] = flatten_gc( dplume_avg, dice_avg,  ...
				    length(dice_avg.time):length(dice_avg.time) );

%save('/home/cvg222/paper_work/mat_files/nov_1_work_dat.mat');
%save('/home/cvg222/paper_work/mat_files/nov_2_nochannels.mat');
load '/home/cvg222/paper_work/mat_files/nov_1_work_dat.mat'
%load '/home/cvg222/paper_work/mat_files/nov_2_nochannels.mat';
fig_dir = '/home/cvg222/paper_work/no_channels';

%plot_amb_water_column(ts,ss,zs,fig_dir);

% plot keel crossing
times = [1];
%plot_keel_crossing(dice_avg,dplume_avg,times,fig_dir);

%plot_draft_sections(dice_avg,2+[1,5,21,41,81,141],fig_dir);

plot_plume2(dplume_avg,dice_avg,-1,fig_dir);

%plot_spectral_evolution(dice_avg,fig_dir);

gc_scatter(flat_ocean,flat_ice,flat_ocean,flat_ice,amb_t,fig_dir);

%plot_geostrophic(dplume_avg,1,amb_t,amb_s)

%plume_momentum_bal(dplume_avg,amb_t,amb_s)


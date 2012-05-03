%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% script to carry out some work for plume paper   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;

setenv('GC_JOBS','/archive/cvg222/gc_output/2011/paper_jobs/');

jobs = getenv('GC_JOBS');
jname = 'central_final_nov17';
%jname = 'oct25_high_min_visc_smooth_5000.0_k12amp_25.0_restart_4';
%jname = 'no_tangle_oct30_perturb_usq_bottom_-0.2';
%jname = 'central_paper';
%jname = 'no_tangle_oct30_perturb_usq_bottom_-0.2';

f_ice = strcat([jobs,'/',jname,'/',jname,'.out.nc']);
f_ocean = strcat([jobs,'/',jname,'/plume.',jname,'.out.nc']);

%istart = 1990;
istart = -1;
%iend = 90;
iend = -1;
stride = 1;
dice = nc_ice_read(f_ice,istart,stride,iend);

%istart = 15;
istart = -1;
%iend = 35;
iend = -1;
stride = 1;
dplume = nc_plume_read(f_ocean,istart,stride,iend);

f_amb = strcat([jobs,'/',jname,'/ambout']);
[amb_t,amb_s] = nc_read_amb(f_amb);
zs = 0:1:800;
ts = amb_t(zs);
ss = amb_s(zs);

dplume_avg = nc_plume_avg(dplume);
dice_avg = nc_ice_avg(dice);

%[flat_ocean_transient,flat_ice_transient] = flatten_gc( dplume, dice, 1:length(dice.time) );
%[flat_ocean,flat_ice] = flatten_gc( dplume_avg, dice_avg,  ...
%				    length(dice_avg.time):length(dice_avg.time) );

%load '/home/cvg222/paper_work/mat_files/nov_15_work_dat.mat';
%save '/home/cvg222/paper_work/mat_files/nov_20_work_dat.mat';
%load '/home/cvg222/paper_work/mat_files/nov_20_work_dat.mat';

fig_dir = '/home/cvg222/paper_work/apr_16_figs/';

%plot_amb_water_column(ts,ss,zs,fig_dir);

% plot keel crossing
%times = [1];
%plot_keel_crossing(dice_avg,dplume_avg,times,fig_dir);

%plot_draft_sections(dice_avg,2+[1,5,21,41,81,141],fig_dir);

plot_plume2(dplume_avg,dice_avg,-1,fig_dir);

%plot_spectral_evolution(dice_avg,fig_dir);

%gc_scatter(flat_ocean,flat_ice,flat_ocean,flat_ice,amb_t,fig_dir);

%plot_geostrophic(dplume_avg,1,amb_t,amb_s)


%res = plume_momentum_bal(dplume_avg,amb_t,amb_s);
%plot_momentum_bal;

%res = plume_momentum_bal2(dplume_avg,amb_t,amb_s);
%plot_momentum_bal2;


%plot_temp_salt_depths;

%plot_plume_section;
%fs = 16;
%text(20.1,-570,'y=1 km','FontSize',fs);
%text(20.1,-510,'y=3 km','FontSize',fs);
%text(20.1,-430,'y=5 km','FontSize',fs);
%text(20.1,-360,'y=7 km','FontSize',fs);
%text(20.1,-290,'y=9 km','FontSize',fs);
%text(20.0,-190,'y=15 km','FontSize',fs);
%text(20.0,-110,'y=40 km','FontSize',fs);%

%plot_lower_surf;

%plume_res = plume_dep_bal(dplume_avg);
%plot_dep_bal;


%plot_k_perturb_total_melt;


%plot_h_diff;

%plot_top_warming;

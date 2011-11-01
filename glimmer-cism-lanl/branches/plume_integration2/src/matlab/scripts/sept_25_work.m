

fs = 14;

jobs = getenv('GC_JOBS');
jname = 'test_k_12_300_detrain_ustar_2_amp_30.0_tau_25';
jname = 'test_k_12_ustar_2_amp_30.0_tau_25';
jname = 'central_final_sept25';
jname = 'central_paper_homtopy_frac_0.5_ramp_1.0_pmin_10.0';

istart = 0;
iend = 374;
stride = 25;

f_ice = strcat([jobs,'/',jname,'/',jname,'.out.nc']);
f_ocean = strcat([jobs,'/',jname,'/plume.',jname,'.out.nc']);
f_amb = strcat([jobs,'/',jname,'/ambout']);

%dice =     nc_ice_read(f_ice,  istart,stride,iend);
%dplume = nc_plume_read(f_ocean,istart,stride,iend);

%dplume2 = nc_plume_read(f_ocean,365,1,375);
%dplume.train_avg = mean(dplume2.train,3);

%[flat_ocean_transient,flat_ice_transient] = flatten_gc( dplume, dice, 1:length(dice.time) );
%[flat_ocean,flat_ice] = flatten_gc( dplume, dice,  length(dice.time):length(dice.time) );

%save('/home/cvg222/matlab_figs/oct_17_work_dat.mat');
load '/home/cvg222/matlab_figs/oct_17_work_dat.mat';

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
%plot_keel_crossing(dice,dplume,[4]);

figure(11);
clf;
%plot_draft_sections(dice,[1,5,17,37,77,157]);
%plot_draft_sections(dice,2+[1,5,21,41,81,141]);
%pause;

%plot_plume2(dplume,dice,-1,'poster_oct_19');

%plot_spectral_evolution(dice);

gc_scatter(flat_ocean,flat_ice,flat_ocean_transient,flat_ice_transient,amb_t);



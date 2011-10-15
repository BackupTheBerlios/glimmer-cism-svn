fs = 14;

jobs = getenv('GC_JOBS');
jname = 'test_k_12_300_detrain_ustar_2_amp_30.0_tau_25';
jname = 'test_k_12_ustar_2_amp_30.0_tau_25';
jname = 'central_final_sept25';

maxslices = 5;
stride = 1;

f_ice = strcat([jobs,'/',jname,'/',jname,'.out.nc']);
f_ocean = strcat([jobs,'/',jname,'/plume.',jname,'.out.nc']);
f_amb = strcat([jobs,'/',jname,'/ambout']);

%dice = nc_ice_read(f_ice,stride,maxslices);
%dplume = nc_plume_read(f_ocean,stride,maxslices);

%[flat_ocean_transient,flat_ice_transient] = flatten_gc( dplume, dice, 1:length(dice.time) );
%[flat_ocean,flat_ice] = flatten_gc( dplume, dice,  length(dice.time):length(dice.time) );

%save('/home/cvg222/matlab_figs/sept_25_work_dat.mat');
load '/home/cvg222/matlab_figs/sept_25_work_dat.mat'

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
%plot_keel_crossing(dice,dplume,[1]);

%figure(11);
%plot_draft_sections(dplume,[1,5,17,37,77,157]);
%pause;

%plot_plume2(dplume,dice,-1,'central');

%gc_scatter(flat_ocean,flat_ice,flat_ocean_transient,flat_ice_transient,amb_t);



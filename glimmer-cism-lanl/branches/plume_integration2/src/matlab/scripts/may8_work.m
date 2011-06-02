f_ice_trial_8 = '/data/gladish/may3_coupled_super_channels_trial_8/may3_coupled_super_channels_trial_8.out.nc';
f_plume_ice_trial_8 = '/data/gladish/may3_coupled_super_channels_trial_8/plume.may3_coupled_super_channels_trial_8.out.nc';

f_ice_trial_90 = '/data/gladish/may3_coupled_super_channels_trial_90/may3_coupled_super_channels_trial_90.out.nc';
f_plume_ice_trial_90 = '/data/gladish/may3_coupled_super_channels_trial_90/plume.may3_coupled_super_channels_trial_90.out.nc';

f_ice_trial_8_solid = '/data/gladish/may3_coupled_super_channels_trial_8_restart_1/may3_coupled_super_channels_trial_8_restart_1.out.nc';
f_plume_trial_8_solid = '/data/gladish/may3_coupled_super_channels_trial_8_restart_1/plume.may3_coupled_super_channels_trial_8_restart_1.out.nc';

f_plume_3_1 = '/data/gladish/may3.1_diff_10.0_cdb_2.5_dtime_14400.0_entype_6_pthk_10.0_tvel_0.0/plume.may3.1_diff_10.0_cdb_2.5_dtime_14400.0_entype_6_pthk_10.0_tvel_0.0.out.nc';
f_ice_3_1 =   '/data/gladish/may3.1_diff_10.0_cdb_2.5_dtime_14400.0_entype_6_pthk_10.0_tvel_0.0/may3.1_diff_10.0_cdb_2.5_dtime_14400.0_entype_6_pthk_10.0_tvel_0.0.out.nc';

%d_plume_trial_8 = nc_plume_read(f_plume_ice_trial_8);
%d_ice_trial_8 = nc_read(f_ice_trial_8);

%d_plume_trial_8_solid = nc_plume_read(f_plume_trial_8_solid);
%d_ice_trial_8_solid = nc_read(f_ice_trial_8_solid);

d_plume_trial_90 = nc_plume_read(f_plume_ice_trial_90);
%d_ice_trial_90 = nc_read(f_ice_trial_90);

%d_plume_3_1 = nc_plume_read(f_plume_3_1);
%d_ice_3_1 = nc_read(f_ice_3_1);

%plot_plume2(d_plume_trial_8_solid,d_ice_trial_8_solid,230,'trial_8_solid');
%plot_plume2(d_plume_trial_8_solid,d_ice_trial_8_solid,1,'trial_8_steady');
%plot_plume2(d_plume_3_1,d_ice_trial_90,230,'3_1_final');


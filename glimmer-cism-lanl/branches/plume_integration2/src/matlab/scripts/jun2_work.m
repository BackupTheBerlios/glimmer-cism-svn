f_ice1 = '/data/gladish/gc_output/igs_jobs/may18.2_coupled_temp_-20.0_tau_25000.0_trial_3/may18.2_coupled_temp_-20.0_tau_25000.0_trial_3.out.nc';
f_plume1 = '/data/gladish/gc_output/igs_jobs/may18.2_coupled_temp_-20.0_tau_25000.0_trial_3/plume.may18.2_coupled_temp_-20.0_tau_25000.0_trial_3.out.nc';

f_ice2 = '/data/gladish/gc_output/igs_jobs/may3_coupled_super_channels_trial_46/may3_coupled_super_channels_trial_46.out.nc';
f_plume2 = '/data/gladish/gc_output/igs_jobs/may3_coupled_super_channels_trial_46/plume.may3_coupled_super_channels_trial_46.out.nc';

f_ice3 = '/data/gladish/gc_output/gfdl_jobs/may3_coupled_super_channels_trial_8/may3_coupled_super_channels_trial_8.out.nc';
f_plume3 = '/data/gladish/gc_output/gfdl_jobs/may3_coupled_super_channels_trial_8/plume.may3_coupled_super_channels_trial_8.out.nc';

f_ice4 = '/data/gladish/gc_output/igs_jobs/may26.T_profile_highres_20min_temp_-10.0_tau_25000.0_diff_5.0_pmin_20.0_trial_10/may26.T_profile_highres_20min_temp_-10.0_tau_25000.0_diff_5.0_pmin_20.0_trial_10.out.nc';
f_plume4 = '/data/gladish/gc_output/igs_jobs/may26.T_profile_highres_20min_temp_-10.0_tau_25000.0_diff_5.0_pmin_20.0_trial_10/plume.may26.T_profile_highres_20min_temp_-10.0_tau_25000.0_diff_5.0_pmin_20.0_trial_10.out.nc';

f_ice5 = '/data/gladish/gc_output/igs_jobs/may26.T_profile_highres_20min_temp_-20.0_tau_25000.0_diff_5.0_pmin_20.0_trial_1/may26.T_profile_highres_20min_temp_-20.0_tau_25000.0_diff_5.0_pmin_20.0_trial_1.out.nc';
f_plume5='/data/gladish/gc_output/igs_jobs/may26.T_profile_highres_20min_temp_-20.0_tau_25000.0_diff_5.0_pmin_20.0_trial_1/plume.may26.T_profile_highres_20min_temp_-20.0_tau_25000.0_diff_5.0_pmin_20.0_trial_1.out.nc';

f_ice6 = '/data/gladish/gc_output/igs_jobs/may29_highres_temp_-10.0_tau_25000.0_diff_10.0_k_4.0_amp_50.0_tempbot_-0.1_pmin_10.0/may29_highres_temp_-10.0_tau_25000.0_diff_10.0_k_4.0_amp_50.0_tempbot_-0.1_pmin_10.0.out.nc';
f_plume6='/data/gladish/gc_output/igs_jobs/may29_highres_temp_-10.0_tau_25000.0_diff_10.0_k_4.0_amp_50.0_tempbot_-0.1_pmin_10.0/plume.may29_highres_temp_-10.0_tau_25000.0_diff_10.0_k_4.0_amp_50.0_tempbot_-0.1_pmin_10.0.out.nc';

dice1 = nc_ice_read(f_ice2,10,64);
dplume1 = nc_plume_read(f_plume2,10,64);

dice2 = nc_ice_read(f_ice2,10,126);
dplume2 = nc_plume_read(f_plume2,10,126);

dice3 = nc_ice_read(f_ice3,10,130);
dplume3 = nc_plume_read(f_plume3,10,130);

dice4 = nc_ice_read(f_ice4,10,35);
dplume4 = nc_plume_read(f_plume4,40,35);

dice5 = nc_ice_read(f_ice5,10,33);
dplume5 = nc_plume_read(f_plume5,40,33);

dice6 = nc_ice_read(f_ice6,10,51);
dplume6 = nc_plume_read(f_plume6,10,51);
jobs = getenv('GC_JOBS');

jname = 'oct25_high_min_visc_smooth_5000.0_k12amp_25.0_restart_4';
f_amb = strcat([jobs,'/',jname,'/ambout']);
[amb_t,amb_s] = nc_read_amb(f_amb);
zs = 0:1:800;
ts = amb_t(zs);
ss = amb_s(zs);


jname2 = 'no_tangle_oct30_perturb_usq_bottom_0.2';
f_amb2 = strcat([jobs,'/',jname2,'/ambout']);
[amb_t,amb_s] = nc_read_amb(f_amb2);
zs2 = 0:1:800;
ts2 = amb_t(zs);
ss2 = amb_s(zs);

jname3 = 'no_tangle_oct30_perturb_usq_top_-1.35';
f_amb3 = strcat([jobs,'/',jname3,'/ambout']);
[amb_t,amb_s] = nc_read_amb(f_amb3);
zs3 = 0:1:800;
ts3 = amb_t(zs);
ss3 = amb_s(zs);

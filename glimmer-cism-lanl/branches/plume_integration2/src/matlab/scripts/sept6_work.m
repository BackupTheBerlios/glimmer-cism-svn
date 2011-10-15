f_ice = '/scratch/cvg222/gcp/GC_jobs_paper/may29_highres_original/ice.out.avg.nc';
f_plume = '/scratch/cvg222/gcp/GC_jobs_paper/may29_highres_original/plume.out.avg.nc';

f_ice = '/scratch/cvg222/gcp/GC_jobs_paper/may29_highres_original/ice.out.nc';
f_plume = '/scratch/cvg222/gcp/GC_jobs_paper/may29_highres_original/plume.out.nc';

f_ice = '/scratch/cvg222/gcp/GC_jobs_paper/may29_newamb/ice.out.nc';
f_plume = '/scratch/cvg222/gcp/GC_jobs_paper/may29_newamb/plume.out.nc';


dice = nc_ice_read(f_ice,1,10);
dplume = nc_plume_read(f_plume,1,10);

[amb_t,amb_s] = nc_read_amb('/scratch/cvg222/gcp/GC_jobs_paper/may29_newamb/ambout');

%plot_plume2(dplume,dice,1,'k4_steady');

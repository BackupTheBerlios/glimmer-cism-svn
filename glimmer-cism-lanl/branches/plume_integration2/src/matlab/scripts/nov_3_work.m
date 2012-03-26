f = '/data/gladish/oct25_high_min_visc_smooth_10000.0_k12amp_25.0/oct25_high_min_visc_smooth_10000.0_k12amp_25.0.out.nc';
g = '/data/gladish/oct25_high_min_visc_smooth_10000.0_k12amp_25.0/plume.oct25_high_min_visc_smooth_10000.0_k12amp_25.0.out.nc';

%dice = nc_ice_read(f, 310,1,315);
dplume=nc_plume_read(g, 1200,1,1250);

contourf(dplume.xgrid/1000,dplume.ygrid/1000,dplume.vorticity(:,:,end));colorbar;


m = 24
n = 34
kinbcw = 2
ramp_width = 0
landw = 0
n_level = 7
hx = 100.0
hy = 100.0
dt = 0.05
tend = 1.0
kx = 1.0
chan_depth = 10.0

ifdep = 600.0
gldep = 600.0
wcdep = 200.0
up_thk = gldep
up_vel = -1000.0


run_name = 'test_restart'

plume_nl_filename = '%s.nl' % run_name
gc_config_filename = '%s.config' % run_name

nc_input_fname = '%s.in.nc' % run_name

# This is to do a restart
input_nc_file = 'restart.nc'
input_nc_file_t_read = 3

nc_regrid_args = [input_nc_file,
                  '%s.in.nc'% run_name, input_nc_file_t_read,
                  m,n,0,0,0,0,0,0,0,0,
                  2,0,0,0]
#nc_regrid_args = []

#nc_gen_input_args = ['gs',nc_input_fname, m,n,kinbcw, n_level,
#                     hx,hy,up_thk, up_vel, ifdep,kx, chan_depth,ramp_width,landw]
nc_gen_input_args = []

plume_vals = {'ifdep' : ifdep, 
              'plume_min_thickness' : 10.0,
              'm_grid' : m,
              'n_grid' : n,
              'hx' : hx,
              'hy' : hy,
              'dt1' : 150.0,
              'phi' : 0.0,
              'horturb' : True,
              'ah' : 1.0,
              'kh' : 1.0,
              'tempbot' : -1.0,
              'temptop' : -1.0,
              'saltbot' : 34.5,
              'salttop' : 34.5,
              'gldep' : gldep,
              'wcdep' : wcdep}

gc_vals = {'plume' : { 'plume_imin' : 1,
 			'plume_imax' : m,
			'plume_kmin' : 5,
	        	'plume_kmax' : n,
                        'plume_steadiness_tol' : 1.0e-3,
                       'plume_write_all_states' : False,
			'plume_output_prefix' : run_name,
                       'plume_min_spinup_time':1,
                       'plume_output_dir' : '/scratch/plume_netcdf_output',
			'plume_output_file' : 'plume.%s.out.nc' % run_name,
	                'plume_nl_file' : plume_nl_filename},
            'time' : {'tend' : tend,
 		      'dt' : dt },
	    'grid' : { 'upn' : n_level,
			'ewn' : m,
			'nsn' : n,
			'dew' : hx,
			'dns' : hy },
           'parameters' : { 'ice_limit' : 10.0 },
           'Petermann shelf' : {'accumulation_rate' : 0.0,
                                'air_temperature' : -5.0},
 	    'CF output' : { 'frequency' : dt,
			    'name' : '%s.out.nc' % run_name},
            'CF input' : { 'name' : nc_input_fname,
	                   'time' : 1 },
	    'CF default' : { 'comment' : '',
			     'title' : run_name },
            'options' : { 'flow_law' : 0,
                          'use_plume' : 0,
                          },

	}


# Sample plume-GC configuration file

# Note: Uses python syntax, so this file must be a valid python file

m = 46
n = 46
ramp_width = 5
n_level = 5
hx = 200.0
hy = 200.0
dt = 0.025
tend = 0.05
kx = 2.0
chan_depth = 50.0

ifdep = 400.0
gldep = 600.0
wcdep = 200.0
up_thk = gldep
up_vel = -10.0

run_name = 'ghost_shelf_test'

plume_nl_filename = '%s.nl' % run_name
gc_config_filename = '%s.config' % run_name
nc_input_fname = '%s.in.nc' % run_name		

nc_regrid_args = []

nc_gen_input_args = ['gs',
                    nc_input_fname,
		    m,n, n_level,hx,hy,
                    up_thk, up_vel, ifdep,
                    kx, chan_depth,ramp_width
                    ]

plume_vals = {'ifdep' : ifdep, 
              'plume_min_thickness' : 10.0,
              'm_grid' : m,
              'n_grid' : n,
              'hx' : hx,
              'hy' : hy,
              'dt1' : 50.0,
              'phi' : 0.0,
              'tempbot' : 0.0,
              'temptop' : 0.0,
              'saltbot' : 34.5,
              'salttop' : 34.5,
              'gldep' : gldep,
              'wcdep' : wcdep}

gc_vals = {'plume' : { 'plume_imin' : 1,
 			'plume_imax' : m,
			'plume_kmin' : 1,
	        	'plume_kmax' : n-2,
                        'plume_steadiness_tol' : 1.0e-4,
			'plume_output_prefix' : run_name,
			'plume_output_file' : 'plume.%s.out.nc' % run_name,
	                'plume_nl_file' : plume_nl_filename},
            'time' : {'tend' : tend,
 		      'dt' : dt },
	    'grid' : { 'upn' : n_level,
			'ewn' : m,
			'nsn' : n,
			'dew' : hx,
			'dns' : hy },
           'Petermann shelf' : {'accumulation_rate' : 0.0,
                                'air_temperature' : -5.0},
 	    'CF output' : { 'frequency' : dt,
			    'name' : '%s.out.nc' % run_name},
            'CF input' : { 'name' : nc_input_fname,
	                   'time' : 1 },
	    'CF default' : { 'comment' : '',
			     'title' : run_name },
           'options' : { 'flow_law' : 2,
                         'use_plume' : 0},

	}


# Sample plume-GC configuration file

# Note: Uses python syntax, so this file must be a valid python file

m = 46
n = 86
hx = 200.0
hy = 200.0
dt = 0.25
tend = 1.0


ifdep = 400.0
gldep = 600.0
wcdep = 200.0

run_name = 'petermann_combined_test'

plume_nl_filename = '%s.nl' % run_name
gc_config_filename = '%s.config' % run_name
nc_input_fname = '%s.in.nc' % run_name		

nc_gen_input_args = ['cs',
                    nc_input_fname,
		    m,n, hx,hy,
                    80, #start of shelf
                    700.0, #grounded ice thickness
                    2, # end of shelf
		    ifdep, #ice front ice thickness
                    -gldep, # ocean topography
                    200.0, # land topography
 	            0.0,0.0,0.0, #no channels
                    3, # 3 cells of kinbc on sides
                    False, # shelf is open to the south
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
	    'grid' : { 'upn' : 5,
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
                         'which_bmlt' : 1},

	}


# Sample plume-GC configuration file

# Note: Uses python syntax, so this file must be a valid python file

m = 20
n = 20
hx = 100.0
hy = 100.0
run_name = 'petermann_combined_test'

plume_nl_filename = '%s.nl' % run_name
gc_config_filename = '%s.config' % run_name
nc_input_fname = '%s.in.nc' % run_name		

nc_gen_input_args = ['cs',
		     nc_input_fname,
		     m,n, hx,hy, 10, 700.0, n-2, 100.0,500.0,200.0,
 	            0.0,0.0,0.0,3,True]

plume_vals = {'ifdep' : 1000.0, 
	'plume_min_thickness' : 10.0,
	 'm_grid' : m,
	 'n_grid' : n,
	 'hx' : hx,
	 'hy' : hy,
	 'dt1' : 50.0,
	 'gldep' : 1000.0,
	 'wcdep' : 500.0}

gc_vals = {'plume' : { 'plume_imin' : 1,
   			'plume_imax' : m,
			'plume_kmin' : 1,
	        	'plume_kmax' : n-2,
			'plume_output_prefix' : 'test_prefix',
			'plume_output_file' : run_name,
	                'plume_nl_file' : 'test_pnl.nl'},
            'time' : {'tend' : 1.0,
 		      'dt' : 0.25 },
	    'grid' : { 'upn' : 11,
			'ewn' : 46,
			'nsn' : 86,
			'dew' : 200.0,
			'dns' : 200.0 },
 	    'CF output' : { 'frequency' : 0.25,
			    'name' : 'test.out.nc' },
            'CF input' : { 'name' : 'test.in.nc',
	                   'time' : 1 },
	    'CF default' : { 'comment' : 'test',
			     'title' : 'test case' },

	}


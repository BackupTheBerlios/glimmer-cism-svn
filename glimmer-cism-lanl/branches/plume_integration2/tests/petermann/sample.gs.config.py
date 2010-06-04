# Sample plume-GC configuration file

# Note: Uses python syntax, so this file must be a valid python file

_run_name = 'c15'

_m = 41
_n = 81
_kinbcw = 2
_ramp_width = 0
_landw = 0
_n_level = 5
_hx = 250.0
_hy = 250.0
_kx = 1.0
_chan_depth = 20.0

_dt = 0.05
_tstart = 0.0
_tend = 1.0

_ifpos = 5
_ifdep = 600.0
_gldep = 600.0
_wcdep = 200.0
_up_thk = _gldep
_up_vel = -1000.0


_nc_input_fname = '%s.in.nc' % _run_name

#############################################################

plume_nl_filename = '%s.nl' % _run_name
gc_config_filename = '%s.config' % _run_name


nc_regrid_args = [None,
                  '%s.in.nc'% _run_name, None,
                  _m,_n,0,0,0,0,0,0,0,0,
                  _kinbcw,0,0,0]

nc_gen_input_args = ['gs',
                     _nc_input_fname,
                     _m,_n,_kinbcw, _n_level,_hx,_hy,
                     _up_thk, _up_vel, _ifdep,
                     _kx, _chan_depth,_ramp_width,_landw,_ifpos
                    ]

input_style = 'gen_input'

checkpoint_duration = 0.10
checkpoint_tstart =  0.0
checkpoint_tend   =  1.0
checkpoint_initial_input_style = 'gen_input'

plume_vals = {'ifdep' : _ifdep, 
              'plume_min_thickness' : 10.0,
              'm_grid' : _m,
              'n_grid' : _n,
              'hx' : _hx,
              'hy' : _hy,
              'ah' : 0.0,
              'kh' : 0.0,
              'dt1' : 50.0,
              'phi' : 0.0,
              'rholinear' : 1,
              'tangle'   : 0,
              'tempbot' : 0.0,
              'temptop' : 0.0,
              'saltbot' : 34.5,
              'salttop' : 34.5,
              'gldep' : _gldep,
              'wcdep' : _wcdep}

gc_vals = {'plume' : { 'plume_imin' : 1,
 			'plume_imax' : _m,
			'plume_kmin' : 1,
	        	'plume_kmax' : _n-2,
                        'plume_steadiness_tol' : 1.0e-4,
			'plume_output_prefix' : _run_name,
			'plume_output_file' : 'plume.%s.out.nc' % _run_name,
	                'plume_nl_file' : plume_nl_filename,
                        'plume_const_bmlt' : 0,
                        'plume_const_bmlt_rate' : 0.0
                       },
           'time' : {'tstart' : _tstart,
                     'tend' : _tend,
 		      'dt' : _dt },
	    'grid' : { 'upn' : _n_level,
			'ewn' : _m,
			'nsn' : _n,
			'dew' : _hx,
			'dns' : _hy },
           'parameters' : { 'ice_limit' : 20.0,
                            'default_flwa' : 1.0e-16,
                            },
           'Petermann shelf' : {'accumulation_rate' : 0.0,
                                'air_temperature' : -5.0},
 	    'CF output' : { 'frequency' : _dt,
                            'start' : _tstart,
                            'stop' : _tend,
			    'name' : '%s.out.nc' % _run_name},
            'CF input' : { 'name' : _nc_input_fname,
	                   'time' : 1 },
	    'CF default' : { 'comment' : '',
			     'title' : _run_name },
           'options' : { 'flow_law' : 2,
                         'use_plume' : 0,
                         'hotstart' : 0,
                         'x_invariant' : 0,
                         },
           'ho_options': {'which_bmelt' : 0,
                          'which_disp' : 0,
                          },


	}


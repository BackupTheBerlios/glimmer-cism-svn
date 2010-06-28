# Sample plume-GC configuration file

# Note: Uses python syntax, so this file must be a valid python file

_run_name = 'steady_test'

_m = 5
_n = 81
_n_level = 3
_hx = 500.0
_hy = _hx

_dt = 0.0125
_tstart = 0.0
_tend =   2.0

_kinbcw = 2
_ramp_width = 0
_landw = 0
_plume_landw = 2
_kx = 0.0
_chan_depth = 0.0
_ifpos = 5
_gldep = 1000.0
_ifdep = _gldep 
_wcdep = 200.0
_up_thk = _gldep
_up_vel = -1000.0

_const_bmlt_rate = 1.0

_nc_input_fname = '%s.in.nc' % _run_name

# This is to do a restart
_restart_nc_file = ''
_restart_nc_file = ''
_restart_nc_file_t_read = -1

_use_plume = 1

#############################################################

plume_nl_filename = '%s.nl' % _run_name
gc_config_filename = '%s.config' % _run_name


nc_regrid_args = [_restart_nc_file,
                  '%s.in.nc'% _run_name, _restart_nc_file_t_read,
                  _m,_n,0,0,0,0,0,0,0,0,
                  _kinbcw,0,0,0]

nc_gen_input_args = ['gs',
                     _nc_input_fname,
                     _m,_n,_kinbcw, _n_level,_hx,_hy,
                     _up_thk, _up_vel, _ifdep,
                     _kx, _chan_depth,_ramp_width,_landw,_ifpos
                    ]

input_style = 'gen_input'
#input_style = 'regrid'

if (input_style == 'regrid'):
    _hotstart = 1
else:
    _hotstart = 0

checkpoint_duration = 5.0
checkpoint_tstart =  0.0
checkpoint_tend   =  25.0
checkpoint_initial_input_style = 'gen_input'
checkpoint_initial_input_style = 'regrid'

plume_vals = {'ifdep' : _ifdep, 
              'plume_min_thickness' : 10.0,
              'm_grid' : _m+2*_plume_landw,
              'n_grid' : _n+1*_plume_landw,
              'hx' : _hx,
              'hy' : _hy,
              'ah' : 100.0,
              'kh' : 100.0,
              'dt1' : 50.0,
              'phi' : 0.0,
              'tangle' : 0,
              'rholinear' : 1,
              'tempbot' : -1.850,
              'temptop' : -1.850,
              'saltbot' : 34.5,
              'salttop' : 34.5,
              'gldep' : _gldep,
              'wcdep' : _wcdep}

gc_vals = {'plume' : { 'plume_imin' : 1,
 			'plume_imax' : _m+2*_plume_landw,
			'plume_kmin' : _ifpos,
	        	'plume_kmax' : _n+1*_plume_landw,
                        'plume_steadiness_tol' : 1.0e-4,
			'plume_output_prefix' : _run_name,
			'plume_output_file' : 'plume.%s.out.nc' % _run_name,
	                'plume_nl_file' : plume_nl_filename,
                        'plume_const_bmlt' : 1,
                        'plume_const_bmlt_rate' : _const_bmlt_rate},
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
			'temperature' : 0,
                         'use_plume' :_use_plume,
                         'x_invariant' : 1,
                         'hotstart' :  _hotstart},
           'ho_options' :{'which_disp' : 0,
                          'which_ho_babc' : 6,   # no basal traction anywhere
                          },

	}


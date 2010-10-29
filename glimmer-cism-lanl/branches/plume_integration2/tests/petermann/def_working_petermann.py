from gcplume import *
import new

j = SteadyShelfJob()

j.name = 'working_petermann'
j.jobDir = '.'

j.default_flwa = 1.0e-16
j.uniform_acab = 0.0

j.m = int(40)
j.n = int(120)
j.nlevel = 3

j.hx = 500.0
j.hy = 500.0

j.tstart = 0.0
j.ice_dt = 0.01
j.tend = 100.0

j.use_plume = 1
j.plume_dt = 60.0 #seconds

j.otopg = -2000.0
j.upthk = 600.0
j.ifthk = 500.0

j.randthk = 0.00
j.kx = 2.0
j.chan_amp = 0.0
j.chan_init_length = 5000.0

j.upvel = -1000.0

# no-slip walls if 0, free-slip if 3

j.total_side_buf = 0
#j.total_side_buf_west = 0
#j.total_side_buf_east = 0

j.plume = { 'plume_min_thickness' : 2.0,
	    'entr_time_const' : 3600.0, # one hours
	    'use_min_plume_thickness' : True,
	    'knfloain' : j.n + 2 - 5 + 1,
	    'knfloein' : j.n + 2 - 2,
	    'infloain' : 3,
	    'infloein' : j.m+2,
	    'depinffix' : 2.0,
	    'entype' : 1,
	    'tangle' : True,
	    'horturb' : True,
	    'nonlin' : True,
	    'rholinear' : True,
            'cdb' : 1.5e-3 * 1.0,  #drag constant
	    'cl' : 1.775e-2 * 1.0,
            'temptop' : 0.0,
            'tempbot' : 0.0,
            'salttop' : 34.765,
            'saltbot' : 34.765,
            'tiuniform' : -20.0,
            'min_melt_depth' : 0.0,
            'phi'     : 80.0,
	    'ah' : 1000.0,
	    'kh' : 1000.0,
	    'plume_southern_bc' : 0,
            'snottim' : 86400.0 / 24.0, #daily
            'lnottim' : 86400.0 / 24.0, #daily
	   }

j.gc = {'options' : {'flow_law' : 0,
                     'temperature' : 0,
                     },

        'boundary condition params' : {'tau_xy_0' : 0.0*1000.0,
                                       'x_invariant' : False,
                                       'use_lateral_stress_bc' : True,
                                       },
        'Petermann shelf' : { 'air_temperature' : -20.0,
                              'accumulation_rate' : 0.0,
                              },
        
        'picard parameters' : {'small_vel' : 0.01,
                               'minres' : 1.0e-6,
                               'y_overrideres' : 1.0e-9,
                               'cvg_accel' : 1.25,
                               },
        'plume' : {'plume_const_bmlt' : False,
                   'plume_steadiness_tol' :
                          10.0e-2*j.plume_dt/(3600.0*24.0*365.25),
		   'plume_speed_steadiness_tol' :
                          10.0e-2*j.plume_dt/(3600.0*24.0*365.25),
                   'plume_write_all_states' : True,
		   'plume_write_every_n' : int(600.0/j.plume_dt),
                   'plume_max_spinup_time' : 50.0,
                   'plume_min_spinup_time' : 5.0e0,
		   'plume_min_subcycle_time' : 1.0,
                   'plume_imin' : 1,
                   'plume_imax' : j.m + 4,
                   'plume_kmin' : 5,
                   'plume_kmax' : j.n-1+2, #+2 is to include some rock
#                   'plume_kmax' : j.n-1, 
                   },
        }

j.assertCanStage()
j.serialize()
j.stage()

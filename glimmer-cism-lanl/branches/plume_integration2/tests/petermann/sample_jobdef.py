from gcplume import *

j = SteadyShelfJob()

j.name = 'sample'
j.jobDir = '.'
j.m = 20
j.n = 60
j.nlevel = 3
j.tstart = 0.0
j.tend =  1.0
j.ice_dt = 0.1
j.hx = 1000.0
j.hy = 1000.0
j.randthk = 0.00
j.plume_dt = 180.0
j.upthk = 600.0
j.otopg = -1200
j.upvel = -1000.0

j.default_flwa = 1.0e-16
j.uniform_acab = -1.0
j.use_plume = 1

j.plume['saltbot'] = 34.5
j.plume['salttop'] = 34.5
j.plume['temptop'] = -1.0
j.plume['tempbot'] = -1.0
j.plume['phi'] = 75.0
j.plume['plume_min_thickness'] = 20.0

j.gc = {'options': {'flow_law' : 0,
	            'temperature' : 0,
	           },
        'boundary condition params' : {'tau_xy_0' : 25.0e+3,
                                       'x_invariant' : False,
                                       'use_lateral_stress_bc' : True,
                                       'use_plastic_bnd_cond' : False,
                                       },
        'picard parameters' : {'small_vel' : 0.01,
                               'minres' : 1.0e-6,
                               'y_overrideres' : 1.0e-9,
                               'cvg_accel' : 1.25,
                               },
        'plume' : {'plume_const_bmlt' : False,
		   'plume_steadiness_tol' : 1.0e-4,
                   },
        }


j.assertCanStage()
j.serialize()


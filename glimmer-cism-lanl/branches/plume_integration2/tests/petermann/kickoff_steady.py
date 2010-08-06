
from gcplume import *
import os

for tau_xy_0 in [0,10, 25, 50]:
    for acab in [-2.0, 0.0, 2.0]:

        j = SteadyShelfJob()
        j.name = 'ssj1_%s_kPa_%s_acab' % (tau_xy_0, acab)
        jd = os.path.expandvars('$GC_JOBS/%s' % j.name)
        if (not(os.path.lexists(jd))):
            os.mkdir(jd)
        j.jobDir = jd

        # varying properties
        j.uniform_acab = acab

        # fixed default properties
        j.m = 20
        j.n = 40
        j.hx = 1000.0
        j.hy = 1000.0
        j.nlevel = 3

        j.default_flwa = 1.0e-16
        j.upvel = -1000.0
        j.upthk = 600.0
        j.randthk = 0.0
        j.otopg = -2000.0

        j.tend = 200.0
        j.tstart = 0.0
        j.ice_dt = 0.1

        j.use_plume = 0
        j.plume_dt = 0.0
        
        j.gc = {'options' : { 'flow_law' : 2,
                              },
                'Petermann shelf' : { 'accumulation_rate' : acab,
             		              'thk_steady_tol' : 1.0e-3*j.ice_dt,
                                     },
                'boundary condition params' : {'tau_xy_0' : tau_xy_0 * 1000.0,
                                               },
	        'picard parameters' : {'small_vel' : 0.0001,
	 			       'minres' : 1.0e-6,
				       'y_overrideres' : 1.0e-9,
					'cvg_accel' : 1.25,
		                      },
                                               
                }
        j.assertCanStage()
        j.serialize()
        j.stage()
        j.run()
    
    

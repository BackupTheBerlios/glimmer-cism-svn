from submit_GC_job import submit_job
from gcplume import *

import os
import sys

def kickoff(unique_str, queue_mode, email, walltime):

    
    for tau_xy_0 in [0,10, 25, 50]:
        for acab in [-2.0, 0.0, 2.0]:
            
            j = SteadyShelfJob()
            j.name = 'ssj1_%s_kPa_%s_acab_%s' % (tau_xy_0, acab,unique_str)
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

            j.tend = 400.0
            j.tstart = 0.0
            j.ice_dt = 0.1

            j.use_plume = 0
            
            j.plume_dt = 0.0

            j.plume = {'salttop' : 0.0,
                       'saltbot' : 0.0,
                       'temptop' : 0.0,
                       'tempbot' : 0.0,
                       'phi' : 0.0}
            
            j.gc = {'options' : { 'flow_law' : 2,
                                  'temperature' : 0,
                                  },
                    'Petermann shelf' : { 'accumulation_rate' : acab,
                                          'thk_steady_tol' : 1.0e-3*j.ice_dt,
                                          },
                    'boundary condition params' : {'tau_xy_0' : tau_xy_0 * 1000.0,
                                                   },
                    'picard parameters' : {'small_vel' : 0.01,
                                           'minres' : 1.0e-6,
                                           'y_overrideres' : 1.0e-9,
                                           'cvg_accel' : 1.25,
                                           },
#                    'plume' : {'plume_const_bmlt' : True,
#                               'plume_const_bmlt_rate' : -1.0 * acab,
#                               },
                    }
            j.assertCanStage()
            j.serialize()
            print( "submitting job %s" % j.name)
            submit_job(j,email,walltime,queue_mode)
    
    

USAGE = 'python kickoff_steady.py <unique_str> <queue_mode>'

if __name__ == '__main__':

    if (len(sys.argv) != 3):
        print (USAGE)

    else:
        kickoff( sys.argv[1], sys.argv[2],'gladish@cims.nyu.edu','24:00:00')
    

    

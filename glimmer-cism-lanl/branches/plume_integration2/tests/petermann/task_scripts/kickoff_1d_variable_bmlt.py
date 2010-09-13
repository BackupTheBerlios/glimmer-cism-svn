import os

from gcplume import *
from submit_GC_job import *

def kickoff( unique_str,queue_mode, email, walltime ):

    j = SteadyShelfJob()
    
    j.default_flwa = 1.0e-16
    j.n = 81
    j.m = 5
    j.nlevel = 3
    j.tend = 200.0
    j.tstart = 0.0
    j.ice_dt = 0.025
    j.hx = 500.0
    j.hy = 500.0
    j.use_plume = 1 # but with const melt rate

    j.plume['salttop'] = 0.0
    j.plume['saltbot'] = 0.0
    j.plume['temptop'] = 0.0
    j.plume['tempbot'] = 0.0
    j.plume['phi']  = 0.0

    j.plume_dt = 0.0
    
    j.otopg = -2000.0
    j.upthk = 1000.0
    j.randthk = 0.0
    j.upvel = -1000.0
    
    j.gc = {'options' : {'flow_law' : 2,     # use default value of A
                         'temperature' : 0,  # set to air temp
                         },
            'boundary condition params' : {'x_invariant' : True,
                                           'use_lateral_stress_bc' : True,
                                           'tau_xy_0' : 0.0,
                                           'use_plastic_bnd_cond' : False,
                                           },
            'picard parameters' : {'small_vel' : 0.01,
                                   'minres' : 1.0e-6,
                                   'y_overrideres' : 1.0e-9,
                                   'cvg_accel' : 1.25,
                                   },
            'plume' : { 'salttop' : 0.0,
                        'saltbot' : 0.0,
                        'temptop' : 0.0,
                        'tempbot' : 0.0,
                        },
            }


    acabs = [0, 10, 25]

    for acab in acabs:

        j.uniform_acab = acab * -1.0
        j.name = '1d_%s_bmlt_%s' % (acab, unique_str)

        j.gc.update( {'plume' : {'plume_const_bmlt' : True,
                                 'plume_const_bmlt_rate' : acab*1.0,
                     
                                 },
                      })

        j.assertCanStage()

        jdir = os.path.join(os.path.expandvars('$GC_JOBS'),
                            j.name)

        if (not(os.path.lexists(jdir))):
            os.mkdir(jdir)
        
        j.jobDir = jdir
        
        j.serialize()
        
        submit_job(j,email, walltime, queue_mode)
                    

USAGE = 'python kickoff_1d_variable_bmlt.py <unique_str> <queue_mode>'

if __name__ == '__main__':

    if (len(sys.argv) != 3 ):
        raise Exception("Call like: \n %s" % USAGE)

    unique_str = sys.argv[1]
    queue_mode = sys.argv[2]
    
    kickoff(unique_str, queue_mode,  'gladish@cims.nyu.edu', '48:00:00')


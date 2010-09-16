import os

from gcplume import *
from submit_GC_job import *

def kickoff( unique_str,queue_mode, email, walltime ):


    thks = [100, 200, 500]
    upvels = [-1000.0]

    for thk in thks:
        for vel in upvels:

            j = LinearShelfJob()

            j.default_flwa = 1.0e-16
            j.n = 100
            j.m = 5

            j.nlevel = 3
            j.tend = 0.001
            j.tstart = 0.0

            j.ice_dt = 0.001

            j.hx = 40000 / j.n
            j.hy = j.hx

            j.use_plume = 0 # but with const melt rate
            j.plume_dt = 0.0

            j.otopg = -2000.0
            j.ifthk = thk * 1.0
            j.upthk = thk * 1.0
            j.randthk = 0.0

            j.upvel = vel

            j.gc = {'options' : {'flow_law' : 2,     # use default value of A
                                 'temperature' : 0,  # set to air temp
                                 },
                    'boundary condition params' : {'x_invariant' : True,
                                                   'use_lateral_stress_bc' : True,
                                                   'tau_xy_0' : 0.0,
                                                   'use_plastic_bnd_cond' : False,
                                                   },
                    'picard parameters' : {'small_vel' : 0.00,
                                           'minres' : 1.0e-9,
                                           'y_overrideres' : 1.0e-9,
                                           },
                    }

            j.uniform_acab = 0.0
            j.name = 'flat_shelf_%s_m_%s' % (thk, unique_str)

            j.assertCanStage()

            jdir = os.path.join(os.path.expandvars('$GC_JOBS'),
                                j.name)


            if (not(os.path.lexists(jdir))):
                os.mkdir(jdir)

            j.jobDir = jdir

            j.serialize()

            submit_job(j,email, walltime, queue_mode)
                    

USAGE = 'python kickoff_flat_shelf.py <unique_str> <queue_mode>'

if __name__ == '__main__':

    if (len(sys.argv) != 3 ):
        raise Exception("Call like: \n %s" % USAGE)

    unique_str = sys.argv[1]
    queue_mode = sys.argv[2]
    
    kickoff(unique_str, queue_mode,  'gladish@cims.nyu.edu', '48:00:00')


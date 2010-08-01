import os

from gcplume import *
from submit_GC_job import *

def kickoff( start_index, email, walltime ):

    j = SteadyShelfJob()
    j.default_flwa = 1.0e-16
    j.uniform_acab = -1.2
    j.n = 70
    j.m = 20
    j.nlevel = 3
    j.tend = 200.0
    j.tstart = 0.0
    j.ice_dt = 0.1
    j.hx = 1000.0
    j.hy = 1000.0
    j.use_plume = 1
    j.plume_dt = 60.0
    j.otopg = -2000.0
    j.upthk = 600.0
    j.randthk = 0.0
    j.gc = {'options' : {'flow_law' : 2,
                         'temperature' : 0,
                         'x_invariant' : False,
                         'use_lateral_stress_bc' : True,
                         },
            'parameters' : { 'tau_xy_0' : 50.0e+3,
                             },
            'plume' : {'plume_const_bmlt' : False,
                       'plume_steadiness_tol' : 1.0e-4,
                       },
            }


    jobindex = start_index

    oceantemps = [-1.0, 0.0, 1.0]
    upvels = [-500.0, -1000.0, -1500.0]
    phis = [0.0, 75.0]

    for t in oceantemps:
        for upvel in upvels:
            for phi in phis:

                j.name = 'job%s' % jobindex
                
                jdir = os.path.join(os.path.expandvars('$GC_JOBS'),
                                    j.name)

                if (not(os.path.lexists(jdir))):
                    os.mkdir(jdir)
        
                j.jobDir = jdir
            
                j.upvel = upvel
                j.plume = {'temptop' : t,
                           'tempbot' : t,
                           'salttop' : 34.5,
                           'saltbot' : 34.5,
                           'plume_min_thickness' : 20.0,
                           'phi'     : phi}

                j.assertCanStage()
                j.serialize()
                submit_job(j,email, walltime,'q')

                jobindex = jobindex+1


if __name__ == '__main__':

    start_index = int(sys.argv[1])
    
    
    kickoff(start_index, 'gladish@cims.nyu.edu', '12:00:00')


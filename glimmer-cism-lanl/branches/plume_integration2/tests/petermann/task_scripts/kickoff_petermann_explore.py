#!/usr/bin/python

import os

from gcplume import *
from submit_GC_job import *

def kickoff( email, walltime, unique_str, queue_mode ):


    oceantemps = [-0.5, 0.0, 0.5]
    upvels = [-900.0, -1000.0, -1100.0]
    phis = [0.0,81.0]
    min_melt_depths = [0.0, 100.0, 200.0]
    cdbs = [1,2,4,5,10]
    taus = [0,10,25,50]

    upvel = -1000.0
    t = 0.0
    #for t in oceantemps:
    for tau in taus:
        #for upvel in upvels:
        for cdb in cdbs:
            for phi in phis:
                for min_melt_depth in min_melt_depths:
                    j = LinearShelfJob()

                    j.default_flwa = 1.0e-16

                    j.uniform_acab = -1.2  #going into nc_gen_input only

                    j.n = 50
                    j.m = 20
                    j.nlevel = 3
                
                    j.hx = 1000.0
                    j.hy = 1000.0
                
                    j.tend = 300.0
                    j.tstart = 0.0
                    j.ice_dt = 0.0125

                    j.use_plume = 1
                    j.plume_dt = 30.0

                    j.otopg = -2000.0
                    j.upthk = 600.0
                    j.ifthk = 550.0
                    j.randthk = 0.0

                    j.upvel = upvel

                    j.plume = { 'plume_min_thickness' : 1.0,
                                'cdb' : 1.5e-3 * (1.0/cdb),  #drag constant
                                'temptop' : t,
                                'tempbot' : t,
                                'salttop' : 34.765,
                                'saltbot' : 34.765,
                                'tiuniform' : -20.0,
                                'min_melt_depth' : min_melt_depth,
                                'phi'     : phi }

                    j.gc = {'options' : {'flow_law' : 0,
                                         'temperature' : 0,
                                         },
                            'boundary condition params' : {'tau_xy_0' : tau*1000.0,
                                                           'x_invariant' : False,
                                                           'use_lateral_stress_bc' : True,
                                                           },
                            'Petermann shelf' : { 'air_temperature' : -20.0,
                                                  'accumulation_rate' : -1.2,
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

                
                    j.name = 'pn_%.0fd_%.0fm_%s_cdb_%s_kPa_%s' % (phi,min_melt_depth,cdb,tau,unique_str)
                    jdir = os.path.join(os.path.expandvars('$GC_JOBS'),
                                        j.name)
                    
                    if (not(os.path.lexists(jdir))):
                        os.mkdir(jdir)
        
                    j.jobDir = jdir
            
                    j.assertCanStage()
                    j.serialize()
                    submit_job(j,email, walltime,queue_mode)
                    

USAGE = 'python kickoff_plastic_central.py <unique_str> <queue_mode>'

if __name__ == '__main__':

    if (len(sys.argv) != 3 ):
        raise Exception("Call like: \n %s" % USAGE)
    unique_str = sys.argv[1]
    queue_mode = sys.argv[2]
    
    kickoff('gladish@cims.nyu.edu', '48:00:00',unique_str, queue_mode)


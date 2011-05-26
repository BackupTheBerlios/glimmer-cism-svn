from submit_GC_job import *
from gcplume import *


f = open('jobs_may21.T-0.5_coupled_unexplained.txt','r')

for l in f:

    (t,ice_t, tauxy0, diff, k, amp, code) = [n.strip() for n in l.split(',')]


    jobname = 'may21.T%s.2_coupled_temp_%s_tau_%s_diff_%s_k_%s_amp_%s' % (t,ice_t, tauxy0, diff, k, amp)
    cold_jobname = 'may21.T%s_coupled_temp_%s_tau_%s_diff_%s_k_%s_amp_%s' % ('-1.0',ice_t, tauxy0, diff, k, amp)

    j = RestartIceJob(cold_jobname,newName=jobname)

    j.tend = 0.0
    j.useMaxRunTimeLimit = True
    j.maxRunTimeLimit = 200.0

    j.doPlumeRestart = False

    j.gc['picard parameters']['cmax'] = 250
    j.gc['picard parameters']['start_umc'] = 300
    
    j.plume['tempbot'] = -0.5

    j.assertCanStage()
    j.serialize()
    print (cold_jobname)
    submit_job(j,'carlgladish@gmail.com','8:00:00','q',1,True)

    
f.close()

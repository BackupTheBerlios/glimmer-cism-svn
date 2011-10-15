#!/usr/bin/python

#
# run Glimmer-cism job
#

import sys
import os.path
import pickle
import subprocess


def readjob(jobfile):

    if (not(os.path.exists(jobfile))):
        raise Exception('ERROR: job file not found')
    g = open(jobfile,'r')
    try:
        j = pickle.load(g)
    finally:
        g.close()
        
    return j


def write_jobscript(j,email,jobfile,walltime,ppn,overwrite):

    #
    # make scripts to run the code
    #

    f = open(os.path.join(os.path.expandvars('$GC_JOB_SCRIPTS'),
                          'run_GC_%s.sh' % j.name),
             'w')

    script = '''#!/bin/sh
                #PBS -l nodes=1:ppn=%s,walltime=%s
                #PBS -N %s
                #PBS -M %s
                #PBS -m abe
                #PBS -e localhost:%s/\${PBS_JOBNAME}.e
                #PBS -o localhost:%s/\${PBS_JOBNAME}.o

                p=$PWD
                cd %s
                export OMP_NUM_THREADS=%s
                python `which run_job.py` %s %s &> %s.out.log
                cd $p
                exit 0;
                EOF
            ''' % (ppn,walltime, j.name, email,
                   os.path.expandvars('$GC_JOB_SCRIPTS/'),
                   os.path.expandvars('$GC_JOB_SCRIPTS/'),
                   os.path.dirname(os.path.abspath(jobfile)),
		   ppn,
                   os.path.abspath(jobfile),
                   (overwrite and '-overwrite' or ''),
                   j.name)

 
    f.write(script)
    f.close()

def queue_job(mode, j):
    #
    # if queueing job, qsub the script
    #

    runcmd =  os.path.join(os.path.expandvars('$GC_JOB_SCRIPTS'),
                           'run_GC_%s.sh' % j.name)

    if (mode != 'q'):
        print("running job script  %s" % runcmd)
        if (mode == 'fake'):
            pass
        else:
            subprocess.call(['sh', runcmd])
    else:
        print('queueing jobscript: %s' % runcmd)
        retcode = subprocess.call(['qsub', runcmd])
	if (retcode != 0):
	    raise Exception("Error running qsub")

def submit_job(job,email,walltime,mode,ppn,overwrite):
    
    if (type(job) is str):
        j = readjob(job)
    else:
        j = job

    j.assertCanStage()
    write_jobscript(j,email,j.serialFile,walltime,ppn,overwrite)
    
    queue_job(mode,j)

    
def main():

    USAGE="submit_GC_job.sh <job_file> <walltime> <mode> [<ppn>] [-overwrite]"

    #
    # set user-specific stuff
    #

    EMAIL='gladish@cims.nyu.edu'
    ppn = 4
    OVERWRITE = False
    
    #
    # trap incorrect number of arguments
    #

    if (len(sys.argv) < 4 ):
        print('ERROR: incorrect number of arguments')
        print(USAGE)
        return

    #
    # set variables
    #

    JOB=sys.argv[1]
    print('%s: %s' % ('jobfile', JOB))
    WALLTIME=sys.argv[2]
    print('%s: %s' % ('walltime',WALLTIME))
    MODE=sys.argv[3]
    print('%s: %s' % ('mode', MODE))
    if (len(sys.argv) >= 5):
        arg5 = sys.argv[4]
        if (arg5.startswith('-')):
            if (arg5 == '-overwrite'):
                OVERWRITE = True
            else:
                raise Exception("Unrecognized option: %s" % arg5)
        else:
            ppn = sys.argv[4]

    if (len(sys.argv) == 6):
        if (sys.argv[5] == '-overwrite'):
            OVERWRITE = True
        else:
            raise Exception("Unrecognized option: %s" % sys.argv[5])
        
    print('%s: %s' % ('ppn', ppn))
    print('%s: %s' % ('overwrite', OVERWRITE))
    
    submit_job(JOB,EMAIL,WALLTIME,MODE,ppn,OVERWRITE)


if __name__ == '__main__':
    main()

    

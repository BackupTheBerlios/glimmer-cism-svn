#!/usr/bin/python

#
# run Glimmer-cism job
#

import sys
import os.path
import pickle
import subprocess


def readjob(jobfile):
    #jobfile should be a .gcpl file 
    j = None
    if (not(os.path.exists(jobfile))):
        print('WARNING: job file not found')
    else:
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
    #job should be the name of a .gcpl file, or else an actually job object
    if (type(job) is str):
        j = readjob(job)
    else:
        j = job
    j.assertCanStage()
    write_jobscript(j,email,j.serialFile,walltime,ppn,overwrite)
    queue_job(mode,j)

#def submit_job(job,email,walltime,mode,ppn,overwrite):
#    print "submitting %s" % str(job)
    
def glob_jobs(joblistfile):
    #joblistfile is the name of a textfile
    joblist = []

    if os.path.exists(joblistfile):
        f = open(joblistfile)
        names = f.readlines()
        f.close()
        names = [os.path.join(os.path.expandvars('$GC_JOBS'),
                              n.strip(),'%s.gcpl' % n.strip()) for n in names]
        if (reduce((lambda x,y:x and y),[os.path.exists(n) for n in names],True)):
            joblist = names

    return joblist
    
def main(job_files=None):

    #
    # trap incorrect number of arguments
    #
    
    USAGE="submit_GC_job.sh <job_file> <walltime> <mode> [<ppn>] [-overwrite]"
    if (len(sys.argv) < 4 ):
        print('ERROR: incorrect number of arguments')
        print(USAGE)
        return

    #
    # set variables
    #
    
    EMAIL='gladish@cims.nyu.edu'
    ppn = 4
    OVERWRITE = False
    
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

    gjobs = glob_jobs(JOB)
    if gjobs:
        for j in gjobs:
            submit_job(j,EMAIL,WALLTIME,MODE,ppn,OVERWRITE)
    else:
        submit_job(JOB,EMAIL,WALLTIME,MODE,ppn,OVERWRITE)


if __name__ == '__main__':
    main()

    

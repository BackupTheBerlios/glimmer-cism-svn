#!/usr/bin/python

#
# run Glimmer-cism job
#

import sys
import os.path
import subprocess

def write_jobscript(jname,jobs_dir,new_tstep,email,walltime):

    #
    # make scripts to run the code
    #

    scriptname = os.path.join(os.path.expandvars('$GC_JOB_SCRIPTS'),
                              'run_filter_%s.sh' % jname)
    f = open(scriptname,'w')

    script = '''#PBS -l nodes=1:ppn=1,walltime=%s
#PBS -N filter_%s
#PBS -M %s
#PBS -m abe
#PBS -e localhost:%s/\${PBS_JOBNAME}.e
#PBS -o localhost:%s/\${PBS_JOBNAME}.o
    
p=$PWD
cd %s
python `which filter_times.py` %s `ls -I*log` &> %s.out.log
cd $p
exit 0;
EOF''' % (walltime, jname, email,
          os.path.expandvars('$GC_JOB_SCRIPTS/'),
          os.path.expandvars('$GC_JOB_SCRIPTS/'),
          jobs_dir,new_tstep,
          jname)

 
    f.write(script)
    f.close()

    return scriptname

def queue_job(mode, scriptname):
    #
    # if queueing job, qsub the script
    #

    if (mode != 'q'):
        print("running job script  %s" % scriptname)
        if (mode == 'fake'):
            pass
        else:
            subprocess.call(['sh', scriptname])
    else:
        print('queueing jobscript: %s' % scriptname)
        retcode = subprocess.call(['qsub', scriptname])
	if (retcode != 0):
	    raise Exception("Error running qsub")

def submit_job(jname,jobs_dir,new_tstep,email,walltime,mode):
    
    script = write_jobscript(jname,jobs_dir,new_tstep,email,walltime)
    
    queue_job(mode,script)

    
def main():

    USAGE="submit_filter_job.sh <jobs_dir> <new_tstep> <walltime> <mode>"

    #
    # set user-specific stuff
    #

    EMAIL='gladish@cims.nyu.edu'
    ppn = 1
    
    #
    # trap incorrect number of arguments
    #

    if (len(sys.argv) < 5 ):
        print('ERROR: incorrect number of arguments')
        print(USAGE)
        return

    #
    # set variables
    #

    JOBSDIR=sys.argv[1]
    print('%s: %s' % ('jobs_dir', JOBSDIR))
    NEWTSTEP=float(sys.argv[2])
    print('%s: %s' % ('new_tstep', NEWTSTEP))    
    WALLTIME=sys.argv[3]
    print('%s: %s' % ('walltime',WALLTIME))
    MODE=sys.argv[4]
    print('%s: %s' % ('mode', MODE))
    
    submit_job(os.path.basename(os.path.abspath(JOBSDIR)),
               os.path.abspath(JOBSDIR),
               NEWTSTEP,
               EMAIL,
               WALLTIME,
               MODE)


if __name__ == '__main__':
    main()

    

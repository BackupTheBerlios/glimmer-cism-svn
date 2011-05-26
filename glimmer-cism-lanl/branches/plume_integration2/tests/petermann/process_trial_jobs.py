#!/usr/bin/env python

import os
import re
import gcplume
import pickle
import subprocess
import sys
from submit_GC_job import *

def main(jobsPrefix,doRestart,pattern=None):
    
    UNEXPLAINED = 0
    PLUME_ERROR = 1
    ICE_ERROR = 2
    STEADY_STATE = 3
    TIMED_OUT = 4
    RESTART = 5

    jobs_dir  = os.path.join(os.path.expandvars('$GC_JOBS'))
    all_jobs = os.path.os.listdir(jobs_dir)
    jobs = [jdir for jdir in all_jobs if jdir.startswith(jobsPrefix)]

    steady_ice_list = []
    plume_error_list = []
    ice_error_list = []
    killed_list = []
    unexplained_list = []
    restart_list = []
    restart_jobs = []
    job_list = []

    pattern = '.*temp_([\-0-9\.]*)_tau_([0-9\.]*)_diff_([0-9\.]*)_trial_([0-9]*)' 
    
    for jname in jobs:
        unexplained = True
        g = re.match(pattern,jname).group

        job_key = (float(g(1)),float(g(2)),float(g(3)),int(g(4)))
        print(job_key)
        job_list.append(job_key)

        job_file = open(os.path.join(os.path.expandvars('$GC_JOBS'),
                                     jname,'%s.gcpl' % jname),'r')
        try:
            j_obj = pickle.load(job_file)
        finally:
            job_file.close()

        # take a look at the ice log file
        log_file = os.path.join(os.path.expandvars('$GC_JOBS'),
                                jname,'%s.config.log' % jname)

        cmd = ['tail','-n20',log_file]
        log_tail_p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
        log_tail_p.wait()
        final_log_lines = log_tail_p.stdout.readlines()

        for l in final_log_lines:
            if l.strip().startswith('Stopped time-stepping'):
                unexplained = False
                steady_ice_list.append(job_key)

            if ('ERROR' in l):
                unexplained = False
                ice_error_list.append(job_key)
                restart_jobs.append(j_obj)

        job_contents_p = subprocess.Popen(['ls',os.path.join(os.path.expandvars('$GC_JOBS'),jname)],
                                          stdout=subprocess.PIPE)
        job_contents_p.wait()
        job_contents = job_contents_p.stdout.readlines()
        plume_output_file = ''
        for job_file in job_contents:
            if job_file.strip().endswith('_output'):
                plume_output_file = os.path.join(os.path.expandvars('$GC_JOBS'),
                                                 jname,
                                                 job_file.strip())

        if plume_output_file:
            cmd = ['tail','-n100',plume_output_file]
            plume_output_tail_p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
            plume_output_tail_p.wait()
            plume_output_final_lines = plume_output_tail_p.stdout.readlines()

            for l in plume_output_final_lines:
                if 'error' in l:
                    unexplained = False
                    plume_error_list.append(job_key)
                    restart_jobs.append(j_obj)
        else:
            print('Warning: did not find _output file for %s' % jname)

        # take a look at the job output log file to see if job was killed
        outlog_file = os.path.join(os.path.expandvars('$GC_JOBS'),
                                   jname,'%s.out.log' % jname)

        cmd = ['tail','-n100',outlog_file]
        outlog_tail_p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
        outlog_tail_p.wait()
        final_outlog_lines = outlog_tail_p.stdout.readlines()

        for l in final_outlog_lines:
            if 'process killed (SIGTERM)' in l:
                unexplained = False
                killed_list.append(job_key)
                restart_jobs.append(j_obj)
                break

            if ('Error running:' in l and doRestart):
                if(unexplained):
                    unexplained = False
                    restart_list.append(job_key)
                    job_gcpl_file = os.path.join(os.path.expandvars('$GC_JOBS'),
                                                 jname,'%s.gcpl' % jname)
                    restart_jobs.append(j_obj)
                    break
                
        if unexplained:
            unexplained_list.append(job_key)

    print('%s in all jobs' % len(jobs))        
    print('%s in steady_ice_list' % len(steady_ice_list))
    print('%s in ice_error_list' % len(ice_error_list))
    print('%s in plume_error_list' % len(plume_error_list))
    print('%s in killed_list' % len(killed_list))
    print('%s in unexplained_list' % len(unexplained_list))
    print('%s in restart_list' % len(restart_list))
    
    #print('intersections')
    #print(len(set(steady_ice_list).intersection(set(ice_error_list))))
    #print(len(set(steady_ice_list).intersection(set(plume_error_list))))
    #print(len(set(steady_ice_list).intersection(set(unexplained_list))))
    #print(len(set(steady_ice_list).intersection(set(killed_list))))
    #print(len(set(ice_error_list).intersection(set(plume_error_list))))
    #print(len(set(ice_error_list).intersection(set(unexplained_list))))
    #print(len(set(ice_error_list).intersection(set(killed_list))))
    #print(len(set(plume_error_list).intersection(set(unexplained_list))))
    #print(len(set(plume_error_list).intersection(set(killed_list))))
    #print(len(set(killed_list).intersection(set(unexplained_list))))

    #print('in restart and killed %s' % len(set(restart_list).intersection(set(killed_list))))
    #print('in restart and steady %s' % len(set(restart_list).intersection(set(steady_ice_list))))
    #print(len(set(restart_list).intersection(set(ice_error_list))))
    #print(len(set(restart_list).intersection(set(plume_error_list))))
    #print(len(set(restart_list).intersection(set(unexplained_list))))

    steady = open('jobs_%s_steady.txt' % (jobsPrefix),'w')
    iceerror=open('jobs_%s_iceerror.txt' % (jobsPrefix),'w')
    perror = open('jobs_%s_perror.txt' % (jobsPrefix),'w')
    killed = open('jobs_%s_killed.txt' % (jobsPrefix),'w')
    unex = open('jobs_%s_unexplained.txt' % (jobsPrefix),'w')
    restart = open('jobs_%s_restart.txt' % (jobsPrefix),'w')
    
    try:
        for job_key in job_list:
            if job_key in steady_ice_list:
                steady.write('%8.1f,%8.1f,%8.1f,%8.1f\n' %
                             (job_key[0],job_key[1],job_key[2],job_key[3]))
            elif job_key in plume_error_list:
                perror.write('%8.1f,%8.1f,%8.1f,%8.1f\n' %
                             (job_key[0],job_key[1],job_key[2],job_key[3]))
            elif job_key in ice_error_list:
                iceerror.write('%8.1f,%8.1f,%8.1f,%8.1f\n' %
                               (job_key[0],job_key[1],job_key[2],job_key[3]))
            elif job_key in killed_list:
                killed.write('%8.1f,%8.1f,%8.1f,%8.1f\n' %
                             (job_key[0],job_key[1],job_key[2],job_key[3]))
            elif job_key in unexplained_list:
                unex.write('%8.1f,%8.1f,%8.1f,%8.1f\n' %
                           (job_key[0],job_key[1],job_key[2],job_key[3]))
            elif job_key in restart_list:
                restart.write('%8.1f,%8.1f,%8.1f,%8.1f\n' %
                           (job_key[0],job_key[1],job_key[2],job_key[3]))
            else:
                print('not in any list: %8.1f,%8.1f,%8.1f,%8.1f\n' %
                      (job_key[0],job_key[1],job_key[2],job_key[3]))

    finally:
        steady.close()
        iceerror.close()
        perror.close()
        killed.close()
        unex.close()

    if (doRestart):
        #for job_gcpl_file in restart_jobs:
        for job_obj in restart_jobs:
            #print('restarting job: %s' % job_gcpl_file)
            #print('restarting %s' % str(job_obj))
            job_obj.gc['picard params']['cmax'] = 200
            submit_job(job_obj,'carlgladish@gmail.com','8:00:00','q',1,True)

if __name__ == '__main__':

    if (len(sys.argv) < 3):
        print ("usage: process_jobs.py <jobprefix> <doRestarts>")
    else:
        doRestart = sys.argv[2].startswith('t') or sys.argv[2].startswith('T')

        main(sys.argv[1],doRestart)

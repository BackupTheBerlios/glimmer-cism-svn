#!/usr/bin/env python

import os
import re
import gcplume
import pickle
import subprocess
import sys
from submit_GC_job import *

def main(jobsPrefix):
    
    jobs_dir  = os.path.join(os.path.expandvars('$GC_JOBS'))
    all_jobs = os.path.os.listdir(jobs_dir)
    jobs = [jdir for jdir in all_jobs if jdir.startswith(jobsPrefix)]

    job_list = []
    
    steady_ice_list = []
    plume_error_list = []
    ice_error_list = []
    job_error_list = []

    not_started_list = []
    still_running_list = []
    killed_list = []
    
    unexplained_list = []

    kv_pattern = '_([a-z,A-Z]+)_([-\.0-9]+)'
    
    for jname in jobs:
        
        unexplained = True
        stillRunning = True
        
        kv_list = re.findall(kv_pattern, jname[len(jobsPrefix):len(jname)])
        job_key = tuple([v for (k,v) in kv_list])
        print(jname)
        
        job_list.append(jname)

        #try:
        #    job_file = open(os.path.join(os.path.expandvars('$GC_JOBS'),
        #                                 jname,'%s.gcpl' % jname),'r')
        #    j_obj = pickle.load(job_file)
        #    job_file.close()
        #except:
        #    print('Could not load job_file %s' % job_file)
        #    return


        # list the files in the job directory
        job_contents_p = subprocess.Popen(['ls',os.path.join(os.path.expandvars('$GC_JOBS'),jname)],
                                          stdout=subprocess.PIPE)
        [stdout,stderr] = job_contents_p.communicate(None)
        job_contents = stdout.split('\n')

        job_log_file = ''
        plume_output_file = ''
        ice_log_file = ''
        
        for job_file in job_contents:
            if job_file.strip().endswith('.out.log'):
                job_log_file = os.path.join(os.path.expandvars('$GC_JOBS'),jname,job_file.strip())
            elif job_file.strip().endswith('_output'):
                plume_output_file = os.path.join(os.path.expandvars('$GC_JOBS'),jname,job_file.strip())
            elif job_file.strip().endswith('.config.log'):
                ice_log_file = os.path.join(os.path.expandvars('$GC_JOBS'),jname,job_file.strip())

        # take a look at the ice log file
        if (ice_log_file):
            cmd = ['tail','-n20',ice_log_file]
            log_tail_p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
            [stdout,stderr] = log_tail_p.communicate(None)
            final_log_lines = stdout.split('\n')

            for l in final_log_lines:
                if l.strip().startswith('Stopped time-stepping'):
                    unexplained = False
                    stillRunning = False
                    steady_ice_list.append(jname)

                if ('ERROR' in l):
                    unexplained = False
                    stillRunning = False
                    ice_error_list.append(jname)
        else:
            not_started_list.append(jname)
            unexplained = False
            stillRunning = False

        if plume_output_file:
            cmd = ['tail','-n100',plume_output_file]
            plume_output_tail_p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
            [stdout,stderr] = plume_output_tail_p.communicate(None)
            plume_output_final_lines = stdout.split('\n')

            for l in plume_output_final_lines:
                if 'error' in l:
                    unexplained = False
                    stillRunning = False
                    plume_error_list.append(jname)
                    break

        else:
            if (jname in not_started_list):
                pass
            else:
                print('Warning: did not find _output file for %s' % jname)
                

        # take a look at the job output log file to see if job was killed
        if (job_log_file):

            cmd = ['tail','-n100',job_log_file]
            outlog_tail_p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
            [stdout,stderr] = outlog_tail_p.communicate(None)
            final_outlog_lines = stdout.split('\n')

            for l in final_outlog_lines:
                if 'process killed (SIGTERM)' in l:
                    unexplained = False
                    stillRunning = False
                    killed_list.append(jname)
                    break

                if ('Error running:' in l):
                    stillRunning = False
                    if not((jname in ice_error_list) or (jname in plume_error_list)):
                        unexplained = True
                    break
        if (stillRunning):
            unexplained = False
            still_running_list.append(jname)                
                
        if (unexplained):
            unexplained_list.append(jname)
    
   
    print('%s in all jobs' % len(jobs))
    
    print('%s in steady_ice_list' % len(steady_ice_list))
    print('%s in ice_error_list' % len(ice_error_list))
    print('%s in plume_error_list' % len(plume_error_list))
    print('%s in killed_list' % len(killed_list))
    print('%s in still_running_list' % len(still_running_list))
    print('%s in not_started_list' % len(not_started_list))
    print('%s in unexplained_list' % len(unexplained_list))
    #print ('unexplained jobs: %s' % unexplained_list)
    #print ('ice error jobs: %s' % ice_error_list)
    #print ('plume error jobs: %s' % plume_error_list)
    steady = open('jobs_%s_steady.txt' % (jobsPrefix),'w')
    iceerror=open('jobs_%s_iceerror.txt' % (jobsPrefix),'w')
    perror = open('jobs_%s_perror.txt' % (jobsPrefix),'w')
    killed = open('jobs_%s_killed.txt' % (jobsPrefix),'w')
    stillrunning = open('jobs_%s_running.txt' % (jobsPrefix),'w')
    notstarted = open('jobs_%s_notstarted.txt' % (jobsPrefix),'w')
    unex = open('jobs_%s_unexplained.txt' % (jobsPrefix),'w')
    
    try:
        for jname in job_list:
            if jname in steady_ice_list:
                steady.write('%s\n' % jname)
            elif jname in plume_error_list:
                perror.write('%s\n' % jname)
            elif jname in ice_error_list:
                iceerror.write('%s\n' % jname)
            elif jname in killed_list:
                killed.write('%s\n' % jname)
            elif jname in unexplained_list:
                unex.write('%s\n' % jname)
            elif jname in still_running_list:
                stillrunning.write('%s\n' % jname)
            elif jname in not_started_list:
                notstarted.write('%s\n' % jname)
            else:
                print('not in any list:%s' % jname)

    finally:
        steady.close()
        iceerror.close()
        perror.close()
        killed.close()
        unex.close()

if __name__ == '__main__':

    if (len(sys.argv) < 2):
        print ("usage: process_jobs.py <jobprefix>")
    else:
        main(sys.argv[1])

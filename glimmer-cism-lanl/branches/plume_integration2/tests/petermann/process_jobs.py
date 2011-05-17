#!/usr/bin/env python


import os
import re
import gcplume
import pickle
import subprocess

UNEXPLAINED = 0
PLUME_ERROR = 1
ICE_ERROR = 2
STEADY_STATE = 3
TIMED_OUT = 4

jobs_dir  = os.path.join(os.path.expandvars('$GC_JOBS'))
all_jobs = os.path.os.listdir(jobs_dir)
apr_jobs = [jdir for jdir in all_jobs if jdir.startswith('apr22.2')]

steady_ice_list = []
plume_error_list = []
ice_error_list = []
killed_list = []
unexplained_list = []
job_list = []

pattern = 'apr22.2_k_([0-9\.]*)_amp_([0-9\.]*)_upvel_([0-9\.]*)_temp_([\-0-9\.]*)_itemp_([\-0-9\.]*)_tvel_([0-9\.]*)'
#pattern = 'may_([0-9\.]*)_amp_([0-9\.]*)_upvel_([0-9\.]*)_temp_([\-0-9\.]*)_itemp_([\-0-9\.]*)_tvel_([0-9\.]*)'

for jname in apr_jobs:

    unexplained = True
    
    g = re.match(pattern,jname).group

    job_key = (float(g(1)),float(g(2)),float(g(3)),
               float(g(4)),float(g(5)),float(g(6)))
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

    cmd = ['tail','-n9',log_file]
    log_tail_p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
    log_tail_p.wait()
    final_log_lines = log_tail_p.stdout.readlines()
    
    for l in final_log_lines:
        if l.strip().startswith('Stopped time-stepping'):
            unexplained = False
            steady_ice_list.append(job_key)

        if 'ERROR' in l:
            unexplained = False
            ice_error_list.append(job_key)
            

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

    if unexplained:
        unexplained_list.append((j_obj,job_key))

print('%s in all apr_jobs' % len(apr_jobs))        
print('%s in steady_ice_list' % len(steady_ice_list))
print('%s in ice_error_list' % len(ice_error_list))
print('%s in plume_error_list' % len(plume_error_list))
print('%s in killed_list' % len(killed_list))
print('%s in unexplained_list' % len(unexplained_list))

print(len(set(steady_ice_list).intersection(set(ice_error_list))))
print(len(set(steady_ice_list).intersection(set(plume_error_list))))
print(len(set(steady_ice_list).intersection(set(unexplained_list))))
print(len(set(steady_ice_list).intersection(set(killed_list))))
print(len(set(ice_error_list).intersection(set(plume_error_list))))
print(len(set(ice_error_list).intersection(set(unexplained_list))))
print(len(set(ice_error_list).intersection(set(killed_list))))
print(len(set(plume_error_list).intersection(set(unexplained_list))))
print(len(set(plume_error_list).intersection(set(killed_list))))
print(len(set(killed_list).intersection(set(unexplained_list))))

#print(unexplained_list)

fjobresult = open('results.txt','w')

try:
    for job_key in job_list:
        if job_key in steady_ice_list:
            fjobresult.write('%s,%s,%s,%s,%s,%s,%s' % (job_key[0],job_key[1],job_key[2],job_key[3],
                                                       job_key[4],job_key[5], STEADY_STATE))
        elif job_key in plume_error_list:
            fjobresult.write('%s,%s,%s,%s,%s,%s,%s' % (job_key[0],job_key[1],job_key[2],job_key[3],
                                                       job_key[4],job_key[5], PLUME_ERROR))
        elif job_key in ice_error_list:
            fjobresult.write('%s,%s,%s,%s,%s,%s,%s' % (job_key[0],job_key[1],job_key[2],job_key[3],
                                                       job_key[4],job_key[5], ICE_ERROR))
        elif job_key in killed_list:
            fjobresult.write('%s,%s,%s,%s,%s,%s,%s' % (job_key[0],job_key[1],job_key[2],job_key[3],
                                                       job_key[4],job_key[5], TIMED_OUT))
        elif job_key in unexplained_list:
            fjobresult.write('%s,%s,%s,%s,%s,%s,%s' % (job_key[0],job_key[1],job_key[2],job_key[3],
                                                       job_key[4],job_key[5], UNEXPLAINED))            
        else:
            fjobresult.write('%s,%s,%s,%s,%s,%s,%s' % (job_key[0],job_key[1],job_key[2],job_key[3],
                                                       job_key[4],job_key[5], -1))

        fjobresult.write('\n')
    
finally:
    fjobresult.close()
    
     
        


#!/usr/bin/env python

import sys
import subprocess
import pickle
import os

def nc_open_jobs(jobnames):

    plist = []

    counter = 1
    for jname in jobnames:

        ncfile = os.path.join(os.path.expandvars('$GC_JOBS'),
                              jname,'%s.out.nc' % jname)

        print('showing %s of %s: %s' % (counter,len(jobnames),ncfile))
        p = subprocess.Popen(['ncview', ncfile])
        p.wait()
        counter = counter+1


if __name__ == '__main__':

    if (len(sys.argv) < 2):
        print("Usage: nc_all.py job1 job2 ... jobn")
    else:
        if (sys.argv[1].endswith('.txt')):
            f = open(sys.argv[1],'r')
            jobs = [j.strip() for j in f.readlines()]
            f.close()
        else:
            jobs = sys.argv[1:]

        print(jobs)
        nc_open_jobs(jobs)

    

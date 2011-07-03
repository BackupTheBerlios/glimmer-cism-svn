#!/usr/bin/env python

import subprocess
import sys
import os
import re
import math

def find_time_step(nc_filename):
    cmd = ['ncks','-dtime,0,1,1','-vtime',nc_filename]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p.wait()
    dump = ''.join(p.stdout.readlines())
    times = re.findall('time\[[0-9]*\]=([0-9\.]*) ',dump)
    return (float(times[1])-float(times[0]))
    
def find_num_time_slices(nc_filename):
    cmd = ['ncdump','-h',nc_filename]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p.wait()
    dump = ''.join(p.stdout.readlines())
    
    return int(re.search('UNLIMITED ; // \(([0-9]*) currently',dump).groups()[0])

def find_last_time(nc_filename):
    
    cmd = ['nc_find_last_nonnan_time',nc_filename]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p.wait()
    return int(p.stdout.readlines()[0].strip())
    
def filter_file(nc_filename,stride,lastslice,outfile):
    
    cmd = ['ncks','-dtime,0,%s,%s' % (lastslice-1,stride),'-O' ,nc_filename,outfile]
    print(' '.join(cmd))
    subprocess.call(cmd)
    

def main(new_tstep, jnames):

    ntotal = len(jnames)
    ncurr = 1

    for jname in jnames:
        
        if (jname.endswith('/')):
            jname = jname[0:-1]

        icename = os.path.join(os.path.expandvars('$GC_JOBS'),
                               '%s' % jname,
                               '%s.out.nc' % jname)
        filtered_icename = os.path.join(os.path.expandvars('$GC_JOBS'),
                               '%s' % jname,
                               '%s.out.nc' % jname)
        plumename = os.path.join(os.path.expandvars('$GC_JOBS'),
                                 '%s' % jname,
                                 'plume.%s.out.nc' % jname)
        filtered_plumename = os.path.join(os.path.expandvars('$GC_JOBS'),
                                          '%s' % jname,
                                          'plume.%s.out.nc' % jname)

        print ("%s of %s" % (ncurr, ntotal))
        sys.stdout.flush()
        ncurr = ncurr+1
        
        if (os.path.exists(icename)):

            try:
                tstep =  find_time_step(icename)
            except IndexError:
                tstep = None

            if (tstep):
                stride = int(math.floor(new_tstep/tstep))
                tindex = find_last_time(icename)
                if (stride > 1):
                    filter_file(icename,stride,tindex,filtered_icename)
            else:
                print ('skipping %s' % icename)
        else:
            print('Did not find %s' % icename)
    
        if (os.path.exists(plumename)):
            try:
                tstep = find_time_step(plumename)
            except:
                tstep = None
                
            if (tstep):
                num_plume = find_num_time_slices(plumename)
                stride = int(math.floor(new_tstep/tstep))
                if (stride > 1):
                    filter_file(plumename,stride,num_plume,filtered_plumename)
            else:
                print ('skipping %s' % plumename)
        else:
            print('Did not find %s' % plumename)
    
if __name__ == '__main__':

    USAGE = "filter_times.py <new_tstep> <jobdir1> <jobdir2> ... <jobdirN>"

    if (len(sys.argv) < 2):
        print (USAGE)

    new_tstep = float(sys.argv[1])
    jnames = sys.argv[2:]
    print ('new time step: %s' % new_tstep)
    print ('applied to: \n%s' % '\n'.join(jnames))

    main(new_tstep, jnames)

#!/usr/bin/python


import sys
import subprocess
import re
import os

USAGE = 'nc_last_time_slice.py <nc_file_in>*'


def process_files(fnames):
    
    for fname in fnames:
        #print('scanning %s' % fname)
        p = subprocess.Popen(['ncdump','-h',fname],stdout=subprocess.PIPE)
        (output,err) = p.communicate()
        m = re.search(r'UNLIMITED ; // \(([0-9]+) currently',
                                   output)
        if (m is None):
            raise ("Could not match for number of time slices")
        
        final_time = int(m.group(1))

        cmd = ['ncra','-dtime,%s,%s' % (final_time-1,
                                                        final_time),
                               fname,
                               './last_time.%s' % os.path.basename(fname)]
        #print(cmd)
        retcode = subprocess.call(cmd)
        if (retcode != 0):
            raise Exception("%s failed" % " ".join(cmd))
        
        
        



if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print (USAGE)
        
    fnames = sys.argv[1:]

    process_files(fnames)
    




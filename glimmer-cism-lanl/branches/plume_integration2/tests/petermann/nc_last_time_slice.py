#!/usr/bin/python


import sys
import subprocess
import re
import os

USAGE = 'nc_last_time_slice.py <nc_file_in>* [-o<outputfile]'


def process_files(fnames,outputfile=None):

    if ((outputfile is not None) and (len(fnames) > 1)):
        raise Exception('Can only process a single file if given outputname')
    
    for fname in fnames:
        #print('scanning %s' % fname)
        p = subprocess.Popen(['ncdump','-h',fname],stdout=subprocess.PIPE)
        (output,err) = p.communicate()
        m = re.search(r'UNLIMITED ; // \(([0-9]+) currently',
                                   output)
        if (m is None):
            raise ("Could not match for number of time slices")
        
        final_time = int(m.group(1))
        if (outputfile is None):
            outputfile = './last_time.%s' % os.path.basename(fname)

        cmd = ['ncra','-A','-dtime,%s,%s' % (final_time-1,
                                             final_time),
                               fname,
                               outputfile]
        #print(cmd)
        retcode = subprocess.call(cmd)
        if (retcode != 0):
            raise Exception("%s failed" % " ".join(cmd))
        
        
if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print (USAGE)

    if (sys.argv[-1].startswith('-o')):
        outputfile = sys.argv[-1].split('-o')[1]
        if (len(sys.argv) < 3):
            print(USAGE)
        fnames = sys.argv[1:-1]
    else:
        fnames = sys.argv[1:]
        outputfile = None

    process_files(fnames,outputfile)
    




#!/usr/bin/env python
import subprocess
import sys

def main(joblower,jobupper,suffix):

    for jobid in range(joblower,jobupper+1):
        jobid_del = '%s%s' % (jobid,suffix)

        print(jobid_del)
        try:
            subprocess.call(['qdel', jobid_del])
        except:
            pass


if  __name__ == '__main__' :

    if (len(sys.argv) < 4):
        print("USAGE: deljob.py <job_id_lower> <job_id_upper> <suffix>")
    else:
        jlow = int(sys.argv[1])
        jup =  int(sys.argv[2])
        suffix = sys.argv[3]
        
        print('deleting %s%s to %s%s' % (jlow,suffix,jup,suffix))
        main(jlow,jup,suffix)


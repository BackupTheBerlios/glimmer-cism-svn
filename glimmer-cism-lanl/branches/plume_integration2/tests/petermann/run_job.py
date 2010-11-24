#!/usr/bin/python

import gcplume
import pickle
import sys
import time

USAGE = "run_job.py <jobfile> [-overwrite]"

if (len(sys.argv) < 2 ):
    print ('Not enough arguments.\n Usage:\n %s \n' % USAGE)
    exit(1)
	
overwrite = False

if (len(sys.argv) == 3):
   if (sys.argv[2] == '-overwrite'):
	overwrite = True
   else:
	raise Exception('Unrecognized parameter: %s' % sys.argv[2])
	
jfile = sys.argv[1]
    
f = open(jfile,'r')

try:
    try:
        j = pickle.load(f)
    finally:
        f.close()

    j.assertCanStage()
    j.resolve()
    j.stage()
    j.timeStartStr = time.ctime()
    j.timeStart = time.time()
    j.started = True
    j.error = False
    j.errorMessage = ''
    j.completed = False
    j.serialize()

    if (overwrite):
        j.run(overwrite=overwrite)
    else:
        j.run()

    j.timeStop = time.time()
    j.timeStopStr = time.ctime()

    j.completed = True
    j.error = False
    j.errorMessage = ""
    
except Exception, e:
    j.completed = False
    j.error = True
    j.errorMessage = str(e)
    print(j.errorMessage)

j.serialize()


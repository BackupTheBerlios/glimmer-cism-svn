#!/usr/bin/python

import gcplume
import pickle
import sys

USAGE = "run_job.py <jobfile>"

if (len(sys.argv) < 2 ):
 print ('Not enough arguments.\n Usage:\n %s \n' % USAGE)
 exit(1)
	
jfile = sys.argv[1]

f = open(jfile,'r')

try:
	j = pickle.load(f)
finally:
	f.close()

j.assertCanStage()
j.stage()
j.run()
j.completed = True

#j.serialize()

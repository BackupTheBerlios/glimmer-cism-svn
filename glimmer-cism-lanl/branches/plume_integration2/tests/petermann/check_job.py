#!/usr/bin/python

import gcplume
import pickle
import sys
import time

USAGE = "check_job.py [--printCompleted] <jobfile>"


if (len(sys.argv) < 2 ):

    raise Exception('Not enough arguments.\n Usage:\n %s \n' % USAGE)

if (sys.argv[1] == '--printCompleted'):
    printCompleted = True
    jfiles = sys.argv[2:]
else:
    printCompleted = False
    jfiles = sys.argv[1:]

for jfile in sys.argv[2:]:

    f = open(jfile,'r')

    try:
        j = pickle.load(f)
    finally:
        f.close()

    msg = ''

    printMsg = True

    try:
        
        if (not(j.started) and not(j.error)):
            msg = "not started"
        elif (j.started and not(j.completed or j.error)):
            msg = "running"
        elif (j.completed and not(j.error)):
            msg = "completed"
            if (not(printCompleted)):
                printMsg = False
        elif (j.error):
            msg = "error: %s" % j.errorMessage

    except AttributeError,e:
	print(e)
        msg = "old job version"
        
    if (printMsg):
        print ("%s: %s" % (j.name, msg))
            

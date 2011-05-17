#!/usr/bin/python

import sys

from gcplume import *


if (len(sys.argv) < 2):
    raise Exception("need to specify job name")

if (sys.argv[0] == 'python'):
    if (len(sys.argv) < 3):
        raise Exception("need to specify job name")
    else:
        jobNames = sys.argv[2:]
else:        
    jobNames = sys.argv[1:]

print (jobNames)

for jobName in jobNames:
    print(jobName)
    try:
        j = FromFilesJob(jobName)
        j.assertCanStage()
        j.resolve()
        j.serialize()
    except:
        print('Could not do job %s' % jobName)


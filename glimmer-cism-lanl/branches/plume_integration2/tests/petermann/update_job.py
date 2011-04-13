#!/usr/bin/python

import sys

from gcplume import *


if (len(sys.argv) < 2):
    raise Exception("need to specify job name")

if (sys.argv[0] == 'python'):
    if (len(sys.argv) < 3):
        raise Exception("need to specify job name")
    else:
        jobName = sys.argv[2]
else:        
    jobName = sys.argv[1]

j = FromFilesJob(jobName)

j.assertCanStage()
j.resolve()
j.serialize()


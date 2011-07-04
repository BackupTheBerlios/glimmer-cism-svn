#!/usr/bin/env python

import sys
import subprocess



files = sys.argv[1:]

for f in files:
    print(f)
    cmd = ['sh','epstopdf',f]
    r = subprocess.call(cmd)
    
    
    

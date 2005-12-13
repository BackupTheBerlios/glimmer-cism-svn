#!/usr/bin/env python

# Copyright 2004, Magnus Hagdorn

"""Create initial sediment distribution"""

import sys,PyCF, Numeric

# creating option parser
parser = PyCF.CFOptParser()
parser.time()
opts = PyCF.CFOptions(parser,2)

# open input file
infile = opts.cffile()
time = opts.times(infile)
topg = infile.getvar('topg')
# creat output file
outfile = infile.clone(opts.args[1])

# create sediment variable
seds = outfile.createVariable('seds3')

# load time slice
topo = topg.get2Dfield(time)
outfile.file.variables['time'][0] = infile.file.variables['time'][time]

# set sediment thickness
data = Numeric.where(topo<1000.,5,0)
data = Numeric.where(topo>-400.,data,0)
seds[0,:,:] = Numeric.transpose(data).astype(Numeric.Float32)  

outfile.close()

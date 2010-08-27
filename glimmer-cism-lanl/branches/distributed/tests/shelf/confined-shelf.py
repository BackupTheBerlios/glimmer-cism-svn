#!/usr/bin/env python
# This script runs a "Confined Shelf Experiment".
# Files are written in the "output" subdirectory.
# The script performs the following three steps:
# 1. Create a netCDF input file for Glimmer.
# 2. Run Glimmer, creating a netCDF output file.
# 3. Move any additional files written by Glimmer to the "scratch" subdirectory.
# Written by Glen Granzow at the University of Montana on April 8, 2010

import sys, os, glob, shutil, numpy
from netCDF import *
from ConfigParser import ConfigParser

# Check to see if a config file was specified on the command line.
# If not, confined-shelf.config is used.
if len(sys.argv) > 1:
  if sys.argv[1][0] == '-': # The filename can't begin with a hyphen
    print '\nUsage:  python confined-shelf.py [FILE.CONFIG]\n'
    sys.exit(0)
  else:
    configfile = sys.argv[1]
else:
  configfile = 'confined-shelf.config'

# Create a netCDF file according to the information in the config file.
parser = ConfigParser()
parser.read(configfile)
nx = int(parser.get('grid','ewn'))
ny = int(parser.get('grid','nsn'))
nz = int(parser.get('grid','upn'))
dx = float(parser.get('grid','dew'))
dy = float(parser.get('grid','dns'))
filename = parser.get('CF input', 'name')

print 'Writing', filename
try:
  netCDFfile = NetCDFFile(filename,'w',format='NETCDF3_CLASSIC')
except TypeError:
  netCDFfile = NetCDFFile(filename,'w')

netCDFfile.createDimension('time',1)
netCDFfile.createDimension('x1',nx)
netCDFfile.createDimension('y1',ny)
netCDFfile.createDimension('level',nz)
netCDFfile.createDimension('x0',nx-1) # staggered grid 
netCDFfile.createDimension('y0',ny-1)

x = dx*numpy.arange(nx,dtype='float32')
y = dx*numpy.arange(ny,dtype='float32')

netCDFfile.createVariable('time','f',('time',))[:] = [0]
netCDFfile.createVariable('x1','f',('x1',))[:] = x.tolist()
netCDFfile.createVariable('y1','f',('y1',))[:] = y.tolist()

netCDFfile.createVariable('x0','f',('x0',))[:] = (dx/2 + x[:-1]).tolist()
netCDFfile.createVariable('y0','f',('y0',))[:] = (dy/2 + y[:-1]).tolist()

# Check to make sure that the flow law parameter in the config file is correct.
default_flwa = float(parser.get('parameters','default_flwa'))
if default_flwa != 4.6e-18:
  print 'WARNING: The parameter default_flwa in',configfile,'should be 4.6e-18'
  print '         Currently it is',default_flwa

# Determine from the config file whether periodic boundary conditions are to be
# imposed in the x direction.  
periodic_ew = int(parser.get('options','periodic_ew'))

# Calculate values for the required variables.
thk  = numpy.zeros([ny,nx],dtype='float32')
beta = numpy.empty([ny-1,nx-1],dtype='float32')
kbc  = numpy.zeros([ny,nx],dtype='int')
zero = numpy.zeros([nz,ny-1,nx-1],dtype='float32')
thk[4:-2,2:-2] = 1000
beta[:,:] = 9.96921e+36
beta[2:nx-2,2:ny-2] = 0
kbc[ny-3:,:]  = 1
if not periodic_ew:
  kbc[:,:3] = 1
  kbc[:,nx-3:] = 1

# Create the required variables in the netCDF file.
netCDFfile.createVariable('thk',      'f',('time','y1','x1'))[:] = thk.tolist()
netCDFfile.createVariable('kinbcmask','i',('time','y1','x1'))[:] = kbc.tolist()
netCDFfile.createVariable('topg',     'f',('time','y1','x1'))[:] = ny*[nx*[-2000]]
netCDFfile.createVariable('beta',     'f',('time','y0','x0'))[:] = beta.tolist()
netCDFfile.createVariable('uvelhom',  'f',('time','level','y0','x0'))[:] = zero.tolist()
netCDFfile.createVariable('vvelhom',  'f',('time','level','y0','x0'))[:] = zero.tolist()

netCDFfile.close()

# Run Glimmer
print 'Running Glimmer/CISM'
os.system('echo '+configfile+' | simple_glide')

# Clean up by moving extra files written by Glimmer to the "scratch" subdirectory
# Look for files with extension "txt", "log", or "nc"
for files in glob.glob('*.txt')+glob.glob('*.log')+glob.glob('*.nc'):
# Delete any files already in scratch with these filenames 
  if files in os.listdir('scratch'):
    os.remove(os.path.join('scratch',files))
# Move the new files to scratch
  shutil.move(files,'scratch')



#!/usr/bin/env python

# File to generate an input netcdf file readable by Glimmer-CISM,
# using a .mat file generated using an appropriate matlab script 
#
# Sytnax: 
# ./brpoc.py bproc.config
#
# Script will parse the grid section of the config file to produce proper output.

from numpy import array,zeros,size
import sys
from scipy.io import loadmat
try:
  from netCDF4 import Dataset
except ImportError:
  from Scientific.IO.NetCDF import NetCDFFile as Dataset

def nc_from_config(configFilename):
  from ConfigParser import ConfigParser
  parser = ConfigParser()
  parser.read(configFilename)
  nx = int(parser.get("grid","ewn"))
  ny = int(parser.get("grid","nsn"))
  nz = int(parser.get("grid","upn"))
  dx = float(parser.get("grid","dew"))
  dy = float(parser.get("grid","dns"))
  filename = parser.get("CF input", "name")

  print "Writing to", filename
  try:
    netCDFfile = Dataset(filename,'w',format='NETCDF3_CLASSIC')
  except TypeError:
    netCDFfile = Dataset(filename,'w')
  netCDFfile.createDimension('time',1)
  netCDFfile.createDimension('x1',nx)
  netCDFfile.createDimension('y1',ny)
  netCDFfile.createDimension('level',nz)
  netCDFfile.createDimension('x0',nx-1) # staggered grid 
  netCDFfile.createDimension('y0',ny-1)
  from numpy import arange
  x = dx*arange(nx,dtype='float32')
  y = dx*arange(ny,dtype='float32')
  netCDFfile.createVariable('time','f',('time',))[:] = [0]
  netCDFfile.createVariable('x1','f',('x1',))[:] = x
  netCDFfile.createVariable('y1','f',('y1',))[:] = y
  netCDFfile.createVariable('x0','f',('x0',))[:] = dx/2 + x[:-1] # staggered grid
  netCDFfile.createVariable('y0','f',('y0',))[:] = dy/2 + y[:-1]
  return netCDFfile, (nx,ny,nz,dx,dy)

#Parse the config file to determine how to set up the netcdf file
nc, shape = nc_from_config(sys.argv[1])

topg     = nc.createVariable('topg',    'f',('time','y1','x1'))
thk      = nc.createVariable('thk',     'f',('time','y1','x1'))
usurf    = nc.createVariable('usurf',   'f',('time','y1','x1'))
artm     = nc.createVariable('artm',    'f',('time','y1','x1'))
bheatflx = nc.createVariable('bheatflx','f',('time','y1','x1'))
beta     = nc.createVariable('beta',    'f',('time','y0','x0'))
tauf     = nc.createVariable('tauf',    'f',('time','y0','x0'))

kinbcmask = nc.createVariable('kinbcmask','f',('time','y0','x0'))
uvelhom   = nc.createVariable('uvelhom',  'f',('time','level','y0','x0'))
vvelhom   = nc.createVariable('vvelhom',  'f',('time','level','y0','x0'))

#Read the data:
d = loadmat('bproc.mat')

#Set the fields
topg[0,:,:] = array(d['topg'],dtype='float32')
thk[0,:,:] = array(d['thck'],dtype='float32')
usurf[0,:,:] = array(d['usrf'],dtype='float32')
artm[0,:,:] = array(d['airt'],dtype='float32')
bheatflx[0,:,:] = array(d['qgeo'],dtype='float32')
beta[0,:,:] = array(d['beta'],dtype='float32')

tauf[0,:,:]   = array(d['tauf'],dtype='float32')
kinbcmask[0,:,:] = array(d['kinbcmask'],dtype='int8')
uvelhom[0,:,:,:] = array(d['uvelhom'],dtype='float32')
vvelhom[0,:,:,:] = array(d['vvelhom'],dtype='float32')

nc.close()

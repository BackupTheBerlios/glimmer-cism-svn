#!/usr/bin/env python

# File to generate an input netcdf file readable by Glimmer-CISM,
# using a .mat file generated using an appropriate matlab script 
#
# Sytnax: 
# ./brpoc.py bproc.config
#
# Script will parse the grid section of the config file to produce proper output.

import glimcdf
from math import sin,cos,tan,pi
from numpy import array,zeros,size
import sys
from scipy.io import loadmat

#Parse the config file to determine how to set up the netcdf file
nc, shape = glimcdf.nc_from_config(sys.argv[1])

#Create variables that will appear:
topg = glimcdf.setup_variable(nc, "topg")
thk = glimcdf.setup_variable(nc, "thk")
usurf = glimcdf.setup_variable(nc, "usrf")
artm = glimcdf.setup_variable(nc, "artm")
bheatflx = glimcdf.setup_variable(nc, "bheatflx")
beta = glimcdf.setup_variable(nc, "beta", staggered=True)

#kinbcmask = glimcdf.setup_variable(nc, "kinbcmask", staggered=True)
#uvelhom = glimcdf.setup_variable(nc, "uvelhom", staggered=True, useZ=True)
#vvelhom = glimcdf.setup_variable(nc, "vvelhom", staggered=True, useZ=True)

#Read the data:
d = loadmat('bproc.mat')

#Set the fields
topg[0,:,:] = array(d['topg'],dtype='float32')

#temp = d['usrf_10km'] - d['topg_10km']
#temp[temp<0.]=0.
#thk[0,:,:] = array(temp,dtype='float32')

thk[0,:,:] = array(d['thck'],dtype='float32')
usurf[0,:,:] = array(d['usrf'],dtype='float32')
artm[0,:,:] = array(d['airt'],dtype='float32')
bheatflx[0,:,:] = array(d['qgeo'],dtype='float32')
beta[0,:,:] = array(d['beta'],dtype='float32')

#kinbcmask[0,:,:] = array(d['kinbcmask'],dtype='int8')
#uvelhom[0,:,:,:] = array(d['uvelhom'],dtype='float32')
#vvelhom[0,:,:,:] = array(d['vvelhom'],dtype='float32')

nc.close()

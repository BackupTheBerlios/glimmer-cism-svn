#!/usr/bin/env python
# This script runs ISMIP-HOM experiments using Glimmer-CISM.
# Output files are written in the "output" subdirectory.
# The script loops over experiments performing the following three steps:
# 1. Create two input files, a configuration file and a netCDF input file.
# 2. Run Glimmer, creating a netCDF output file.
# 3. Write a standard ISMIP-HOM file from the data in the netCDF output file.
# When finished, any additional files written by Glimmer are moved to the "scratch" subdirectory.

# After running this script, run plotISMIPHOM.py to plot the results.
# See the accompanying README file for more information.
# To see all command line options run: python runISMIPHOM.py --help
# Written March 2, 2010 by Glen Granzow at the University of Montana.
# Update for the distributed branch (parallel option) on Jaguar 
# April 6, 2011 by Kate Evans at Oak Ridge Nat'l Lab

# OptionParser callback that splits a comma-separated list.
def appendToList(option,str_opt,value,parser):
  listOfValues = getattr(parser.values,option.dest)
  if listOfValues == None: listOfValues = list()
  listOfValues += value.split(',')
  setattr(parser.values,option.dest,listOfValues)

defaultExperiments = ['a']      # ['a','b','c','d']
defaultSizes = ['10','20','40'] # ['5','10','20','40','80','160']
defaultProcs = ['0'] # ['12']

if __name__ == '__main__':
  import os
  import glob
  import shutil
  import ConfigParser
  from optparse import OptionParser
  from math import tan, sin, pi, exp
  from netCDF import *

# Parse the command line arguments
  parser = OptionParser()
  parser.add_option('-e','--exp',dest='experiments',type='string',action='callback',callback=appendToList,help='Specify ISMIP-HOM experiments to run')
  parser.add_option('-s','--size',dest='sizes',type='string',action='callback',callback=appendToList,help='Specify domain sizes to run')
  parser.add_option('-g','--grid-size',dest='horizontal_grid_size',type='int',help='(overrides ewn and nsn in ishom.a.config)')
  parser.add_option('-v','--vert-grid-size',dest='vertical_grid_size',type='int',help='(overrides upn in ishom.a.config)')
  parser.add_option('-d','--diagnostic-scheme',dest='diagnostic_scheme',help='(overrides ishom.a.config)')
  parser.add_option('-r','--run',dest='executable',default='simple_glide',help='Set path to the GLIMMER executable (defaults to simple_glide)')
  parser.add_option('-n','--procs',dest='procs',help='Set the number of processors to run simple_glide. (default is 0, serial run, setting j=1 means MPI on but only 1 processor)')
  parser.add_option('-p','--prefix',dest='prefix',default='glm1',help='Prefix to use for model output files (defaults to glm1)')
  parser.add_option('-f','--format-only',dest='format_only',action='store_true',help='Generate the config and NetCDF input files only')
  options, args = parser.parse_args()
# If the user didn't specify a list of experiments or domain sizes, run the whole suite
  if options.experiments == None: options.experiments = defaultExperiments
  if options.sizes == None: options.sizes = defaultSizes
  if options.procs == None: options.procs = defaultProcs

# Loop over the experiments requested on the command line
  for experiment in options.experiments:

#   Loop over the sizes requested on the command line
    for size in map(int,options.sizes):

      if experiment == 'f': 
        if size != 100 or len(options.sizes) > 1:
          print 'NOTE: Experiment f uses a domain size of 100 km only'
        size = 100

#     Create a configuration file by copying an existing example
      configParser = ConfigParser.SafeConfigParser()
      configParser.read('ishom.'+experiment+'.config')
#     Set or get the number of horizontal grid points in each direction
      if options.horizontal_grid_size == None:
        nx = int(configParser.get('grid','ewn'))
        ny = int(configParser.get('grid','nsn'))
      else:
        nx = ny = options.horizontal_grid_size
        configParser.set('grid','ewn',str(nx))
        configParser.set('grid','nsn',str(ny))
#     Set the grid spacing in the horizontal directions
#     Glimmer's periodic boundary conditions seem to be messed up so
#     the grids used here are something of a hack.
#     They were copied from a previous script (verify.py)
      dx = size*1000/(nx-1)
      dy = size*1000/(ny-1)
      configParser.set('grid','dew',str(dx))
      configParser.set('grid','dns',str(dy))
#     Specify the netCDF input and output filenames
#     All files will be written in the "output" subdirectory
      filename = os.path.join('output','ishom.'+experiment+'.'+str(size)+'km')
      configParser.set('CF input', 'name',filename+'.nc')
      configParser.set('CF output','name',filename+'.out.nc')
#     Make additional changes if requested on the command line
      if options.vertical_grid_size != None:
        configParser.set('grid','upn',str(options.vertical_grid_size))
      if options.diagnostic_scheme != None:
        configParser.set('ho_options','diagnostic_scheme',options.diagnostic_scheme)
#     Write the new configuration file
      configFile = open(filename+'.config','w')
      configParser.write(configFile)
      configFile.close()

#     Create the netCDF input file needed by Glimmer
      if netCDF_module == 'netCDF4':
        netCDFfile = NetCDFFile(filename+'.nc','w',format='NETCDF3_CLASSIC')
      else:
        netCDFfile = NetCDFFile(filename+'.nc','w')
      netCDFfile.createDimension('time',1)
      netCDFfile.createDimension('x1',nx)   # unstaggered grid
      netCDFfile.createDimension('y1',ny)
      netCDFfile.createDimension('x0',nx-1) # staggered grid 
      netCDFfile.createDimension('y0',ny-1)
      time = netCDFfile.createVariable('time','f',('time',))
      x1   = netCDFfile.createVariable('x1','f',('x1',))
      y1   = netCDFfile.createVariable('y1','f',('y1',))
      x0   = netCDFfile.createVariable('x0','f',('x0',))
      y0   = netCDFfile.createVariable('y0','f',('y0',))
      thk  = netCDFfile.createVariable('thk' ,'f',('time','y1','x1'))
      topg = netCDFfile.createVariable('topg','f',('time','y1','x1'))
      if experiment in ('c','d'):
        beta = netCDFfile.createVariable('beta','f',('time','y0','x0'))
      time[0] = 0
      x1[:] = [i*dx for i in range(nx)] # unstaggered grid
      y1[:] = [j*dy for j in range(ny)]
      x0[:] = [(i+0.5)*dx for i in range(nx-1)] # staggered grid 
      y0[:] = [(j+0.5)*dy for j in range(ny-1)]

#     Generate the ice thickness, bed topography, and (sometimes) 
#     basal friction coefficient for the experiment


      thickness  = list()
      topography = list()
      basalFriction = list()

      if experiment in ('a','b'):
        xx = [i*dx for i in range(nx)]
        yy = [j*dy for j in range(ny)]
        xx = [x *float(nx-1)/float(nx-2) for x in xx]
        yy = [y *float(ny-1)/float(ny-2) for y in yy]
        alpha = 0.5 * pi/180
        zz = [4000-i*dx*tan(alpha) for i in range(nx)]
      elif experiment in ('c','d'):
        xx = [i*dx for i in range(nx)]
        yy = [j*dy for j in range(ny)]
        xx = [x *float(nx-1)/float(nx-3) for x in xx]
        yy = [y *float(ny-1)/float(ny-3) for y in yy]
        alpha = 0.1 * pi/180
        zz = [1000-i*dx*tan(alpha) for i in range(nx)]
      elif experiment == 'f':
        xx = [i*dx for i in range(nx)]
        yy = [j*dy for j in range(ny)]
        alpha = 3.0 * pi/180
        zz = [6000-i*dx*tan(alpha) for i in range(nx)]
        xc = (xx[0]+xx[-1])/2
        yc = (yy[0]+yy[-1])/2
        a0 = 100
        sigma2 = 10000**2

      omega = 2*pi / (size*1000)
      for y in yy:
        row = list()
        for x in xx:
          if experiment == 'a':
            row.append(1000 - 500*sin(omega*x)*sin(omega*y))
          elif experiment == 'b':
            row.append(1000 - 500*sin(omega*x))
          elif experiment == 'c':
            row.append(1000 + 1000*sin(omega*x)*sin(omega*y))
          elif experiment == 'd':
            row.append(1000 + 1000*sin(omega*x))
          elif experiment == 'f':
            row.append(1000 - a0*exp(-((x-xc)**2+(y-yc)**2)/sigma2))
        if experiment in ('a','b','f'):
          thickness.append(row)
          topography.append([z-t for (z,t) in zip(zz,row)])
        else:
          basalFriction.append(row[:-1])
 
      if experiment in ('a','b','f'):
        thk [:] = thickness
        topg[:] = topography
      elif experiment in ('c','d'):
        thk [:] = ny*[nx*[1000]]
        topg[:] = ny*[zz]
        beta[:] = basalFriction[:-1]
      netCDFfile.close()

      if not options.format_only:

#       Run Glimmer
        print 'Running',options.executable,'for experiment',experiment.upper(),'with domain size',size,'km'
#       Serially
      if options.procs == '0': 
        exitCode = os.system('echo '+filename+'.config'+' | '+options.executable)
      else:
# In parallel, right now can only handle aprun command
        exitCode = os.system('time aprun -n'+options.procs+' '+options.executable+' '+filename+'.config')

        if exitCode == 0:
#         Extract the output data for comparison to the other models
          if experiment == 'a': # Get the following (variable,level)'s
            variables = [('uvelhom',0),('vvelhom',0),('wvel',0),('tau_xz',-1),('tau_yz',-1)]
          if experiment == 'b':
            variables = [('uvelhom',0),('wvel',0),('tau_xz',-1)]
          if experiment == 'c':
            variables = [('uvelhom',0),('vvelhom',0),('wvel',0),('uvelhom',-1),('vvelhom',-1),('tau_xz',-1),('tau_yz',-1)]
          if experiment == 'd':
            variables = [('uvelhom',0),('wvel',0),('uvelhom',-1),('tau_xz',-1)]
          if experiment == 'f':
            variables = [('usurf',None),('uvelhom',0),('vvelhom',0),('wvel',0)]
#         Open the netCDF file that was written by Glimmer
          netCDFfile = NetCDFFile(filename+'.out.nc','r')
          data = [(netCDFfile.variables[v[0]],v[1]) for v in variables]

#         Write a "standard" ISMIP-HOM file (example file name: "glm1a020.txt") in the "output" subdirectory 
          ISMIP_HOMfilename = os.path.join('output',options.prefix+experiment+'%03d'%size+'.txt')
          ISMIP_HOMfile = open(ISMIP_HOMfilename,'w')
          for i in range(nx-2):
            x = float(i)/(nx-3)   # In a more perfect world: x = (i+0.5)/(nx-2)
            for j in range(ny-2):
              y = float(j)/(ny-3) # In a more perfect world: y = (j+0.5)/(ny-2)
              if netCDF_module == 'Scientific.IO.NetCDF':
                if experiment in ('a','c'):
                  ISMIP_HOMfile.write('\t'.join(map(str,[x,y]+[v[0,level,j,i][0] for (v,level) in data]))+'\n')
                if experiment in ('b','d','e'):
                  ISMIP_HOMfile.write('\t'.join(map(str,[x]+[v[0,level,ny/2,i][0] for (v,level) in data]))+'\n')
                elif experiment == 'f':
                  ISMIP_HOMfile.write('\t'.join(map(str,[x,y,data[0][0][-1,j,i][0]]+[v[-1,level,j,i][0] for (v,level) in data[1:]]))+'\n')
              else:
                if experiment in ('a','c'):
                  ISMIP_HOMfile.write('\t'.join(map(str,[x,y]+[v[0,level,j,i] for (v,level) in data]))+'\n')
                if experiment in ('b','d','e'):
                  ISMIP_HOMfile.write('\t'.join(map(str,[x]+[v[0,level,ny/2,i] for (v,level) in data]))+'\n')
                elif experiment == 'f':
                  ISMIP_HOMfile.write('\t'.join(map(str,[x,y,data[0][0][-1,j,i]]+[v[-1,level,j,i] for (v,level) in data[1:]]))+'\n')
          ISMIP_HOMfile.close()
          netCDFfile.close()

        else:

          print 'OH NO! Glimmer seems to have failed. (I hate it when that happens.)'

#     Experiment f should be run for one size (100 km) only
      if experiment == 'f': break

# Clean up by moving extra files written by Glimmer to the "scratch" subdirectory
# Look for files with extension "txt", "log", or "nc"
  for files in glob.glob('*.txt')+glob.glob('*.log')+glob.glob('*.nc'):
#   Delete any files already in scratch with these filenames 
    if files in os.listdir('scratch'):
      os.remove(os.path.join('scratch',files))
#   Move the new files to scratch
    shutil.move(files,'scratch')


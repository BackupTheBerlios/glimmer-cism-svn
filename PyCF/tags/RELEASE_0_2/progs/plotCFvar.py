#! /usr/bin/env python

# Copyright 2004, Magnus Hagdorn
# 
# This file is part of GLIMMER.
# 
# GLIMMER is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# GLIMMER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with GLIMMER; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

"""A simple plot of CF fields."""

import PyGMT,PyCF,sys

# creating option parser
parser = PyCF.CFOptParser()
parser.variable()
parser.profile_file()
parser.time()
parser.region()
parser.plot()
opts = PyCF.CFOptions(parser,-2)

if opts.options.level == None:
    level = 0
else:
    level = opts.options.level

count = 0
numplots = 0
if opts.nvars>1:
    count = count + 1
    numplots = opts.nvars
if opts.ntimes>1:
    count = count + 1
    numplots = opts.ntimes
if opts.nfiles>1:
    count = count + 1
    numplots = opts.nfiles
    
if count > 1:
    sys.stderr.write('Error, can only have either more than one time slice or more than one variable or more than one file!\n')
    sys.exit(1)

if opts.options.profname!=None:
    infile = opts.cfprofile()
else:
    infile = opts.cffile()
var  = opts.vars(infile)
time = opts.times(infile)
# get number of plots
deltax = 1.
deltay = 1.
sizex = opts.options.width+deltax
sizey = infile.aspect_ratio*opts.options.width+deltay

plot=None
if numplots > 1:
    numx = int((opts.papersize[0])/(sizex))
    numy = int((opts.papersize[1])/(sizey))
    numpages = int(float(numplots-0.1)/float(numx*numy))
    p=-1
    for i in range(0,numplots):
        if i%(numx*numy)==0:
            # need to open a new plot file
            if plot!=None:
                plot.close()
            if numpages>0:
                p=p+1
                plot = opts.plot(number=p)
            else:
                p=0
                plot = opts.plot()
            bigarea = PyGMT.AreaXY(plot,size=opts.papersize)

        if opts.nfiles>1:
            if opts.options.profname!=None:
                infile = opts.cfprofile(i)
            else:
                infile = opts.cffile(i)
        if opts.nvars>1:
            var  = opts.vars(infile,i)
        else:
            var  = opts.vars(infile)
        if opts.ntimes>1:
            time = opts.times(infile,i)
        else:
            time = opts.times(infile)
        x = i%numx
        y = int((i-p*(numx*numy))/numx)
        area = PyCF.CFArea(bigarea,infile,pos=[x*sizex,opts.papersize[1]-(y+1)*sizey],size=sizex-deltax)
        if opts.options.land:
            area.land(time)
        area.image(var,time,clip = opts.options.clip,level=level,mono=opts.options.mono)
        area.coastline()
        try:
            thk = infile.getvar('thk')
            area.contour(thk,[0.1],'-W2/0/0/0',time)
        except:
            pass
        if parser.profile!=None:
            area.profile(args='-W5/0/0/0')
        area.axis='wesn'
        area.coordsystem()
        area.printinfo(time)
else:
    plot = opts.plot()
    area = PyCF.CFArea(plot,infile,pos=[0.,3.],size=sizex-deltax)
    if opts.options.land:
        area.land(time)
    area.image(var,time,clip = opts.options.clip,level=level,mono=opts.options.mono)
    if var.name=="vel":
        area.velocity_field(time,level=level)
    area.coastline()
    try:
        thk = infile.getvar('thk')
        area.contour(thk,[0.1],'-W2/0/0/0',time)
    except:
        pass
    if var.name == 'is' or var.name == 'thk':
        area.contour(var,[500,1000,2500,3000],'-W1/255/255/255',time)
    if parser.profile!=None:
        area.profile(args='-W5/0/0/0')
    area.coordsystem()
    area.printinfo(time)
    if opts.options.dolegend:
        PyGMT.colourkey(area,var.colourmap.cptfile,title=var.long_name,pos=[0,-2])
    
plot.close()

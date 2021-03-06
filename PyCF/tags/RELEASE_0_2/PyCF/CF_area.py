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

"""Class for plotting ISM grid files."""

__all__ = ['CFArea']

import PyGMT,Numeric,math
from CF_loadfile import CFvariable
from CF_colourmap import CFcolours
from StringIO import StringIO

class CFArea(PyGMT.AreaXY):
    """CF grid plotting area."""

    def __init__(self,parent,cffile,pos=[0.,0.],size=10):
        """Initialising ISM area.

        parent: can be either a Canvas or another Area.
        cffile: CF file used for setting up projection and area
        pos: position of area relative to the parent
        size: size of GMT area
        """

        # initialising geographic area
        self.file = cffile
        if cffile.projection == 'lin':
            self.geo = PyGMT.AreaXY(parent,pos=pos,size=[size,cffile.aspect_ratio*size])
            self.geo.setregion([cffile.ll_geo[0]/cffile.deltax,cffile.ll_geo[1]/cffile.deltay],
                               [cffile.ur_geo[0]/cffile.deltax,cffile.ur_geo[1]/cffile.deltay])
        else:
            self.geo = PyGMT.AreaGEO(parent,cffile.projection.getGMTprojection(mapwidth=size).upper(), pos=pos, size=size)
            self.geo.setregion(cffile.ll_geo, cffile.ur_geo)
        # initialising paper area
        self.paper = PyGMT.AreaXY(parent,pos=pos,size=self.geo.size)

        # initialising XY area
        PyGMT.AreaXY.__init__(self,parent,pos=pos,size=self.geo.size)
        self.setregion(cffile.ll_xy, cffile.ur_xy)

    def coordsystem(self,grid=True):
        """Plot coordinate system.

        grid: if true plot lat/long grid"""
        self.geo.axis = self.axis
        self.geo.coordsystem(grid=grid)

    def coastline(self,args='-W -A0/1/1'):
        """Plot coastline.

        args: arguments passed on to pscoast."""
        try:
            self.geo.coastline(args)
        except:
            pass

    def stamp(self,text):
        """Print text in lower left corner."""

        self.paper.text([0.15,0.15],text,textargs='8 0 0 LB',comargs='-W255/255/255o')

    def printinfo(self,time):
        """Print a data name and time slice."""

        self.stamp('%s   %.2fka'%(self.file.title,self.file.time(time)))
        
    def image(self,var,time,level=0,clip=None,mono=False):
        """Plot a colour map.

        var: CFvariable
        time: time slice
        level: horizontal slice
        clip: only display data where clip>0.
        mono: convert colour to mono
        """
        
        clipped = False
        if clip in ['topg','thk','usurf'] :
            cvar = CFvariable(var.cffile,clip)
            self.clip(cvar.getGMTgrid(time),0.1)
            clipped = True
        if mono:
            args="-M"
        else:
            args=""
        PyGMT.AreaXY.image(self,var.getGMTgrid(time,level=level),var.colourmap.cptfile,args=args)
        if clipped:
            self.unclip()

    def velocity_field(self,time,level=0,mins=10.):
        """Plot vectors of velocity field

        time: time slice
        level: horizontal slice
        mins: do not plot vectors below this size (default 0.01)
        """

        # get velocity components
        velx = self.file.getvar('uvel')
        vely  = self.file.getvar('vvel')
        vel = self.file.getvar('vel')

        data = vel.get2Dfield(time,level=level)
        datax = velx.get2Dfield(time,level=level)
        datay = vely.get2Dfield(time,level=level)

        # calculate node spacing
        vector_density = 0.2
        x_spacing = int(((self.ur[0]-self.ll[0])/(self.size[0]*self.file.deltax))*vector_density)+1
        y_spacing = int(((self.ur[1]-self.ll[1])/(self.size[1]*self.file.deltay))*vector_density)+1

        fact = 360./(2*math.pi)
        outstring = StringIO()
        for i in range(int(self.ll[0]/self.file.deltax),int(self.ur[0]/self.file.deltax),x_spacing):
            for j in range(int(self.ll[1]/self.file.deltay),int(self.ur[1]/self.file.deltay),y_spacing):
                if (data[i,j]>mins):
                    outstring.write('%f %f %f %f\n'%(velx.xdim[i], velx.ydim[j], fact*math.atan2(datay[i,j],datax[i,j]), 0.3))

        self.canvascom('psxy',' -G0 -Sv0.01/0.1/0.1',indata=outstring.getvalue())
        
    def contour(self,var,contours,args,time,level=0):
        """Plot a contour map.

        var: CFvariable
        contours: list of contour intervals
        args: further arguments
        time: time slice
        level: horizontal slice."""

        PyGMT.AreaXY.contour(self,var.getGMTgrid(time,level=level),contours,args)

    def land(self,time,grey='240'):
        """Plot area above sea level.

        time: time slice
        grey: grey value."""

        cvar = CFvariable(self.file,'topg')
        self.clip(cvar.getGMTgrid(time),0.1)
        self.line('-W -L -G%s'%grey,[self.ll[0],self.ur[0],self.ur[0],self.ll[0]],[self.ll[1],self.ll[1],self.ur[1],self.ur[1]])
        self.unclip()

    def profile(self,args='-W1/0/0/0'):
        """Plot profile if present in file.

        args: pen arguments default -W1/0/0/0"""

        if hasattr(self.file,'interpolated'):
            self.line(args,self.file.interpolated[0,:].tolist(),self.file.interpolated[1,:].tolist())

    def rsl_locations(self,rsldb,dataset=None):
        """Plot RSL locations.

        rsldb: RSL database data set
        dataset: if not None, list of dataset ids to be plotted (when None, plot all)."""

        rsldata = rsldb.getLocationRange(self.file.minmax_long,self.file.minmax_lat)
        for loc in rsldata:
            self.geo.plotsymbol([loc[3]],[loc[4]],size=0.2,symbol='a',args='-G%s'%CFcolours[loc[1]])

if __name__ == '__main__':
    from CF_options import *
    from CF_IOrsl import *
    
    parser = CFOptParser()
    parser.variable()
    parser.profile(vars=False)
    parser.time()
    parser.region()
    parser.plot()
    opts = CFOptions(parser,2)
    if parser.profile!=None:
        infile = opts.cfprofile()
    else:
        infile = opts.cffile()
    var = opts.vars(infile)
    ts = opts.times(infile)
    plot = opts.plot()
    area = CFArea(plot,infile)
    if opts.options.land:
        area.land(ts)
    area.image(var,ts,clip = opts.options.clip)
    area.coastline()
    if parser.profile!=None:
        area.profile(args='-W5/0/0/0')
    #rsl = CFRSL('/home/magi/Development/src/PyCF/pelt.dat')
    #area.rsl_locations(rsl)
    area.coordsystem()
    area.printinfo(ts)
    plot.close()

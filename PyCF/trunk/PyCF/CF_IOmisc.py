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

"""Miscellaneous I/O operations."""

__all__ = ['CFreadlines','CFcontours','CFTimeSeries','CFEIStemp','CFEpoch']

import numpy, math

def CFreadlines(fobject, comment='#'):
    """Strip files from comments.

    fobject: file object to be read
    comment: string indicating comment.
    """

    lines = []
    for l in fobject.readlines():
        # ignore comments and empty lines
        l = l.strip()
        pos = l.find(comment)
        if pos>-1:
            l = l[:pos]
        if len(l)==0:
            continue
        lines.append(l)
    return lines

class CFcontours(list):
    """Read an ASCII file containing contours and store them in a list."""

    def __init__(self,fobject):
        """initialise.

        fobject: file object to be read
        """

        list.__init__(self)

        # read data
        cur=None
        for l in CFreadlines(fobject):
            pos = l.find(':')
            if pos>-1:
                cur = float(l[:pos])
                vertices = []
                self.append({'val':cur,'vert':vertices})
                continue
            if cur==None:
                raise RuntimeError, 'No contour value selected'
            values = l.split()
            try:
                values = [float(values[0]), float(values[1])]
            except:
                raise RuntimeError, 'Cannot read line:\n%s'%l
            vertices.append(values)

class CFTimeSeries(object):
    """Handling time series."""

    def __init__(self,fname,sep=None,timescale = 0.001):
        """Initialise

        fname: name of file to read from.
        sep:   separator (whitespace if set to none
        timescale: scale time"""

        infile = open(fname,'r')
        self.timescale = 0.001
        stime = []
        sdata = []
        self.numt = 0
        self.__index = 0
        for line in infile.readlines():
            pos = line.find('#')
            if pos>-1:
                line = line[:pos]
            line = line.strip()
            if len(line) == 0:
                continue
            if sep == None:
                l = line.split()
            else:
                l = line.split(sep)
            stime.append(float(l[0])*self.timescale)
            data = []
            for d in l[1:]:
                d = float(d)
                data.append(d)
            sdata.append(data)
            self.numt = self.numt + 1
        self.time = numpy.array(stime)
        self.data = numpy.array(sdata)

    def get_index(self,time):
        """Find index i so that t[i] <= time < t[i+1].

        time: time to find."""

        # BC
        if time <= self.time[0]:
            return 0
        if time >= self.time[self.numt-1]:
            return self.numt-1
        # check if time is around index
        if (time>=self.time[self.__index] and time<self.time[self.__index+1]):
            return self.__index
        # no, let's try the next interval
        self.__index = self.__index + 1
        if (time>=self.time[self.__index] and time<self.time[self.__index+1]):
            return self.__index
        # no, ok let's search for it using Newton
        lower = 1
        upper = self.numt-1
        while True:
            self.__index = lower + int((upper-lower)/2)
            if (time>=self.time[self.__index] and time<self.time[self.__index+1]):
                return self.__index
            if (time > self.time[self.__index]):
                lower = self.__index
            else:
                upper = self.__index
                
    def step(self,time):
        """Interpolate time series by stepping.

        time: time to get
        index: extract particular index, return array if index==None."""

        i = self.get_index(time)
        return self.data[i,:]

    def linear(self,time,index=None):
        """Linearly interpolating time series.

        time: time to get
        index: extract particular index, return array if index==None."""
        
        i = self.get_index(time)
        if i==0 or i==self.numt-1:
            data = self.data[i]
        else:
            d1 = self.data[i,:]
            d2 = self.data[i+1,:]
            factor = (time-self.time[i])/(self.time[i+1]-self.time[i])
            data = []
            for j in range(0,len(d1)):
                data.append(d1[j]+factor*(d2[j]-d1[j]))
        return numpy.array(data)


class CFEIStemp(CFTimeSeries):
    """Handling EIS temperature forcing."""

    def __init__(self,fname,lat=60.,sep=None,timescale = 0.001,temp_type='poly',lat0 = 44.95):
        """Initialise

        fname: name of file to read from.
        lat: latitude
        sep:   separator (whitespace if set to none
        timescale: scale time
        temp_type: select temp distribution function, can be poly (default) or exp
        lat0: parameter used for exponential function (default: 44.95)"""

        CFTimeSeries.__init__(self,fname,sep=sep,timescale=timescale)
        sdata = []
        for j in range(0,len(self.data)):
            val = self.data[j]
            t = 0.
            if temp_type == 'poly':
                l = 1.
                for i in range(0,len(val)):
                    t = t + l*val[i]
                    l = l*lat
            elif temp_type == 'exp':
                t = val[0]+val[1]*math.exp(val[2]*(lat-lat0))
            else:
                raise RuntimeError, 'No handle for temperature calculations type=\'%s\''%temp_type
            sdata.append([t])
        self.data = numpy.array(sdata)

class CFEpoch(object):
    """Handle epochs."""

    def __init__(self,fname):
        """Initialise epoch data from file.

        fname: name of file containing epoch data.

        file format: comments start with #
        each row contains 4 columns separate by commas: name, start time, end time and GMT RGB string."""

        self.data = []
        self.timescale = 0.001
        f = open(fname,'r')
        for l in CFreadlines(f):
            l = l.split(',')
            self.data.append({'name' : l[0], 'start':float(l[1]), 'end':float(l[2]),'colour':l[3].strip()})
        f.close()

        self.__current = 0

    def get_epoch(self,t):
        """Return the name of the epoch given a time.

        t: time in ka"""

        time = t/self.timescale

        # first try if current points to the right epoch
        if self.data[self.__current]['start']<=time and self.data[self.__current]['end']>=time:
            return self.__current

        for c in range(0,len(self.data)):
            if self.data[c]['start']<=time and self.data[c]['end']>=time:
                self.__current = c
                return self.__current

        return None

    def get_colour(self,t):
        """Return the RGB string associated with time.

        t: time in ka"""

        c = self.get_epoch(t)
        if c != None:
            return self.data[c]['colour']
        else:
            return '255/255/255'

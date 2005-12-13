#! /usr/bin/env python

"""Plot sediment maps"""

import PyGMT,PyCF,sys,os,string

# creating option parser
parser = PyCF.CFOptParser()
parser.var_options()
parser.profile_file(plist=True)
parser.time()
parser.region()
parser.plot()
opts = PyCF.CFOptions(parser,-2)

count = 0
numplots = 0
if opts.ntimes>1:
    count = count + 1
    numplots = opts.ntimes
if opts.nfiles>1:
    count = count + 1
    numplots = opts.nfiles
    
if count > 1:
    sys.stderr.write('Error, can only have either more than one time slice or more than one file!\n')
    sys.exit(1)

infile = opts.cffile()

time = opts.times(infile)
# get number of plots
deltax = 1.
deltay = 1.
sizex = opts.options.width+deltax
sizey = infile.aspect_ratio*opts.options.width+deltay

plot=None
if numplots > 1:
    pass
else:
    plot = opts.plot()
    area = PyCF.CFArea(plot,infile,pos=[0.,3.],size=sizex-deltax)
    if opts.options.land:
        area.land(time,illuminate=opts.options.illuminate)

    init_seds = infile.getvar("seds3").getGMTgrid(0)
    seds1 = infile.getvar("seds1").getGMTgrid(time)
    seds2 = infile.getvar("seds2").getGMTgrid(time)
    seds3 = infile.getvar("seds3").getGMTgrid(time)

    seds3.data = seds3.data+ seds1.data + seds2.data - init_seds.data

    area.raw_image(infile,time,seds3,os.path.join(PyCF.CFdatadir,"erosion.cpt"),isvelogrid=True,clip = opts.options.clip,mono=opts.options.mono,illuminate=opts.options.illuminate)

    area.coastline()
    try:
        thk = infile.getvar('thk')
        area.contour(thk,[0.1],'-W2/0/0/0',time)
    except:
        pass
    area.coordsystem()
    area.printinfo(time)
    if opts.options.profname !=None:
        i=0
        for pn in opts.options.profname:
            pdata = PyCF.CFprofile(infile,interval=opts.options.interval)
            pdata.coords_file(pn,opts.options.prof_is_projected)
            area.profile(args='-W5/0/0/0',prof=pdata,slabel=string.ascii_uppercase[i])
            i=i+1    
    if opts.options.dolegend:
        PyGMT.colourkey(area,var.colourmap.cptfile,title=var.long_name,pos=[0,-2])
    
plot.close()

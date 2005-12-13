#! /usr/bin/env python

"""Plot sediment profiles."""

import PyGMT,PyCF,sys,os, Numeric, Polygon

def munch_erosion(cfprofile,epochs,etime,yrange=None):
    """Produce erosion surfaces.

    cfprofile: CF file containing erosion.
    epochs: time slices at which to create new polynoms.
    etime: last time slice"""

    time = [0,etime]
    times = cfprofile.time(time)

    seds1 = cfprofile.getprofile('seds1').getProfileTS(time=time).data
    seds2 = cfprofile.getprofile('seds2').getProfileTS(time=time).data
    seds3 = cfprofile.getprofile('seds3').getProfileTS(time=time).data
    
    init_seds = seds1[:,0] + seds2[:,0] + seds3[:,0]
    hb = cfprofile.getprofile('erosion').getProfile(time=etime) - init_seds

    upper_x = Numeric.zeros(len(init_seds)+2,Numeric.Float)
    upper_x[:-2] = Numeric.array(cfprofile.xvalues)
    upper_x[-2] = upper_x[-3]
    upper_x[-1] = upper_x[0]
    lower_x = Numeric.zeros(len(init_seds)+2,Numeric.Float)
    lower_x[2:] = Numeric.array(cfprofile.xvalues)
    lower_x[0] = lower_x[-1]
    lower_x[1] = lower_x[2]

    # loop over time
    for j in range(0,len(seds1[0,:])):
        seds1[:,j] = seds1[:,j] + seds2[:,j] + seds3[:,j] - init_seds[:]

    minseds = PyGMT.round_down(min(hb))
    if yrange!=None:
        minseds = min(minseds,yrange[0])
    maxseds = PyGMT.round_up(max(Numeric.ravel(seds1)))

    upper_y = Numeric.zeros(len(init_seds)+2,Numeric.Float)
    upper_y[-2:] = maxseds
    lower_y = Numeric.zeros(len(init_seds)+2,Numeric.Float)
    lower_y[:2] = minseds

    polygons = []
    colours = []
    stages = []
    # initial sediment distribution
    polygons.append( Polygon.Polygon( Numeric.transpose(Numeric.array([lower_x,lower_y])) ) )
    colours.append('255/255/0')
    stages.append('orig seds')
    current_epoch = epochs.get_epoch(times[0])
    for j in range(0,len(seds1[0,:])):
        e = epochs.get_epoch(times[j])
        if e != current_epoch:
            # add a new stage
            lower_y[2:] = seds1[:,j]
            polygons.append( Polygon.Polygon( Numeric.transpose(Numeric.array([lower_x,lower_y])) ) )
            current_epoch = e
            colours.append(epochs.get_colour(times[j]))
            stages.append(current_epoch)
        
        upper_y[:-2] = seds1[:,j]
        upoly = Polygon.Polygon( Numeric.transpose(Numeric.array([upper_x,upper_y])) )
        for p in range(0,len(polygons)):
            polygons[p] = polygons[p] - upoly

    # hard bedrock erosion
    lower_y[2:] = hb[:]
    polygons.insert(0, Polygon.Polygon( Numeric.transpose(Numeric.array([lower_x,lower_y])) ) )
    colours.insert(0, '202/119/74')
    stages.insert(0, 'Bedrock')

    coords = []
    for p in range(0,len(polygons)):
        c = []
        for i in range(0,len(polygons[p])):
            c.append({'x':[],'y':[]})
            for j in range(0,len(polygons[p][i])):
                c[i]['x'].append(polygons[p][i][j][0])
                c[i]['y'].append(polygons[p][i][j][1])
        coords.append(c)
        
    return (coords,colours,stages)
    
if __name__ == '__main__':
    parser = PyCF.CFOptParser()
    parser.width=15.
    parser.time()
    parser.epoch()
    parser.profile_file()
    parser.region1d()
    parser.plot()
    opts = PyCF.CFOptions(parser,2)

    try:
        epoch = PyCF.CFEpoch(opts.options.epoch)
    except:
        parser.error('Could not load file containing glacial stages: %s'%opts.options.epoch)
    infile = opts.cfprofile()

    time = opts.times(infile,0)

    # load data
    (polygons,colours,stages) = munch_erosion(infile,epoch,time,yrange=opts.options.yrange)

    #print colours,stages

    # setup plot
    plot = opts.plot()
    plot.defaults['LABEL_FONT_SIZE']='12p'
    plot.defaults['ANOT_FONT_SIZE']='10p'
    bigarea = PyGMT.AreaXY(plot,size=opts.papersize)

    key = PyGMT.KeyArea(bigarea,pos=[0.,0.],size=[opts.options.width,2])
    key.num = [4,3]
    for i in range(0,len(stages)):
        key.plot_box(stages[i],'%s'%colours[i])

    area = PyCF.CFProfileMArea(bigarea,pos=[0.,3.5],size=[opts.options.width,opts.options.width/5.])


    prof_seds = infile.getprofile('seds1')
    seds_area = area.empty()
    seds_area.ylabel = 'sediments [m]'
    for p in range(len(polygons)-1,-1,-1):
        for i in (range(0,len(polygons[p]))):
            seds_area.line('-G%s  -W1/0/0/0'%(colours[p]),polygons[p][i]['x'],polygons[p][i]['y'])
    if opts.options.yrange != None:
        seds_area.ll[1] = opts.options.yrange[0]
        seds_area.ur[1] = opts.options.yrange[1]

    # plot ice surface
    prof_is = infile.getprofile('is')
    area.newprof(prof_is,time)

    area.finalise()
    area.coordsystem()

    plot.close()

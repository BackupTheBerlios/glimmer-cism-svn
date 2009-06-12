#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import os.path
import math
from matplotlib import pyplot
import matplotlib.lines
import matplotlib.patches
import matplotlib.figure
import matplotlib.font_manager
import numpy
from getopt import gnu_getopt

import intercomparison


def plotAggregatedFlowlines(axis, meanFlowline, stdevFlowline, dataMember, color, labelPrefix):
    xs = meanFlowline.getPointLocations()
    means = meanFlowline.getDependantVariable(dataMember)
    stdevs = stdevFlowline.getDependantVariable(dataMember)

    #Make sure the color is a tuple, not a list
    color = tuple(color)

    #Plot the mean
    axis.plot(xs, means, ":", color=color, label=labelPrefix + " Mean")
    #Shade the standard deviation around the mean
    #Set up a list of points on the polygon to shade.
    #First, we get a list of the upper and lower lines
    meanPlusSd =  [m+s for m,s in zip(means,stdevs)]
    meanMinusSd = [m-s for m,s in zip(means,stdevs)]
    
    #axis.plot(xs,meanPlusSd)
    #axis.plot(xs,meanMinusSd)
    
    #Now, we set up the list of points that form the polygon
    #This funky list reversing is so that we go around the polygon counter-clockwise;
    #that is, we specify in one list the upper bound going left-to-right, then the
    #lower bound going right-to-left.
    polyX = xs + list(reversed(xs))
    polyY = meanPlusSd + list(reversed(meanMinusSd))
    #Plot the polygon as the same color but at only 1/4 opacity
    axis.fill(polyX, polyY, facecolor=color,edgecolor=color,alpha=.25, label=labelPrefix + " Std. Dev.")

def maxObservedError(meanFlowline, stdevFlowline, myFlowline):
    means = meanFlowline.getDependantVariable(0)
    stdevs = stdevFlowline.getDependantVariable(0)
    computed = myFlowline.getDependantVariable(0)
    currentMax = -float("inf")
    matchingStdevError = None
    for mean, stdev, val in zip(means, stdevs, computed):
        err = (val-mean)/mean
        if err > currentMax:
            currentMax = err
            matchingStdevError = stdev/mean
    return currentMax, matchingStdevError
            

#Returns a tuple with two tuples suitable for creating a legend.  The first tuple contains
#the lines and patches, the second contains the names
def createPlot(experiment, domainSizeKm, fig, subplotNum, notFullStokesModelType):
    axis = fig.add_subplot(2,3,subplotNum)
    axis.set_title(str(domainSizeKm) + " km", size="medium")
    #Read the data that came from the Glimmer run
    glimFileName = glimPrefix + experiment + "%03d"%domainSizeKm + ".txt"
    glimFlowline = intercomparison.grabFlowline(glimFileName)

    axis.plot(glimFlowline.getPointLocations(), glimFlowline.getDependantVariable(0), color=(0,0,0))


    for isFullStokes in [True, False]:
        if isFullStokes:
            modelType = "full-stokes"
        else:
            modelType = notFullStokesModelType
        models = intercomparison.getDataFiles("ismip_all", experiment, domainSizeKm, modelType)
        
        flowlines = []
        for f in models:
            try:
                fl = intercomparison.grabFlowline(f)
                fl.fixRange()
                flowlines.append(fl)
            except Exception, e:
                print f, "failed to load:", e
                raise

        #Compute the mean and standard deviations of the experiments
        mean, stdev = intercomparison.aggregateExperimentFlowlines(flowlines, glimFlowline.getPointLocations())
        
        if not isFullStokes:
            myErr,oneStdev = maxObservedError(mean, stdev, glimFlowline)
            print '\t'.join([str(domainSizeKm) + " km", str(myErr), str(oneStdev)])

        #Plot the mean and std. dev.
        if isFullStokes:
            plotAggregatedFlowlines(axis, mean, stdev, 0, (1,0,0), "Full Stokes")
        else:
            plotAggregatedFlowlines(axis, mean, stdev, 0, (0,0,1), "First Order")

        #Modify the axes tick lables to have smaller fonts
        for tick in axis.xaxis.get_major_ticks():
            tick.label1.set_fontsize("xx-small")
        for tick in axis.yaxis.get_major_ticks():
            tick.label1.set_fontsize("xx-small")
        

if __name__ == "__main__":
    optlist, args = gnu_getopt(sys.argv[1:], 
                               intercomparison.ExperimentOptionsShort, 
                               intercomparison.ExperimentOptionsLong + ["lmla", "prefix=", "subtitle="])
    optdict = dict(optlist)
   
    experiments, domainSizes = intercomparison.getExperimentsToRun(optlist)

    #Convert domain sizes from meters to km
    domainSizes = [d/1000 for d in domainSizes]

    if "--lmla" in optdict:
        notFullStokesModelType = "lmla"
    else:
        notFullStokesModelType = "partial-stokes"

    if "--prefix" in optdict:
        glimPrefix = optdict["--prefix"]
    else:
        glimPrefix = "glm1"

    if "--subtitle" in optdict:
        subtitle = optdict["--subtitle"]
    else:
        subtitle = None

    for experiment in experiments:
        print "ISMIP-HOM",experiment.upper()

        fig = pyplot.figure(subplotpars=matplotlib.figure.SubplotParams(top=.85,bottom=.15))
        for i, domainSize in enumerate(domainSizes):
            createPlot(experiment, domainSize, fig, i+1, notFullStokesModelType)

        #Create the legend!  This is overly complicated because I want the legend to be
        #in multiple columns at the bottom of the visual, and that's impossible with my
        #version of matplotlib (it's apparently coming though?)
        #So, I am creating a different legend for each column
        #l.draw_frame(False) turns off the borders.
        l = fig.legend([matplotlib.lines.Line2D([1,2],[1,2],color=(0,0,0))],
                   ['Model Output'] ,
                   loc=(.1, .05),prop=matplotlib.font_manager.FontProperties(size="x-small"))
        l.draw_frame(False)
        l = fig.legend([matplotlib.lines.Line2D([1,2],[1,2],ls=':',color=(1,0,0)),
                    matplotlib.patches.Patch(edgecolor=None, facecolor=(1,0,0), alpha=.25)],
                    ['Full Stokes Mean', 'Full Stokes Std. Dev.'], 
                   loc=(.3, .02),prop=matplotlib.font_manager.FontProperties(size="x-small"))
        l.draw_frame(False)
        l = fig.legend([matplotlib.lines.Line2D([1,2],[1,2],ls=':',color=(0,0,1)),
                    matplotlib.patches.Patch(edgecolor=None, facecolor=(0,0,1), alpha=.25)],
                    ['First Order Mean', 'First Order Std. Dev.'],
                   loc=(.55, .02),prop=matplotlib.font_manager.FontProperties(size="x-small"))
        l.draw_frame(False)

        #Add an overall title to the figure
        fig.text(.5, .92, "ISMIP-HOM Experiment " + experiment.upper(), horizontalalignment='center', size='large')
        if subtitle:
            fig.text(.5, .89, subtitle, horizontalalignment='center', size='small')
        
        #Add one axis label per side (the labels are the same because of the small multiple format,
        #so we only want to repeat them once in the whole figure!)
        fig.text(.5, .1, "Normalized X coordinate", horizontalalignment='center',size='small')
        fig.text(.06, .5, "Ice Speed (m/a)", rotation="vertical", verticalalignment='center') 
        #Save the figure
        filename = "ISMIP-HOM-" + experiment.upper() + "-" + glimPrefix + ".png" 
        pyplot.savefig(filename)


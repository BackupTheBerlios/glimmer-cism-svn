#!/usr/bin/env python
#
# Magnus Hagdorn
#
# process results from verif tests

import Numeric, Scientific.IO.NetCDF, pygsl, sys
from PyGMT import round_up,round_down
from pygsl import histogram

def process(fname):
    """Process verif output file.

    fname: name of file

    return: tuple of (exp_name, solver, dx, dt, dome_e, max_e, min_e, mean_e, sd_e)"""

    # extract errors
    ncf = Scientific.IO.NetCDF.NetCDFFile(fname)
    diff = ncf.variables['thke'][-1,:,:] - ncf.variables['thk'][-1,:,:]
    centre = (Numeric.shape(diff)[0]-1)/2
    dome_e = diff[centre,centre]
    diff = Numeric.ravel(diff)
    max_e = max(diff)
    min_e = min(diff)
    hist = histogram.histogram(100)
    hist.set_ranges_uniform(round_down(min_e),round_up(max_e))
    for e in diff.tolist():
        hist.increment(e)
    mean_e = hist.mean()
    sd_e = hist.sigma()

    config = ncf.title.split(',')
    exp_name = config[0][-1]
    solver = config[1].strip()
    dx = float(config[2].strip()[:-2])
    dt = float(config[3].strip()[:-1])
    
    ncf.close()

    return (exp_name, solver, dx, dt, dome_e, max_e, min_e, mean_e, sd_e)


if __name__ == '__main__':

    # create output file
    out = open('verif_results.data','w')
    out.write('#exp\tsolver\tdx\t\tdt\t\te_dome\t\tmax_e\t\tmean_e\t\tsd_e\n')
    for f in sys.argv[1:]:
        (exp_name, solver, dx, dt, e_dome, max_e, min_e, mean_e, sd_e) = process(f)
        out.write("%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\n"%(exp_name,solver,dx,dt,e_dome,max_e,mean_e,sd_e))
    out.close()
    
    

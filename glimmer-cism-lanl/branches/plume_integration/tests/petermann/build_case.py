#!/usr/bin/python

#  This script is used to build the following:
#         namelist file for plume model (if required)
#         netcdf input file for glimmer-cism
#         config file for glimmer-cism

# It works by taking values from an input file, which is itself
# a python module, and overwriting the default values with the
# provided values.

import sys
import subprocess

class FortranConversionException(Exception):
    def __init__(self, msg, param_name):
        Exception.__init__(self,msg)
        self.param_name = param_name

def fortran_style(key, val):
    if (type(val) == type(1)):
        return '%d' % val
    elif (type(val) == type(1.0)):
        return 'd'.join(('%e' % val).split('e'))
    elif (type(val) == type(True)):
        return val and 'T' or 'F'
    elif (type(val) == type('happy')):
          return val
    elif val is None:
        raise FortranConversionException(key)
    else:
        raise FortranConversionException(key)

class PlumeNamelist(object):

    def __init__(self):

        # Default values are defined here.
        
        # Where it wouldn't really make sense to have a default value
        # a value of None is assigned.  This must be updated by a non-None
        # value before writing out then namelist contents
        
        self.vals = {'mixlayer' : False,
                     'in_glimmer' : True,
                     'restart' : False,
                     'frazil' : False,
                     'nonlin' : True,
                     'horturb' : False,
                     'entrain' : True,
                     'entype' : 1,
                     'basmelt' : True,
                     'rholinear' : True,
                     'thermobar' : False,
                     'intrace' : False,
                     'vardrag' : False,
                     'topedit' : False,
                     'tangle' : False, 
                     'negfrz' : False,
                     'use_min_plume_thickness' : True,
                     'tottim'  : 0.0,
                     'outtim'  : 0.5,
                     'labtim'  : 0.25,
                     'snottim' : 0.25,
                     'lnottim' : 1.0,
                     'dt1'     : None,
                     'm_grid' : None,
                     'n_grid' : None,          
                     'hx' : None,      
                     'hy' : None,
                     'gldep' : None,
                     'ifdep' : None,
                     'wcdep' : None,
                     'plume_min_thickness' : None,
                     'bsmoothit' : 0,
                     'salttop' : 34.5,
                     'saltbot' : 34.5,
                     'temptop' : -1.000,
                     'tempbot' : -1.000,
                     'phi' : 0.0,
                     'ah' : 0.0,   
                     'kh' : 0.0,
                     'cdb' : 2.5e-3,    
                     'cl' : 1.775e-2,   
                     'ef' : 5.0e-1,     
                     'nus' : +1.0,     
                     'nbar' : 1.0e3,
                     'nice' : 10,
                     'seedtype' : 2,
                     'cseedfix' : 1.0e-7,
                     'cinffix'  : 0.0
                     }
                     
    def update_vals(self, specificVals):
       #next we overwrite default values with any values provided

       self.vals.update(specificVals)

    def produce_namelist_contents(self):
        
        # return a string which can be used as the contents of
        # a plume namelist file
        
        lines = []
        vItems = self.vals.items()
        vItems.sort()
        for (k,v) in vItems:
            try:
                lines.append(', %s = %s\n' % (k,fortran_style(k,v)))
            except FortranConversionException:
                raise Exception('%s \n %s \n %s' % 
                                ('Could not write dict to Fortran namelist.',
                                 'Missing definitions of:',
                                ' '.join([k for (k,v) in
                                          self.vals.items() if (v is None)])))
        return '&plume_nml\n %s /' % ' '.join(lines)


class GCConfig(object):

    def __init__(self):
        
        self.vals = { 'parameters' : { 'geothermal' : -42.0e-3,
                                       'default_flwa' : 4.6e-18,
                                       'flow_factor' : 1,
                                       #'ice_limit' : 10.0,
                                       #'marine_limit' : 0.0
                                       #'calving_fraction' : 0.0,
                                       #'hydro_time' : 0.0,
                                       # 'basal_tract' : ,
                                       # 'basal_tract_const' : ,
                                       'log_level' : 6,
                                       },
                      'Petermann shelf' : {  'air_temperature' : -5.0,
                                             'accumulation_rate' : 10.0,
                                             'eustatic_sea_level' : 0.0 },
                      'options' : {'flow_law' : 1,
                                   'evolution' : 3,
                                   'temperature' : 0,
                                   'vertical_integration' : 1,
                                   'marine_margin' : 0,
                                   'topo_is_relaxed' : 1,
                                   'slip_coeff' : 1,
                                   'sliding_law' : 4,
                                   'stress_calc' : 2,
                                   'periodic_ew' : 0,
                                   'periodic_ns' : 0,
                                   'hotstart' : 0,
                                   'basal_water' : 2,
                                   'which_bmlt' : 1},
                      'grid' : { 'sigma_builtin' : 1,
                                 'upn' : None,
                                 'ewn' : None,
                                 'nsn' : None,
                                 'dew' : None,
                                 'dns' : None },

                      'ho_options' : { 'which_ho_sparse_fallback' : -1,
                                       'basal_stress_input' : 3,
                                       'which_ho_resid' : 0,
                                       'which_ho_babc' : 9,
                                       'guess_specified' : 1,
                                       'which_ho_sparse' : 0,
                                       'diagnostic_scheme' : 3,
                                       'include_thin_ice' : 0,
                                       'which_ho_efvs' : 0},
                      'CF default' : { 'comment' : None,
                                       'title' : None,
                                       'institution' : 'NYU'},
                      'CF input' : { 'name' : None,
                                     'time' : None },
                      'CF output' : { 'variables' : ' '.join(['lsurf','usurf',
                                                              'thk','bmlt',
                                                              'acab',
                                                              'uvelhom',
                                                              'vvelhom',
                                                              'uvelhom_srf',
                                                              'uvelhom_bas',
                                                              'vvelhom_srf',
                                                              'vvelhom_bas',
                                                              'thkmask','topg',
                                                              'kinbcmask',
                                                              'beta','btrc']),
                                      'frequency' : None,
                                      'name' : None },
                      'time' : { 'tstart' : 0.0,
                                 'tend' : None,
                                 'dt' : None,
                                 'niso' : 1.0,
                                 'ntem' : 1.0,
                                 'nvel' : 1.0 },
                      'plume' : { 'plume_nl_file' : None,
                                  'plume_output_file' : None,
                                  'suppress_ascii_output' : True,
                                  'suppress_logging' : False,
                                  'plume_output_prefix' : None,
                                  'plume_output_dir' : './',
                                  'plume_write_all_states' : False,
                                  'plume_min_spinup_time' : 5.0,
                                  'plume_min_subcycle_time' : 0.5,
                                  'plume_steadiness_tol' : 1.0e-6,
                                  'plume_imin' : None,
                                  'plume_imax' : None,
                                  'plume_kmin' : None,
                                  'plume_kmax' : None }
                      }
        
    def update_vals(self, specificVals):

        for (k,v) in specificVals.items():
            if (not (k in self.vals.keys())):
                raise Exception ('Unknown section name: %s' % k)
            self.vals[k].update(v)
        
    def produce_config_file(self):

        sections = []

        vItems = self.vals.items()
        vItems.sort()
        
        for (k,vdict) in vItems:
            try:
                vdictItems = vdict.items()
                vdictItems.sort()
                section_vals = [' %s = %s' %
                                (k_inner,fortran_style(k_inner,v))
                                for (k_inner,v) in vdictItems]
                
            except FortranConversionException, (e):
                raise Exception('Failed to convert %s' %
                                e.param_name)
            
            sections.append('[%s]\n%s\n' % (k, '\n'.join(section_vals)))
                        
        return '\n'.join(sections)
            
def main(config_filename):

    newVals = {}
    
    execfile(config_filename, globals(), newVals)

    req_vars = ['plume_vals',
                'gc_vals',
                'plume_nl_filename',
                'gc_config_filename',
                'nc_gen_input_args',
                ]

    for req_var in req_vars:
        if (not req_var in newVals):
            raise Exception('Did not find %s defined in %s' % (req_var,
                                                               config_filename))
    pnl = PlumeNamelist()    
    pnl.update_vals(newVals['plume_vals'])

    gc_config = GCConfig()
    gc_config.update_vals(newVals['gc_vals'])

    pnl_name = newVals['plume_nl_filename']
    gc_config_name = newVals['gc_config_filename']
    nc_gen_input_args = newVals['nc_gen_input_args']
    
    f = open(pnl_name,'w')
    try:
        f.write(pnl.produce_namelist_contents())
    finally:
        f.close()

    g = open(gc_config_name,'w')
    try:
        g.write(gc_config.produce_config_file())
    finally:
        g.close()

    cmd = ['nc_gen_input']
    cmd.extend(nc_gen_input_args)
    cmd = [fortran_style('nc_gen_input', c) for c in cmd]
    
    retcode = subprocess.call(cmd)
    if (retcode != 0):
        raise Exception('Error running nc_gen_input')

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        raise Exception('Need to call with config file as last argument')
    
    config_filename = sys.argv[-1] #last argument
    main(config_filename)

    
    

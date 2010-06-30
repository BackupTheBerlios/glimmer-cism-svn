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
    def __init__(self, param_name):
        Exception.__init__(self,'Error converting key %s' % param_name)
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
                     'horturb' : True,
                     'entrain' : True,
                     'entype' : 1,
                     'basmelt' : True,
                     'rholinear' : True,
                     'thermobar' : False,
                     'intrace' : False,
                     'vardrag' : False,
                     'topedit' : False,
                     'tangle' : None, 
                     'negfrz' : False,
                     'use_min_plume_thickness' : True,
                     'tottim'  : 0.0,
                     'outtim'  : 1.0,  #every day
                     'labtim'  : 0.0,
                     'snottim' : 0.0,
                     'lnottim' : 0.0,
                     'dt1'     : None,
                     'm_grid' : None,
                     'n_grid' : None,          
                     'hx' : None,      
                     'hy' : None,
                     'gldep' : None,
                     'ifdep' : None,
                     'wcdep' : None,
                     'plume_min_thickness' : None,
                     'bsmoothit' : 0,   #iterations of smoothing
                     'salttop' : None,
                     'saltbot' : None,
                     'temptop' : None,
                     'tempbot' : None,
                     'phi' : None,
                     'ah' : None,   
                     'kh' : None,
                     'cdb' : 2.5e-3,  #bottom drag coeff    
                     'cl' : 1.775e-2, # entrainment coeff
                     'ef' : 5.0e-1,   # entrainment factor
                     'context' : '""',  #NB: don't remove this
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
                lines.append(' %s = %s\n' % (k,fortran_style(k,v)))
            except FortranConversionException:
                raise Exception('%s \n %s \n %s' % 
                                ('Could not write dict to Fortran namelist.',
                                 'Missing definitions of:',
                                ' '.join([k for (k,v) in
                                          self.vals.items() if (v is None)])))
        return '&plume_nml\n %s /' % ' ,'.join(lines)


class GCConfig(object):

    def __init__(self):
        
        self.vals = { 'parameters' : { 'geothermal' :   0.0,      # geothermal heat flux 
                                                                  # model%paramets%geot
                                       'default_flwa' : None,     # Glen's law A to use in isothermal case 
                                                                  # model%paramets%default_flwa
                                       'flow_factor' : 1,         # enhancement factor for glen's A 
                                                                  # model%paramets%flow_factor
                                       'ice_limit' : None,        # minimum thickness for running ice dynamics
                                                                  # model%numerics%thklim
                                       'marine_limit' : 0.0,      # When to chop off marine ice
                                                                  # NB: not used if marine_margin = 0 is used
                                                                  # model%numerics%mlimit
                                       'calving_fraction' : 0.0,  # fraction of ice to remove from floating ice
                                                                  # when using marine_margin = 3
                                                                  # model%numerics%calving_fraction
                                       'hydro_time' : 0.0,        # time constant for basal hydrology
                                                                  # model%paramets%hydtim
                                       # 'basal_tract' :          # 5-value parameters - obsolete
                                       'basal_tract_const' : 0.0, # basal_tract_const used over whole ice base
                                                                  # model%paramets%btrac_const
                                       # 'basal_tract_max' :      # model%paramets%btrac_max
                                       # 'basal_tract_slope' :    # dependence of basal traction on basal water depth
                                                                  # model%paramets%btrac_slope
                                       # 'stressin':              # initial backstress in areas of positive thickness
                                                                  # NB: only used in case marine_margin = 5
                                                                  # model%climate%stressin
                                       # 'stressout' :            # initial backstress assigned to other areas
                                                                  # model%climate%stressout
                                                                  # NB: only used in case slip_coeff = 5
                                       # 'sliding_constant' : =   # model%climate%slidconst
                                       'log_level' : 6,
                                       },
                      'Petermann shelf' : {  'air_temperature' : None,     # Temp assigned to ice in isothermal case
                                             'accumulation_rate' : None,   # In meters per year
                                             'eustatic_sea_level' : 0.0 }, # Height of sea level relative to initial height
                      'options' : {'flow_law' : None,     # flow_law = 2 means constant A
                                                       #          = 0 means calculate A from temperature
                                                       # model%options%whichflwa
                                   'evolution' : 3,    # evolution = 0 means pseudo-diffusion
                                                       #           = 1 means ADI scheme
                                                       #           = 2 means iterated diffusion
                                                       #           = 3 means LANL incrementral remapping method
                                                       #           = 4 upwind advection
                                                       # model%options%whichevol
                                   'temperature' : None,            # temperature = 1 means full 3D thermal evolution
                                                                 #             = 0 means set to air temperature
                                                                 # model%options%whichtemp
                                   'vertical_integration' : 1,  #vertical_integration = 1 constrained to obey kinematic BC
                                                                #NB: only used in full-temperature cases
                                                                # model%options%whichwvel
                                   'marine_margin' : 0,         # marine_margin 
                                                                # = 0 means ignore marine margin (no chopping)
                                                                # model%options%whichmarn
                                   'topo_is_relaxed' : 1,       # topo_is_relaxed = 1 means the init.l topography is relaxed
                                                                # model%options%whichrelaxed

                                   'slip_coeff' : 1,   # slip_coeff = 0 means set equal to zero everywhere
                                                       #          = 1 means basal traction is constant (basal_tract_const)
                                                       # = 2 means Set to (non--zero) constant where where temperature
                                                       #  is at pressure melting point of ice, otherwise to zero
                                                       #= 3 means function of basal water depth 
                                                       # model%options%whichbtrc
                                   'periodic_ew' : 0,  # model%options%periodic_ew
                                   'periodic_ns' : 0,  # model%options%periodic_ns
	                           'x_invariant' : None,  # = 1 means variables don't change in x direction
                                                          # model%options%x_invariant
                                   'diagnostic_run' : 0,  # = 1 makes glide stop after diagnosing velocities
                                                          # model%options%diagnostic_run
                                   'hotstart' : None,   # are we doing a restart (1 = yes, 0 = no)
                                                        # model%options%hotstart
                                   'basal_water' : 3,   # basal_water = 0 means calc from local basal water balance
                                                        #  = 1 means compute basal water flux then find depth via cal
                                                        #  = 2 means no basal water
                                                        # model%options%whichbwat

                                   'use_plume' : None,  # use_plume = 0 means use usual Glimmer-CISM method to calc bmlt
                                                        #  = 1 means use plume model to calculate bmlt above floating ice
                                                        # model%options%use_plume
                                   # 'ioparams' :                #i/o parameters file
                                   },
                      'grid' : { 'sigma_builtin' : 1,   # sigma_builtin = 0 means use the default sigma levels 
                                                        #	        = 1 means evenly-spaced levels
                                                        # model%options%which_sigma_builtin
                                 'upn' : None,      # number of vertical levels : model%general%upn
                                 'ewn' : None,      # number of grid nodes in east-west direction : model%general%ewn
                                 'nsn' : None,      # number of grid nodes in north-south direction : model%general%nsn
                                 'dew' : None,      # east-west grid spacing : model%numerics%dew
                                 'dns' : None       # north-south grid spacing : model%numerics%dns
                                 },
                      #'sigma' : { # 'sigma_levels' : 0.0 0.5 1.0    # model%numerics%sigma
                                   # 'sigma_file' : }
                      'ho_options' : { 'which_ho_sparse_fallback' : -1, # fallback solver method:
                                                                        # which_ho_sparse_fallback = -1 means no fallback
                                                                        # model%options%which_ho_sparse_fallback
                                       'basal_stress_input' : 2,    # how to compute beta:
                                                                    # = 2 means basal traction is given directly
                                                                    #NB: this all doesn't matter if using
                                                                    #    which_ho_babc = 6 (so beta ~= 0.0 everywhere)
                                                                    # model%options%which_ho_beta_in
                                       'basal_stress_type' : 0,      # basal_stress_type = 0 means linear
                                                                    #                     1 means plastic
                                                                    # NB: only used in velo_hom_pattyn and veloc2
                                                                    # model%options%which_ho_bstess
                                       'which_ho_source' : 0,   # how to compute source term for an ice shelf: 
                                                                # = 0 means vertically averaged
                                                                # = 1 means pressure dependent on depth
                                                                # = 2 means shelf front disabled  
                                                                # NB: only used inside Pattyn's veloc2 subroutine
                                                                # model%options%which_ho_source
                                       'which_disp' : None,         # dissipation
                                                                    # model%options%which_disp
                                       'which_ho_resid' : 0,  # method of calculating residual: 
                                                              # = 0 means use max value
                                                              # NB: only used in glam_velo_fordsiapstr
                                                              # model%options%which_ho_resid
                                       'which_bmelt' : 0,   # basal melting
                                                              # model%options%which_bmelt
                                       'which_ho_babc' : None,       # basal boundary condition: 
                                                                  # = 5 means simple ice-shelf
                                                                  # = 6 means floating ice everywhere, no traction
                                                                  # = 3 means circular ice-shelf
                                                                  # NB: only used when using Payne-Price diagnostic scheme
                                                                  # model%options%which_ho_babc
                                       'guess_specified' : 1,     # model%velocity_hom%is_velocity_valid
                                       'which_ho_sparse' : 0,     # which sparse solver to use:
                                                                  #  = 0 means biCG with incomplete LU precond. 
                                                                  #  = 1 means GMRES
                                                                  #  = 2 means UMF (?)
                                                                  #  = 3 means PARADISO (?)
                                                                  # model%options%which_ho_sparse
                                       'diagnostic_scheme' : 3,  # which higher-order diagnostic scheme to use:
                                                                 # = 3 means Payne-Price scheme
                                                                 # = 2 means Pattyn staggered ?            
                                                                 # model%options%which_ho_diagnostic
                                       'prognostic_scheme' : 0,   # which higher-order prognostic scheme to use
                                                                  # thickness evolution (only used in thick_nonline_evolve
                                                                  # and thick_lin_evolve, but not in incremental remapping)
                                                                  # model%options%which_ho_prognostic
                                       'include_thin_ice' : 0,   # whether or not to include thin ice in HO calculation
                                                                 # 1 means true
                                                                 # model%options%ho_include_thinice
                                       'which_ho_efvs' : 0     # ho effective viscosity 
                                                               # = 0 means calculate from strain rate
                                                               # model%options%which_ho_efvs
                                       },
                      'CF default' : { 'comment' : '',
                                       'title' : None,
                                       'institution' : 'NYU',
                                       'references' : ''
                                       },

                      'CF input' : { 'name' : None,     #name of netcdf file containing input data fields
                                     'time' : None      #which time slice to read input data from
                                     },
                      'CF output' : { 'variables' : ' '.join(['lsurf','usurf',
                                                              'thk','bmlt',
                                                              'acab',
                                                              'uvelhom',
                                                              'vvelhom',
                                                              'uvelhom_srf',
                                                              'vvelhom_srf',
                                                              'thkmask','topg',
                                                              'kinbcmask',
                                                              'beta','btrc',
                                                              'temp',
                                                              'tau_hom_xx','tau_hom_yy',
                                                              'tau_hom_xz','tau_hom_yz','tau_hom_xy']),
                                      ### NB: there is a (250) character limit on line length!!!

                                      # the following is a list of all possible output variables:
                                      # level lithoz x0 x1 y0 y1 acab acab_tavg age artm backstress
                                      # beta bheatflx bmlt bmlt_tavg btemp btrc bwat bwatflx calving 
                                      # diffu dusrfdtm eus flwa gl_ew gl_ns gline_flux iarea ivol
                                      # kinbcmask lat litho_temp lon lsurf relx slc soft surfvel tau_xz tau_yz
                                      # taux tauy temp thk thkmask topg ubas ubas_tavg uflx usurf 
                                      # uvel uvelhom vbas vbas_tavg velnormhom vflx vvel vvelhom wgrd wvel

                                      'frequency' : None,  # time in between writing state to output file (in years)
                                      'name' : None,       # name of output file
                                      'start' : None,
                                      'stop'  : None,
                                      },
                      'time' : { 'tstart' : None,   # start time of model run (years) : model%numerics%tstart
                                 'tend' : None,     # end time of model run (years) : model%numerics%tend
                                 'dt' : None,      # time step (years) : model%numerics%tinc
                                 'niso' : 1.0,     # isostasy dt factor
                                 'ntem' : 1.0,     # thermal dt factor : model%numerics%ntem
                                 'nvel' : 1.0,      # velocity dt factor : model%numerics%nvel
                                 # 'ndiag' :        # diagnostic frequency : model%numerics%ndiag
                                 # 'profile' :      # profile period : model%numerics%profile_period
                                 },
                      'plume' : { 'plume_nl_file' : None, # path to plume namelist file
                                  'plume_output_file' : None,   # netcdf file with plume output
                                  'suppress_ascii_output' : True,  #  suppress all old-style ASCII data output
                                  'suppress_logging' : False,    # suppress all screen and file output (logging)
                                  'plume_output_prefix' : None, # prefix to put on the old style ASCII output
                                  'plume_output_dir' : './',     # where to write the output files
                                  'plume_write_all_states' : False,   # option to write out all states (all timesteps)
                                                                      # NB it is very storage hungry

                                  'plume_min_spinup_time' : 5.0,    # minimum time to spinup the plume, in days
                                  'plume_min_subcycle_time' : 0.5,   # minimum subcycle time, in days
                                  'plume_steadiness_tol' : 1.0e-6,  # plume steadiness tolerance 
                                                                    # (max relative change in bmelt and speed)
                                  'plume_imax' : None,
                                  'plume_kmin' : None,
                                  'plume_kmax' : None,
                                  'plume_const_bmlt' : False,       # Apply a uniform melt rate under floating ice
                                  'plume_const_bmlt_rate' : 0.0     # At given rate in meters per year
                                  }
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
                print('\n\nFailed to convert %s -> %s\n\n' % (k,e.param_name))
                raise Exception('Failed to build case')
            
            sections.append('[%s]\n%s\n' % (k, '\n'.join(section_vals)))
                        
        return '\n'.join(sections)


def get_config(config_filename):

    configVals = {}
    
    execfile(config_filename, globals(), configVals)

    req_vars = ['plume_vals',
                'gc_vals',
                'plume_nl_filename',
                'gc_config_filename',
                'nc_gen_input_args',
                'nc_regrid_args',
                'input_style'
                ]

    for req_var in req_vars:
        if (not req_var in configVals):
            raise Exception('Did not find %s defined in %s' % (req_var,
                                                               config_filename))

    return configVals

def build_case(configVals):

    pnl = PlumeNamelist()    
    pnl.update_vals(configVals['plume_vals'])

    gc_config = GCConfig()
    gc_config.update_vals(configVals['gc_vals'])

    pnl_name = configVals['plume_nl_filename']
    gc_config_name = configVals['gc_config_filename']
    nc_gen_input_args = configVals['nc_gen_input_args']
    nc_regrid_args = configVals['nc_regrid_args']
    
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

    inputStyle = configVals['input_style']
    
    if (inputStyle == 'gen_input'):
        cmd = ['nc_gen_input']
        cmd.extend(nc_gen_input_args)
        cmd = [fortran_style('nc_gen_input', c) for c in cmd]
    elif (inputStyle == 'regrid'):
        cmd =['nc_regrid']
        cmd.extend(nc_regrid_args)
        cmd = [fortran_style('nc_regrid', c) for c in cmd]
    elif (inputStyle == 'given'):
        cmd = ['echo','Reading input from a given file']
#        raise Exception('No arguments to generate netcdf input')
    else:
        raise Exception('Invalid input_style: %s' % inputStyle)
    
    retcode = subprocess.call(cmd)
    if (retcode != 0):
        raise Exception('Error running:\n %s' % ' '.join(cmd))

def main():

    if (len(sys.argv) < 2):
        raise Exception('Need to call with config file as last argument')
    
    config_filename = sys.argv[-1] #last argument
    configVals = get_config(config_filename)

    build_case(configVals)
        

if __name__ == '__main__':
    main()
    
    

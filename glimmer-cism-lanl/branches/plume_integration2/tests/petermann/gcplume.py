import sys
import os
import subprocess
import pickle
import copy
import re

class defaultdict(dict):

    def __init__(self,defaultVal):
        dict.__init__(self)
        self._defaultVal = defaultVal

    def __getitem__(self,k):
        if (not(self.__contains__(k))):
            self[k] = copy.copy(self._defaultVal)
        return dict.__getitem__(self,k)            

class FortranConversionException(Exception):
    def __init__(self, param_name):
        Exception.__init__(self,'Error converting key %s' % param_name)
        self.param_name = param_name

def _fortran_style(key, val):
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
                     'restart_data_filename' : '',
                     'frazil' : False,
                     'nonlin' : True,
                     'horturb' : True,
                     'entrain' : True,
                     'entype' : 1,
                     'entype2' : 1,
                     'switch_entype_n' : 0,
                     'C_i' : 20.0,
                     'C_n' : 100.0,
                     'C_s' : 100.0,
                     'nk_m' : 0.45,
                     'nk_n' : 0.2,
                     'basmelt' : True,
                     'rholinear' : True,
                     'thermobar' : False,
                     'intrace' : False,
                     'vardrag' : False,
                     'topedit' : False,
                     'tangle' : False, 
                     'negfrz' : False,
                     'use_min_plume_thickness' : True,
                     'use_neutral_salinity' : False,
                     'depinffix' : 0.0,
                     'depinit' : 0.0,
                     'meltinf' : 0.0,
                     'sgd_type' : -1,
                     'sgd_flux' : 0.0,
                     'plume_southern_bc' : 0,       # 0 means d/dy = 0
                                                    # 1 means linear extrapolation
                                                    # 2 means quadratic extrapolation
                                                    # 3 means advective bc
                     'tottim'  : 0.0,
                     'labtim'  : 1.0,
                     'snottim' : 24.0*3600.0,  #once per day
                     'lnottim' : 100.0*365.25*86400.0, #once per 100 years
                     'dt1'     : None,
                     'm_grid' : None,
                     'n_grid' : None,
                     'namb' : 301,
                     'hx' : None,      
                     'hy' : None,
                     'gldep' : None,
                     'ifdep' : 0.0,
                     'wcdep' : 1000.0,
                     'plume_min_thickness' : 10.0,
                     'plume_max_thickness' : 100.0,
                     'u_star_offset' : 0.0,
                     'tidal_velocity' : 0.0,
                     'infloain' : 0.0,
                     'infloein' : 0.0,
                     'knfloain' : 0.0,
                     'knfloein' : 0.0,
                     'entr_time_const' : None,
                     'detrain_time_const' : None,
                     'bsmoothit' : 0,   #iterations of smoothing
                     'salttop' : None,
                     'saltbot' : None,
                     'temptop' : None,
                     'tempbot' : None,
                     'n_amb_ctl_pt' : 0,
                     'amb_temp_ctl_pt' : '0.0',
                     'amb_salt_ctl_pt' : '0.0',
                     'amb_depth_ctl_pt': '0.0',
                     'tiuniform' : None,
                     'min_melt_depth' : 0.0,
                     'phi' : None,
                     'ah' : 100.0,   
                     'kh' : 100.0,
                     'cdb' : 2.5e-3,  #bottom drag coeff    
                     'cl' : 1.775e-2, # entrainment coeff
                     'ef' : 5.0e-1,   # entrainment factor
                     'context' : '""',  #NB: don't remove this
                     'bathtype' : 13,
                     'gasp_m1' : 0.45,
                     'gasp_m2' : 2.6,
                     'gasp_m3' : 1.9,
                     'gasp_m4' : 2.3,
                     'gasp_m5' : 0.6,
                     }

    def produce_namelist_contents(self):
        
        # return a string which can be used as the contents of
        # a plume namelist file
        
        lines = []
        vItems = self.vals.items()
        vItems.sort()
        for (k,v) in vItems:
            try:
                lines.append(' %s = %s\n' % (k,_fortran_style(k,v)))
            except FortranConversionException:
                raise Exception('%s \n %s \n %s' % 
                                ('Could not write dict to Fortran namelist.',
                                 'Missing definitions of:',
                                ' '.join([k for (k,v) in
                                          self.vals.items() if (v is None)])))
        return '&plume_nml\n%s/' % ''.join(lines)


class GCConfig(object):

    def __init__(self):

        self.vals = \
        { 'parameters' :
                      { 'geothermal' :   0.0,      # geothermal heat flux 
                                                   # model%paramets%geot
                        'default_flwa' : None,     # Glen's law A to use in isothermal case 
                                                   # model%paramets%default_flwa
                        'flow_factor' : 1,         # enhancement factor for glen's A 
                                                   # model%paramets%flow_factor
                        'ice_limit' : 50.0,        # minimum thickness for running ice dynamics
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
          'Petermann shelf' :
                   {  'air_temperature' : -15,     # Temp assigned to ice in isothermal case
                      'accumulation_rate' : 0.0,   # In meters per year
                      'eustatic_sea_level' : 0.0, # Height of sea level relative to initial height
                      'check_for_steady': True,
                      'thk_steady_tol' : 1.0e-2,  #relative change in thickess,
                                                  # below which we assume the shelf is steady
                      'mean_thk_steady_tol' : 1.0e-6,#steadiness requires ave value of
                                                     # |d/dt (ln thk)| less  than this
                                          
                      },
          'options' :
            {'flow_law' : 0,     # flow_law = 2 means constant A
                                              #          = 0 means calculate A from temperature
                                              # model%options%whichflwa
             'evolution' : 3,    # evolution = 0 means pseudo-diffusion
                                           #           = 1 means ADI scheme
                                           #           = 2 means iterated diffusion
                                           #           = 3 means LANL incrementral remapping method
                                           #           = 4 upwind advection
                                           # model%options%whichevol
             'temperature' : None,        # temperature = 1 means full 3D thermal evolution
                                                    #           = 0 means set to air temperature
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
                                           #    = 1 means basal traction is constant (basal_tract_const)
                                           # = 2 means Set to (non--zero) constant where where temperature
                                           #  is at pressure melting point of ice, otherwise to zero
                                           #= 3 means function of basal water depth 
                                           # model%options%whichbtrc
             'periodic_ew' : 0,  # model%options%periodic_ew
             'periodic_ns' : 0,  # model%options%periodic_ns
             'diagnostic_run' : 0,  # = 1 makes glide stop after diagnosing velocities
                                              # model%options%diagnostic_run
             'hotstart' : None,   # are we doing a restart (1 = yes, 0 = no)
                                            # model%options%hotstart
             'basal_water' : 3,   # basal_water = 0 means calc from local basal water balance
                                            #  = 1 means compute basal water flux then find depth via cal
                                            #  = 2 means no basal water
                                            # model%options%whichbwat

             'use_plume' : None,  # use_plume = 0 means Glimmer-CISM method to calc bmlt
                                        #  = 1 means use plume model to calculate bmlt above floating ice
                                           # model%options%use_plume
             # 'ioparams' :                #i/o parameters file
             },
         
          'grid' :
              { 'sigma_builtin' : 1,   # sigma_builtin = 0 means use the default sigma levels 
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
                           'which_disp' : 0,         # dissipation
                                                    # model%options%which_disp
                           'which_ho_resid' : 0,  # method of calculating residual: 
                                                  # = 0 means use max value
                                                  # NB: only used in glam_velo_fordsiapstr
                                                 # model%options%which_ho_resid
                           'which_bmelt' : 0,   # basal melting
                                                # model%options%which_bmelt
                           'which_ho_babc' : 6,       # basal boundary condition: 
                                                   # = 5 means simple ice-shelf
                                                    # = 6 means floating ice everywhere, no traction
                                                   # = 3 means circular ice-shelf
                                                 # NB: only used when using Payne-Price diagnostic scheme
                                                  # model%options%which_ho_babc
                           'which_ho_efvs' : 0,   #  0 = calculate eff visc from strain rates
                                                  #  1 = use constant value of 0.001                                               
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
                           'which_ho_efvs' : 0,    # ho effective viscosity 
                                                   # = 0 means calculate from strain rate
                                                   # model%options%which_ho_efvs
                           },
          'picard parameters' : { 'minres' : 1.0e-5,
                                  'switchres' : 1.0e-2,
                                  'x_overrideres' : 0.0,
                                  'y_overrideres' : 1.0e-8 ,
                                  'cmax' : 3000,
                                  'cmin' : 5,
                                  'cswitch' : 100,
                                  'cvg_accel' : 1.5,
                                  'small_vel' : 0.001,
                                  'start_umc' : 3,
                                  'x_invariant' : 0,  # = 1 means variables don't change in x direction
                                                      # model%picard_params%x_invariant
                                  },
          'boundary condition params' :
              {'use_lateral_stress_bc' : True,
               'use_plastic_bnd_cond' : False,
               'tau_xy_0' : None,
               'annual_percent_var' : 0.0,
               'use_shelf_bc_1' : False,
               'use_sticky_wall' : False, # create a 'sticky spot' along wall or not
               'sticky_length' : 0,
               },
          'CF default' : { 'comment' : '',
                           'title' : None,
                           'institution' : 'NYU',
                           'references' : ''
                           },
          
          'CF input' : { 'name' : None,     #name of netcdf file containing input data fields
                         'time' : None      #which time slice to read input data from
                         },
          'CF output' : { 'variables' : ' '.join(['acab','lsurf','usurf',
                                                  'thk','thk_t','bmlt',
                                                  'kinbcmask',
                                                  'uflx_conv','vflx_conv','flx_conv',
                                                  'uvelhom',
                                                  'vvelhom',
                                                  'thkmask','topg',
                                                  'temp','efvs',
                                                  'tau_hom_xx','tau_hom_yy',
                                                  'tau_hom_xz','tau_hom_yz','tau_hom_xy']),
                          ### NB: there is a (250) character limit on line length!!!
                          # the following is a list of all possible output variables:
                                   # level lithoz x0 x1 y0 y1 acab acab_tavg age artm backstress
                                   # beta bheatflx bmlt bmlt_tavg btemp btrc bwat bwatflx calving 
                                   # diffu dusrfdtm eus flwa gl_ew gl_ns gline_flux iarea ivol
                                   # kinbcmask lat litho_temp lon lsurf relx slc soft surfvel tau_xz tau_yz
                                   # taux tauy temp thk thkmask topg ubas ubas_tavg uflx usurf
                                # uvel uvelhom vbas vbas_tavg velnormhom vflx vvel vvelhom wgrd wvel efvs
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
                      'plume_do_cross_shelf_avg' : False,
                                                          # NB it is very storage hungry
                      'plume_write_every_n' : 1,
                      'plume_output_frequency' : 0.0,
                      'plume_min_spinup_time' : 5.0,    # minimum time to spinup the plume, in days
                      'plume_max_spinup_time' : 100.0,   # maximum time to spinup, in days
                      'plume_min_subcycle_time' : 0.5,   # minimum subcycle time, in days
                      'plume_max_subcycle_time' : 40.0,  # maximum subcycle time, in days
                      'plume_steadiness_tol' : 1.0e-6,  # plume steadiness tolerance
                      'plume_speed_steadiness_tol' : 1.0e-6,  # max relative change in speed
                      'plume_imin' : None,
                      'plume_imax' : None,
                      'plume_kmin' : None,
                      'plume_kmax' : None,
                      'plume_initial_bmlt' : False,     # Apply the initial bmelt field for all time
                      'plume_const_bmlt' : False,       # Apply a uniform melt rate under floating ice
                      'plume_const_bmlt_rate' : 0.0     # At given rate in meters per year
                      }
          }
        
        
    def produce_config_file(self):

        sections = []

        vItems = self.vals.items()
        vItems.sort()
        
        for (k,v) in vItems:
            try:
                vdictItems = v.items()
                vdictItems.sort()
                section_vals = [' %s = %s' %
                                (k_inner,_fortran_style(k_inner,v))
                                for (k_inner,v) in vdictItems]
                
            except FortranConversionException, (e):
                raise Exception(
                    'Could not produce config file.  Failed to convert %s -> %s\n\n' % (k,e.param_name))
            
            sections.append('[%s]\n%s\n' % (k, '\n'.join(section_vals)))
                        
        return '\n'.join(sections)

############################################################
### classes that represent glimmer-cism-plume jobs #########
############################################################
    
class _BaseJob(object):
    '''To define a job, the typical calling sequence should be:
    j = JobClass()
    j.field1 = val1
    ...
    j.fieldN = valN

    j.assertCanStage()
    j.serialize()

    To run a job, just call:
        j.resolve()
        j.stage()
        j.run()
    '''

    def __init__(self):
        ''' A virtual class that defines the basic interface for a job.
        It has the public functions: assertCanStage, resolve, stage, serialize, run.
        '''
        self._name = None
        self._driver = 'shelf_driver'

        self.completed = False
        self.started = False
        self.timeStop = 0.0
        self.timeStart = 0.0
        self.timeStartStr = ''
        self.timeStopStr = ''
        self.error = False
        self.errorMessage = ''
        
        #any constants that should basically never change
        self.kinbcw = 2   #width of band on which ice vel is specified
        self.plume_landw = 2 #width of land cell padding for plume grid
        self.ice_zero_thk_buf = 1 # width of zero thickness on sides of ice domain
        
        self.total_side_buf = self.ice_zero_thk_buf
        self.total_side_buf_east = self.total_side_buf
        self.total_side_buf_west = self.total_side_buf

        self._pnl = PlumeNamelist()
        self._gcconfig = GCConfig()

        #self.plume and self.gc are used to directly specify parameter values that
        # should override the values already in self._pnl and self._gcconfig
        self.plume = {}
        self.gc = defaultdict({})

    def _getJobDir(self):
        jd = os.path.join(os.path.expandvars('$GC_JOBS'),self.name)
        if (not(os.path.exists(jd))):
            os.makedirs(jd)
        return os.path.abspath(jd)
    def _get_input_cmds_log(self):
        return os.path.join(self.jobDir,'input_cmds.log')
    def _getserialfile(self):
        return os.path.join(self.jobDir,'%s.gcpl' % self.name)
    def _getUnderlyingJob(self):
        raise Exception('must override this method')
    def _setname(self,n):
        self._name = n
    def _getname(self):
        if (self._name is None):
            raise Exception("self.name is not defined")
        return self._name
    def _assertName(self):
        if (self._name is None):
            raise Exception("Name not defined yet")
    def _get_inputfile(self):
        self._assertName()
        return os.path.join(self.jobDir,'%s.in.nc' % self.name)
    def _get_outputfile(self):
        self._assertName()
        return os.path.join(self.jobDir,"%s.out.nc" % self.name)
    def _get_gc_config_file(self):
        self._assertName()
        return os.path.join(self.jobDir,"%s.config" % self.name)
    def _get_plume_nl_file(self):
        self._assertName()
        return os.path.join(self.jobDir,"%s.nl" % self.name)
    def _get_plume_outputfile(self):
        self._assertName()
        return os.path.join(self.jobDir,"plume.%s.out.nc" % self.name)
    def _get_use_plume(self):
        if ('options' in self.gc):
            if ('use_plume' in self.gc['options']):
                return self.gc['options']['use_plume']
        return self._gcconfig.vals['options']['use_plume']
                
    name = property(fset=_setname,fget=_getname)
    jobDir = property(fget=_getJobDir)
    underlyingJob = property(fget=_getUnderlyingJob)
    serialFile = property(fget=_getserialfile)
    outputfile = property(fget=_get_outputfile)
    inputfile = property(fget=_get_inputfile)
    gc_config_file = property(fget=_get_gc_config_file)
    input_cmds_log = property(fget=_get_input_cmds_log)
    plume_output_file = property(fget=_get_plume_outputfile)
    plume_nl_file = property(fget=_get_plume_nl_file)
    use_plume = property(fget=_get_use_plume)
    
    def assertCanStage(self):
        '''Check that all the fields of this job have
        been assigned a value.'''
        notFound = []
        for (k,v) in self.__dict__.items():
            if (k.startswith('_')):
                continue
            if (v is None):
                notFound.append(k)
        if len(notFound):
            raise Exception("Failed to find values for %s" % ','.join(notFound))

    def resolve(self,gc={},plume={}):
        raise Exception("must override")
    
    def stage(self,genInput=False):

        if (genInput):
            raise Exception("Don't know how to genInput in _BaseJob")
        
        #update the config dictionaries
        for (section_name,section_data) in self.gc.items():
            if (not (section_name in self._gcconfig.vals)):
                raise Exception("Non-standard section name %s" % section_name)
            else:
                template_section_data = self._gcconfig.vals[section_name]
                for (key,val) in section_data.items():
                    if (not (key in template_section_data)):
                        raise Exception("Non-recognized key %s in section %s" %
                                        (key, section_name))
                    
                self._gcconfig.vals[section_name].update(section_data)

        for (key,val) in self.plume.items():
            if (not (key in self._pnl.vals)):
                raise Exception("Non-recognized key %s in plume config" % key)
        self._pnl.vals.update(self.plume)

        if (self.use_plume):
            # then write out the plume namelist file using this job's name
            f = open(self.plume_nl_file,'w')
            try:
                f.write(self._pnl.produce_namelist_contents())
            finally:
                f.close()

        # and write out the gc config file
        g = open(self.gc_config_file,'w')
        try:
            g.write(self._gcconfig.produce_config_file())
        finally:
            g.close()
            

    def serialize(self):

        try:
            f = open(self.serialFile, 'w')
        except:
            raise Exception("Couldn't create file %s to store serialized job" % self.serialFile)
        try:
            pickle.dump(self, f)
        finally:
            f.close()

    def run(self,overwrite=False):
        if (not(overwrite) and \
            os.path.exists(self.outputfile)):
            raise Exception("Outputfile %s already exists" %
                            self.outputfile)

        if (overwrite):
            os.chmod(self.jobDir,0770)
        
        pwd = os.path.abspath(os.curdir)
        os.chdir(self.jobDir)
        cmd = [self._driver, self.gc_config_file]
        try:
            _check_calls([cmd], self.input_cmds_log)
            #os.chmod(self.jobDir,0550)
            
        finally:
            os.chdir(pwd)

    
class _GenInputJob(_BaseJob):
    ''' To define a _GenInputJob:
    j = _GenInputJob()
    j.name = 'whatever'
    j.jobDir = 'whereever'
    
    #can use
    j.plume['dt1'] = 10.0
    j.plume['salttop'] = 34.5
    ...
    
    #check if staging will work
    j.assertCanStage()

    # write out the job to a file in the jobDir
    j.serialize() 

    Note, that when j.run() is called, it will invoke j.resolve(),
    in which values in j.plume and j.gc override values provided  using 'shortcut'
    properties like j.plume_const_bmlt

    '''
    def __init__(self):
        _BaseJob.__init__(self)

        #fields with default values
        #self.ifpos = 5
#        self.ifpos = self.plume_landw + self.n - 4
#        self.ifpos = None
        self.rhoi = 910.0
        self.rhoo = 1028.0
        self.kx = 0.0
        self.chan_amp = 0.0
        self.chan_init_length = 5000.0

        
        #fields that must be defined 
        self.m = None
        self.n = None
        self.nlevel = None
        self.hx = None
        self.hy = None
        self.tstart = None
        self.tend = None
        self.otopg = None
        self.upthk = None
        self.upvel = None
        self.plume_dt = None
        self.ice_dt = None
        self.default_flwa = None
        self.uniform_acab = None
        self.randthk = None

    def _genInputCmds(self):
        #generate the netcdf input command, assuming that this 
        #job has already been resolve'ed
        raise Exception("must override")
        

    def _getifpos(self):
        if (self.n is None):
            raise Exception('n is not defined yet')
        else:
            return self.n - 4
    ifpos = property(fget=_getifpos)
    
    def _getglpos(self):
        if (self.n is None):
            raise Exception('n is not defined yet')
        else:
            return 2
    glpos = property(fget=_getglpos)
        
    def stage(self,genInput=True):
        _BaseJob.stage(self,genInput=False)
        if (genInput):
            #figure out what command is needed to generate the netcdf input file
            cmds = self._genInputCmds()
            # and then do it
            _check_calls(cmds, self.input_cmds_log)
        
    def resolve(self,gc_override={},plume_override={}):

        self.gc['CF output'].update({'name' : self.outputfile})
        self.gc['CF input'].update({'name' : self.inputfile,
                                    'time' : 1 })
        self.gc['CF default'].update({'title' : self.name })
        
        self.gc['plume'].update( {'plume_nl_file' : self.plume_nl_file,
                                  'plume_output_file' : self.plume_output_file,
                                  'plume_output_prefix' : self.name } )
        
        self.plume['hx'] = self.hx
        self.plume['hy'] = self.hy
        self.plume['gldep'] = self.upthk
        self.plume['m_grid'] = self.m + (2*self.plume_landw)
        self.plume['n_grid'] = self.n + (1*self.plume_landw)
        self.plume['dt1'] = self.plume_dt

        self.gc['options'].update({    'temperature' : 0, #isothermal
                                       'use_plume' : self.use_plume,
                                       'hotstart' : 0,
                                       })
        self.gc['parameters'].update({ 'default_flwa' : self.default_flwa,
                                                  })
        self.gc['CF output'].update({'start' : self.tstart,
                                     'stop' : self.tend,
                                     })

        self.gc['time'].update({ 'dt' : self.ice_dt,
                                            'tend' : self.tend,
                                            'tstart' : self.tstart,
                                            })
        self.gc['grid'].update({'dew' : self.hx,
                                'dns' : self.hy,
                                'ewn' : self.m,
                                'nsn' : self.n,
                                'upn' : self.nlevel
                                
                                })

        self.gc['plume'].update(
            {'plume_imax' : self.m + 2*(self.plume_landw) - self.total_side_buf_east,
             'plume_imin' : 1+                              self.total_side_buf_west,
             #'plume_kmin' : self.ifpos,
             #'plume_kmax' : self.n + self.plume_landw,
             'plume_kmin' : 1,
             'plume_kmax' : self.ifpos + self.plume_landw,
             })


        self.plume.update(plume_override)

        for k in gc_override.keys():
            self.gc[k].update(gc_override[k])

    
class LinearShelfJob(_GenInputJob):

    def __init__(self):
        _GenInputJob.__init__(self)
        self.ifthk = None
        self.inflow_a = None
        self.noslip = False
        
    def _genInputCmds(self):
        cmd = ['nc_gen_input']
        cmd.extend(['ls', self.inputfile,
                    int(self.m), int(self.n), int(self.nlevel),
                    float(self.hx), float(self.hy),
                    float(self.upthk), float(self.upvel),
                    float(self.inflow_a), int(self.ifpos), float(self.ifthk),
                    float(self.otopg), int(self.kinbcw), self.noslip])
        cmd = [_fortran_style('nc_gen_input', c) for c in cmd]
        return [cmd]

    def _getUnderlyingJob(self):
        return self
    underlyingJob = property(fget=_getUnderlyingJob)
    
class SandersonShelfJob(LinearShelfJob):

    def __init__(self):
        LinearShelfJob.__init__(self)

        self.ifthk = 0.0
        self.annual_percent_var = 0.0
        self.tauxy0 = None


    def resolve(self,gc_override={},plume_override={}):
        
        #shelf_length = (self.n-self.kinbcw-(self.ifpos-1))*self.hy
        shelf_length = (self.ifpos-(self.kinbcw-1))*self.hy
        widthY = (self.m- 2*1)*self.hx/2.0
        g = 9.81
        slope = 2*abs(self.tauxy0)/(self.rhoi*g*widthY*(1-self.rhoi/self.rhoo))
        self.ifthk = self.upthk - slope*shelf_length

        self.gc.update(  {'boundary condition params' : {'tau_xy_0' : self.tauxy0,
                                                         'annual_percent_var' : self.annual_percent_var,
                                                         },
                          } )

        LinearShelfJob.resolve(self,gc_override,plume_override)
        
class SteadyShelfJob(_GenInputJob):
    
    def __init__(self):
        _GenInputJob.__init__(self)

    def _genInputCmds(self):
        cmd = ['nc_gen_input']
        cmd.extend(['ss', self.inputfile,
                    int(self.m),int(self.n), int(self.nlevel),float(self.hx),float(self.hy),
                    float(self.upthk),float(self.upvel),int(self.ifpos),
                    float(self.otopg), int(self.kinbcw), float(self.uniform_acab), float(self.rhoi),
                    float(self.rhoo), float(self.default_flwa),float(self.randthk),
                    float(self.kx), float(self.chan_amp),float(self.chan_init_length)])
        cmd = [_fortran_style('nc_gen_input', c) for c in cmd]
        return [cmd]


class RestartIceJob(_BaseJob):
    
    def __init__(self,initJobName,initJobDir=None,newName=None):
        _BaseJob.__init__(self)

        #first locate the initial job directory
        if (initJobDir is None):
            _initJobDir = os.path.join(os.path.expandvars('$GC_JOBS'),
                                       initJobName)
        else:
            _initJobDir = initJobDir
            
        initJobFile = os.path.join(_initJobDir,
                                   '%s.gcpl' % initJobName)
        f = open(initJobFile,'r')
        try:
            self._initJob = pickle.load(f)
        finally:
            f.close()
        
        _initJobName = self._initJob.name

        self._initJobInputNcFile = os.path.join(_initJobDir,
                                               os.path.basename(self._initJob.outputfile))

        # figure out new name for this job
        if (newName is None):
            if (len(_initJobName.split('_restart_')) > 1):
                try:
                    newindex = int(_initJobName.split('_restart_')[1])+1
                    self.name = "%s_restart_%s" % (_initJobName.split('_restart_')[0],newindex)
                except ValueError:
                    self.name = _initJobName
            else:
                self.name = "%s_restart_%s" % (_initJobName, 1)
        else:
            self.name = newName

        #parse the output of the old job to figure out the last time and index
        p = subprocess.Popen(['ncdump','-v', 'time', self._initJobInputNcFile],stdout=subprocess.PIPE)
        times = p.stdout.read()

        self._inputJobLastTimeIndex = int(times.split('(')[1].split('currently')[0].strip())
        self.tstart = float(times.split('time =')[-1].split()[-3])
        self.tend = self._initJob._gcconfig.vals['time']['tend'] # default value
        
        self.gc = defaultdict({})
        self.plume = {}

        self.doPlumeRestart = True
        self.newJob = copy.copy(self._initJob)
        self.newJob.name = self.name
        
        self._plume_restart_file = os.path.join(self.jobDir,
                                                'last_time.%s' %
                                                os.path.basename(self._initJob.plume_output_file))
        self.useMaxRunTimeLimit = False
        self.maxRunTime = 0.0

        self.restartIceIndex = -1
        
    def _getUnderlyingJob(self):
        return self.newJob.underlyingJob
    underlyingJob = property(fget=_getUnderlyingJob)

    def _get_m(self):
        return self.newJob.m
    m = property(fget=_get_m)

    def _get_n(self):
        return self.newJob.n
    n = property(fget=_get_n)

    def resolve(self,gc_override={},plume_override={}):
        #resolve the contained job
        self.gc['time']['tstart'] = self.tstart
        self.gc['CF output']['start'] = self.tstart

        if (self.useMaxRunTimeLimit):
            self.gc['CF output']['stop'] = self.tstart + self.maxRunTimeLimit
            self.gc['time']['tend'] = self.tstart + self.maxRunTimeLimit
        else:
            self.gc['CF output']['stop'] = self.tend
            self.gc['time']['tend'] =     self.tend

        if (self.restartIceIndex > -1):
            if (self.restartIceIndex > self._inputJobLastTimeIndex):
                raise Exception("Ice index %s exceed max available index of %s" % (self.restartIceIndex,
                                                                                   self._inputJobLastTimeIndex))
        else:
            self.restartIceIndex = self._inputJobLastTimeIndex
            
        if ('doPlumeRestart' in self.__dict__):
            if (self.doPlumeRestart):
                self.plume['restart'] = True
                try:
                    self.plume['restart_data_filename'] = '"%s"' % self._plume_restart_file
                except:
                    pass
            else:
                self.plume['restart'] = False
        else:
            self.plume['restart'] = False
            
        for k in gc_override.keys():
            self.gc[k].update(gc_override[k])

        self.plume.update(plume_override)
        
        self.newJob.name = self.name
        self.underlyingJob.name = self.name

        self.newJob.resolve(self.gc,self.plume)

    def assertCanStage(self):
        self.newJob.assertCanStage()
        _BaseJob.assertCanStage(self)

    def stage(self,genInput=True):

        self.assertCanStage()
        self.newJob.stage(genInput=False)
        if (genInput):
            cmds = self._genInputCmds()
            _check_calls(cmds, self.input_cmds_log)
    
    def _genInputCmds(self):
        cmd =['nc_regrid']
        cmd.extend([self._initJobInputNcFile,
                    self.inputfile,
                    self.restartIceIndex,
                    self.m,self.n,-1,
                    0,0,0,0, 
#                    0,4,1,1,
#                    self.newJob.kinbcw, 4,0,0,
                    4,2,1,1,
                    0,self.newJob.kinbcw,0,0,
                    0.0,           0.0])             
        cmd1 = [_fortran_style('nc_regrid', c) for c in cmd]
        
        cmds = [cmd1]

        if (self.doPlumeRestart):
            cmd2 = ['nc_last_time_slice.py', self._initJob.plume_output_file,
                    '-o%s' % self._plume_restart_file.strip()]
            cmds.append(cmd2)
        return cmds

class IntroGLPerturbJob(RestartIceJob):

    def __init__(self, initJobName, initJobDir=None, newJobName=None):
        RestartIceJob.__init__(self,initJobName,initJobDir,newJobName)
        
        self.k = None
        self.amp = None
        self.ramp_len = None
        self.inflow_a = None
        
    def _genInputCmds(self):

        cmd = ['nc_regrid']
        cmd.extend([self._initJobInputNcFile,
                    self.inputfile,
                    self.restartIceIndex,
                    self.m,self.n,-1,
                    0,0,0,0, 
                    4,2,1,1,
                    0,2,0,0,
                    self.inflow_a,
                    self.k, self.amp, self.ramp_len])
        
        cmd = [_fortran_style('nc_perturb_gl',c) for c in cmd]

        return [cmd]

class ListPerturbJob(RestartIceJob):

    def __init__(self, initJobName, initJobDir=None, newJobName=None):
        RestartIceJob.__init__(self,initJobName,initJobDir,newJobName)

        #perturb_list should have the form of a list of quadruples
        self.perturb_list = []
        self.inflow_a = None
        self.vvelhom_new_val = None
        
    def _genInputCmds(self):

        cmd = ['nc_regrid']
        cmd.extend([self._initJobInputNcFile,
                    self.inputfile,
                    self.restartIceIndex,
                    self.m,self.n,-1,
                    0,0,0,0, 
                    4,2,1,1,
                    0,2,0,0,
                    self.inflow_a,
                    self.vvelhom_new_val*1.0,
                    ])
        for (k,amp,phase,len) in  self.perturb_list:
            cmd.extend([k,amp,phase,len])
        
        cmd = [_fortran_style('nc_perturb_gl',c) for c in cmd]

        cmds = [cmd]
        
        if (self.doPlumeRestart):
            cmd2 = ['nc_last_time_slice.py', self._initJob.plume_output_file,
                    '-o%s' % self._plume_restart_file.strip()]
            cmds.append(cmd2)

        return cmds
    
class RegridListPerturbJob(ListPerturbJob):

    def __init__(self, initJobName, initJobDir=None, newName=None):
        ListPerturbJob.__init__(self,initJobName,initJobDir, newName)
        self.new_m = None
        self.new_n = None
        self.new_level = None

    def resolve(self,gc_override={},plume_override={}):
        gc_override.update( {'grid' : {'upn' : self.new_level,
                                       }} )
        ListPerturbJob.resolve(self,gc_override,plume_override)
        
    def _genInputCmds(self):

        cmd =['nc_regrid']
        cmd.extend([self._initJobInputNcFile,
                    self.inputfile,
                    self.restartIceIndex,
                    self.new_m,self.new_n,self.new_level,
                    0,0,0,0,
                    4,2,1,1,
                    0,2,0,0,
                    0.0,0.0])

        for (k,amp,phase,len) in  self.perturb_list:
            cmd.extend([k,amp,phase,len])
        
        cmd = [_fortran_style('nc_perturb_gl_regrid',c) for c in cmd]
        cmds = [cmd]
        
        if (self.doPlumeRestart):
            cmd2 = ['nc_last_time_slice.py', self._initJob.plume_output_file,
                    '-o%s' % self._plume_restart_file.strip()]
            cmds.append(cmd2)

        return cmds
                 
class RegridIceJob(RestartIceJob):

    def __init__(self,initJobFile,newName=None):
        RestartIceJob.__init__(self,initJobFile,initJobDir=None,newName=newName)
        self.new_m = None
        self.new_n = None
        self.new_level = None
        self.restartIceIndex = None
        
    def _genInputCmds(self):

        cmd =['nc_regrid']
        cmd.extend([self._initJobInputNcFile,
                    self.inputfile,
                    self.restartIceIndex,
                    self.new_m,self.new_n,self.new_level,
                    0,0,0,0,
                    #2,4,1,1,  #n,s,e,w thickness buffers
                    #4,0,1,1,
                    4,2,1,1,
                    0,2,0,0,
                    0.0,0.0])             
        cmd = [_fortran_style('nc_regrid', c) for c in cmd]
        return [cmd]
    
class FixedBasalMeltJob(RestartIceJob):

    def __init__(self,initJobName,newName=None):
        RestartIceJob.__init__(self,initJobName,newName=newName)

        self.underlyingJob.gc['options']['use_plume'] = 1
        self.plume_dt = 1.0
        self.plume.update( {'entr_time_const' : 600.0,
                            'phi' : 0.0,
                            'salttop' : 0.0,
                            'saltbot' : 0.0,
                            'tempbot' : 0.0,
                            'temptop' : 0.0,
                            'tiuniform' : 0.0,
                            'restart' : True,
                            'restart_data_filename' : '"%s"' % self.plume_input,
                            })
                      
        self.gc.update( {'plume' : { 'plume_const_bmlt' : False,
                                     'plume_initial_bmlt' : True,
                                     },
                         } )

        self.plume_m = None
        self.plume_n = None
        self.plume_hx = None
        self.plume_hy = None
        
        self._plume_restart_file = '"plume_input.out.nc"'
        
    def resolve(self,gc_override={},plume_override={}):
        plume_override.update( {'restart' : True,
                                'restart_data_filename' : '"%s"' % self.plume_input})
        RestartIceJob.resolve(self,gc_override,plume_override)
        

    def _getPlumeInputFile(self):
        return os.path.join(self.jobDir,
                            'plume_input.out.nc')
    plume_input = property(fget=_getPlumeInputFile)

#    def resolve(self,gc_override={},plume_override={}):
#        self.gc.update(gc_override)
#        self.plume.update(plume_override)
#        self.plume.update( {'restart_data_filename' : '"%s"' % self.plume_input,
#                            'restart' : True,
#                            } )

#        RestartIceJob.resolve(self,self.gc,self.plume)
        
    def _genPlumeInputCmd(self):
        cmd1 = ['nc_gen_plume']
        #print(self.newJob.gc)
        cmd1.extend(['fb', self.plume_input,
                    int(self.plume_m),
                    int(self.plume_n),
                    float(self.plume_hx),
                    float(self.plume_hy)])
        cmd1 = [_fortran_style('nc_gen_plume', c) for c in cmd1]
        
        #cmd2 = LinearShelfJob._genInputCmds(self)
        #cmd2 = [[_fortran_style('nc_gen_plume', c) for c in cmd] for cmd in cmd2]

        #cmds = [cmd1]
        #cmds.extend( cmd2 )
        #return cmds
        return cmd1

    def stage(self,genInput=True):
        RestartIceJob.stage(self,genInput)
        cmd = self._genPlumeInputCmd()
        _check_calls([cmd], self.input_cmds_log)


def _check_calls(cmds,logfile=None):

    if (logfile is not None):
        f = open(logfile,'a')
        try:
            f.write('\n')
            f.write('\n\n'.join([' '.join(cmd) for cmd in cmds]))
            f.write('\n\n')
            
        finally:
            f.close()

    for cmd in cmds:

        retcode = subprocess.call(cmd)
        if (retcode != 0):
            raise Exception("Error running:\n       %s" % ' '.join(cmd))
    
    
            
class FromFilesJob(_BaseJob):

    def __init__(self,jobDir):
        _BaseJob.__init__(self)

        if (not(os.path.exists(jobDir))):
            jobDir2 = os.path.join(os.path.expandvars('$GC_JOBS'),
                                   jobDir)
            if (not(os.path.exists(jobDir2))):
                raise Exception("Can not find %s" % jobDir)
            else:
                jobDir = jobDir2

        config_files = [x for x in os.listdir(jobDir) if
                        (x.endswith('.config') and not('captured' in x))]

        if (len(config_files) != 1):
            raise Exception("Did not find a unique .config file")


        self.name = config_files[0].split('.config')[0]

        # parse the .config file
        f_config = open(os.path.join(jobDir,config_files[0]),'r').read()

        sections = re.findall('\\[([a-z,A-Z,0-9,_,\s]*)\\]([0-9,.,_,+,\-,\s,a-z,A-Z,=,/]*)',f_config)

        d = {}
        
        for section in sections:
            heading = section[0]
            d[heading] = {}
            
            body = section[1]
            kvpairs = re.findall('([^\n]+=[^\n]+)',body)
            for kvp in kvpairs:
                k = kvp.split('=')[0].strip()
                v = kvp.split('=')[1].strip()
                try:
                    v = int(v)
                except:
                    try:
                        v = real(v)
                    except:
                        pass
                d[heading][k] = v

        for k in d.keys():
            self._gcconfig.vals[k].update(d[k])

        if (self.use_plume):
            # now parse the nl file
            nl_fname =os.path.join(self.jobDir,
                                   config_files[0].split('.config')[0]+'.nl')
            if (not os.path.exists(nl_fname)):
                raise Exception("Could not find expected namelist file %s" % nl_fname)

            nl_file = open(nl_fname,'r').read()

            kvpairs = re.findall('([^\n]+)=([^\n]+)',nl_file)
            for (k,v) in kvpairs:
                v = v.strip()
                try:
                    v = int(v)
                except:
                    try:
                        v = real(v)
                    except:
                        pass

                self._pnl.vals[k.strip()] = v

    def _genInputCmds(self):
        return []

    def _get_m(self):
        return int(self._gcconfig.vals['grid']['ewn'])
    m = property(fget=_get_m)

    def _get_n(self):
        return int(self._gcconfig.vals['grid']['nsn'])
    n = property(fget=_get_n)

    def _getUnderlyingJob(self):
        return self
    underlyingJob = property(fget=_getUnderlyingJob)
    
    def resolve(self,gc_override={},plume_override={}):

        self.gc['CF output'].update({'name' : self.outputfile})
        self.gc['CF input'].update({'name' : self.inputfile,
                                    'time' : 1 })
        self.gc['CF default'].update({'title' : self.name })
        
        self.gc['plume'].update( {'plume_nl_file' : self.plume_nl_file,
                                  'plume_output_file' : self.plume_output_file,
                                  'plume_output_prefix' : self.name } )
        
        for k in gc_override.keys():
            self.gc[k].update(gc_override[k])
        self.plume.update(plume_override)
        
    def stage(self, genInput=False):
        if (genInput):
            raise Exception("can not generate input")

        _BaseJob.stage(self,genInput=False)


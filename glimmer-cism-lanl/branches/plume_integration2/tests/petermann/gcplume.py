import sys
import os
import subprocess
import pickle

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
                     'tangle' : False, 
                     'negfrz' : False,
                     'use_min_plume_thickness' : True,
                     'tottim'  : 0.0,
                     'labtim'  : 1.0,
                     'snottim' : 24.0*3600.0,  #once per day
                     'lnottim' : 5.0*86400.0, #once per five days
                     'dt1'     : None,
                     'm_grid' : None,
                     'n_grid' : None,          
                     'hx' : None,      
                     'hy' : None,
                     'gldep' : None,
                     'ifdep' : 0.0,
                     'wcdep' : 1000.0,
                     'plume_min_thickness' : 10.0,
                     'bsmoothit' : 0,   #iterations of smoothing
                     'salttop' : None,
                     'saltbot' : None,
                     'temptop' : None,
                     'tempbot' : None,
                     'phi' : None,
                     'ah' : 100.0,   
                     'kh' : 100.0,
                     'cdb' : 2.5e-3,  #bottom drag coeff    
                     'cl' : 1.775e-2, # entrainment coeff
                     'ef' : 5.0e-1,   # entrainment factor
                     'context' : '""',  #NB: don't remove this
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
        return '&plume_nml\n %s /' % ' ,'.join(lines)


class GCConfig(object):

    def __init__(self):

        self.vals = { 'parameters' : { 'geothermal' :   0.0,      # geothermal heat flux 
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
                      'Petermann shelf' : {  'air_temperature' : -5,     # Temp assigned to ice in isothermal case
                                             'accumulation_rate' : 0.0,   # In meters per year
                                             'eustatic_sea_level' : 0.0, # Height of sea level relative to initial height
                                             'check_for_steady': True,
                                             'thk_steady_tol' : 1.0e-5,  #relative change in thickess,
                                                                         # below which we assume the shelf is steady
                                             },
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
                                              },
                      'boundary condition params' : {'use_lateral_stress_bc' : True,
                                                     'use_plastic_bnd_cond' : False,
                                                     'tau_xy_0' : 50.0e3,
                                                     'use_shelf_bc_1' : False,
                                                     'use_sticky_wall' : False, # create a 'sticky spot' along wall or not
                                                     'sticky_length' : 0,
                                                     'x_invariant' : 0,  # = 1 means variables don't change in x direction
                                                                         # model%picard_params%x_invariant

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
                                                              #uvelhom_srf',
                                                              #'vvelhom_srf',
                                                              'thkmask','topg',
                                                              'kinbcmask',
                                                              'temp',
							      'efvs',
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
                                                                      # NB it is very storage hungry

                                  'plume_min_spinup_time' : 5.0,    # minimum time to spinup the plume, in days
                                  'plume_max_spinup_time' : 100.0,   # maximum time to spinup, in days
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

class _HasJobDir(object):
    
    def __init__(self):
        self._jobDir = None

    def _setJobDir(self,jd):
        jd = os.path.expandvars(jd)
        if (not(os.path.exists(jd))):
            raise Exception("Directory %s does not exist" % jd)
        self._jobDir = os.path.abspath(jd)
    def _getJobDir(self):
        if (self._jobDir is None):
            raise Exception("Jobdir was not specified")
        return self._jobDir
    jobDir = property(fset=_setJobDir,fget=_getJobDir)

class _IO(object):
    ''' A mixin class that adds properties:
    inputfile, outputfile,gc_config_file,plume_nl_file,plume_output_file
    '''
    def _assertName(self):
        if (self.name is None):
            raise Exception("Name not defined yet")

    def get_inputfile(self):
        self._assertName()
        return os.path.join(self.jobDir,
                            '%s.in.nc' % self.name)
    inputfile = property(fget=get_inputfile) 

    def get_outputfile(self):
        self._assertName()
        return os.path.join(self.jobDir,
                            "%s.out.nc" % self.name)
    outputfile = property(fget=get_outputfile)
    
    def get_gc_config_file(self):
        self._assertName()
        return os.path.join(self.jobDir,
                            "%s.config" % self.name)
    gc_config_file = property(fget=get_gc_config_file)

    def get_plume_nl_file(self):
        self._assertName()
        return os.path.join(self.jobDir,
                            "%s.nl" % self.name)
    plume_nl_file = property(fget=get_plume_nl_file)

    def get_plume_outputfile(self):
        self._assertName()
        return os.path.join(self.jobDir,
                            "plume.%s.out.nc" % self.name)
    plume_output_file = property(fget=get_plume_outputfile)


class _BaseJob(_HasJobDir):
    '''To define a job, the typical calling sequence should be:
    j = JobClass()
    j.field1 = val1
    ...
    j.fieldN = valN

    j.assertCanStage()
    j.serialize()

    To run a job, just call:
        j.stage()
        j.run()
    '''

    def __init__(self):
        ''' A virtual class that defines the basic interface for a job.
        It has the public functions: assertCanStage, stage, serialize, run.
        '''
        _HasJobDir.__init__(self)

        self._name = None
        
        self.completed = False
        self.started = False
        self.timeStop = 0.0
        self.timeStart = 0.0
        self.timeStartStr = ''
        self.timeStopStr = ''
        self.error = False
        self.errorMessage = ''
        
        #any constants that should basically never change
        self.kinbcw = 2
        self.plume_landw = 2


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

    def stage(self,genInput=True):
        raise Exception("must override")

    def run(self):
        raise Exception("must override")

    def _getserialfile(self):
        return os.path.join(self.jobDir,
                            '%s.gcpl' % self.name)

    serialFile = property(fget=_getserialfile)

    def serialize(self):
        if (self.name is None):
            raise Exception("No name was assigned to this job")


        if (not( os.path.lexists(self.jobDir))):
            os.mkdir(self.jobDir)
        try:
            f = open(self.serialFile, 'w')
        except:
            raise Exception("Couldn't create file %s to store serialized job" % self.serialFile)
        try:
            pickle.dump(self, f)
        finally:
            f.close()

    def _setname(self,n):
        self._name = n
    def _getname(self):
        if (self._name is None):
            raise Exception("self.name is not defined")
        return self._name
    name = property(fset=_setname,fget=_getname)

    
class _AtomicJob(_BaseJob,_IO):
    ''' To define an _AtomicJob:
    j = _AtomicJob()
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

    Note, that when j.run() is called, it will invoke j._resolveJob(),
    in which values in j.plume and j.gc override values provided  using 'shortcut'
    properties like j.plume_const_bmlt

    '''
    def __init__(self):
        _BaseJob.__init__(self)

        self._driver = 'shelf_driver'
        
        self._pnl = PlumeNamelist()
        self._gcconfig = GCConfig()
        self.plume = {}
        self.gc = {}

    def _genInputCmd(self):
        #generate the netcdf input command, assuming that this 
        #job has already been _resolveJob'ed
        raise Exception("must override")
        
    def run(self):
        pwd = os.path.abspath(os.curdir)
        os.chdir(self.jobDir)
        cmd = [self._driver, self.gc_config_file]
        _check_call(cmd)
        os.chdir(pwd)

    def stage(self,genInput=True):

        #first give the job a chance to populate 
        # self._pnl and self._gcconfig with values
        # derived from the job definition
        self._resolveJob()

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
            
        if (genInput):
            #figure out what command is needed to generate the netcdf input file
            cmd = self._genInputCmd()
            # and then do it
            _check_call(cmd)

class _GenInputJob(_AtomicJob):

    def __init__(self):
        _AtomicJob.__init__(self)

        #fields with default values
        self.ifpos = 5
        self.rhoi = 910.0
        self.rhoo = 1028.0

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
        self.use_plume = None

    def _resolveJob(self):
        self._pnl.vals['hx'] = self.hx
        self._pnl.vals['hy'] = self.hy
        self._pnl.vals['gldep'] = self.upthk
        self._pnl.vals['m_grid'] = self.m + (2*self.plume_landw)
        self._pnl.vals['n_grid'] = self.n + (1*self.plume_landw)
        self._pnl.vals['dt1'] = self.plume_dt

        self._gcconfig.vals['options'].update({    'temperature' : 0, #isothermal
                                                   'use_plume' : self.use_plume,
                                                   'hotstart' : 0,
                                                   })
        self._gcconfig.vals['parameters'].update({ 'default_flwa' : self.default_flwa,
                                                  })
        self._gcconfig.vals['CF output'].update({'frequency' : self.ice_dt,
                                                 'name' : self.outputfile,
                                                 'start' : self.tstart,
                                                 'stop' : self.tend,
                                                })
        self._gcconfig.vals['CF input'].update({'name' : self.inputfile,
                                                'time' : 1,
                                               })
        self._gcconfig.vals['CF default'].update({'title' : self.name,
                                                 })
        self._gcconfig.vals['time'].update({ 'dt' : self.ice_dt,
                                            'tend' : self.tend,
                                            'tstart' : self.tstart,
                                            })
        self._gcconfig.vals['grid'].update({'dew' : self.hx,
                                           'dns' : self.hy,
                                           'ewn' : self.m,
                                           'nsn' : self.n,
                                           'upn' : self.nlevel
                                           })
        self._gcconfig.vals['plume'].update({'plume_imax' : self.m + 2*self.plume_landw,
                                             'plume_imin' : 1,
                                             'plume_kmin' : self.ifpos,
                                             'plume_kmax' : self.n + 1*self.plume_landw,
                                             'plume_nl_file' : self.plume_nl_file,
                                             'plume_output_file' : self.plume_output_file,
                                             'plume_output_prefix' : self.name,
                                             })


class LinearShelfJob(_GenInputJob):

    def __init__(self):
        _GenInputJob.__init__(self)
        self.ifthk = None

    def _genInputCmd(self):
        cmd = ['nc_gen_input']
        cmd.extend(['ls', self.inputfile,
                    int(self.m), int(self.n), int(self.nlevel), float(self.hx), float(self.hy),
                    float(self.upthk), float(self.upvel), int(self.ifpos), float(self.ifthk),
                    float(self.otopg), int(self.kinbcw)])
        cmd = [_fortran_style('nc_gen_input', c) for c in cmd]
        return cmd

class SteadyShelfJob(_GenInputJob):
    
    def __init__(self):
        _GenInputJob.__init__(self)

    def _genInputCmd(self):
        cmd = ['nc_gen_input']
        cmd.extend(['ss', self.inputfile,
                    int(self.m),int(self.n), int(self.nlevel),float(self.hx),float(self.hy),
                    float(self.upthk),float(self.upvel),int(self.ifpos),
                    float(self.otopg), int(self.kinbcw), float(self.uniform_acab), float(self.rhoi),
                    float(self.rhoo), float(self.default_flwa),float( self.randthk)])
        cmd = [_fortran_style('nc_gen_input', c) for c in cmd]
        return cmd


ssj = SteadyShelfJob()
ssj.name = 'test'
jd = os.path.join(os.path.expandvars('$GC_JOBS'),
                          'test')
if (not(os.path.lexists(jd))):
    os.mkdir(jd)
ssj.jobDir = jd

ssj.m = 20
ssj.n = 20
ssj.upvel = -1000.0
ssj.upthk = 600.0
ssj.use_plume = 1
ssj.plume_dt = 60.0
ssj.ice_dt = 0.1
ssj.default_flwa = 1.0e-16
ssj.uniform_acab = 0.0
ssj.nlevel = 3
ssj.tend = 100.0
ssj.tstart = 0.0
ssj.hx = 1000.0
ssj.hy = 1000.0
ssj.otopg = -1200.0
ssj.randthk = 0.0
ssj.plume.update({'saltbot' : 34.5,
                'salttop' : 34.5,
                'temptop' : 0.0,
                'tempbot' : 0.0,
                  'phi'   : 0.0,
                  })
ssj.gc.update({'options' : {'flow_law' : 0,
                            },
               })
            

class RegridJob(_BaseJob,_IO):
    
    def __init__(self,initJobFile):
        _BaseJob.__init__(self)

        self.initJobFile = initJobFile


        f = open(initJobFile,'r')
        try:
            self.initJob = pickle.load(f)
        finally:
            f.close()

        _initJobDir = os.path.dirname(os.path.abspath(initJobFile))
        _initJobName = self.initJob.name

        self.inputNcFile = os.path.join(_initJobDir,
                                        self.initJob.outputfile)


        # figure out new name for this job
        if (len(_initJobName.split('_restart_')) > 1):
            newindex = int(_initJobName.split('_restart_')[1])+1
            self.name = "%s_restart_%s" % (_initJobName.split('_restart_')[0],newindex)
        else:
            self.name = "%s_restart_%s" % (_initJobName, 1)


        #parse the output of the old job to figure out the last time and index
        p = subprocess.Popen(['ncdump','-v', 'time', self.inputNcFile],stdout=subprocess.PIPE)
        times = p.stdout.read()
        self._inputNcTimeIndex = int(times.split('(')[1].split('currently')[0].strip())
        self.tstart = float(times.split('time =')[-1].split()[-3])
        self._gcconfig = self.initJob._gcconfig
        self._pnl = self.initJob._pnl
        self.gc = self.initJob.gc
        self.plume = self.initJob.plume
        self._jobDir = None
        self.tend = None

    def _setJobDir(self,jd):
        realDir = os.path.expandvars(jd)
        if (not(os.path.lexists(realDir))):
            raise Exception("jobDir %s does not exist" % realDir)
        else:
            self.initJob.jobDir = realDir
            self._jobDir = realDir
    def _getJobDir(self):
        if (self._jobDir is None):
            raise Exception("jobDir is not defined")
        return self._jobDir
    jobDir = property(fget=_getJobDir,fset=_setJobDir)

    def _setm(self,m):
        self.initJob.m = m
    def _getm(self):
        return self.initJob.m
    m = property(fget=_getm,fset=_setm)

    def _setn(self,n):
        self.initJob.n = n
    def _getn(self):
        return self.initJob.n
    n = property(fget=_getn,fset=_setn)

    def _settstart(self,t):
        self.initJob.tstart = t
    def _gettstart(self):
        return self.initJob.tstart
    tstart = property(fset=_settstart,fget=_gettstart)
    def _settend(self,t):
        self.initJob.tend = t
    def _gettend(self):
        return self.initJob.tend
    tend = property(fset=_settend,fget=_gettend)

    def _setname(self,n):
        self._name = n
        self.initJob.name = n
    def _getname(self):
        if (self._name  is None):
            raise Exception("self.name is not defined")
        return self._name
    name = property(fset=_setname,fget=_getname)

    def _setuse_plume(self,u):
        self.initJob.use_plume = u
    def _getuse_plume(self):
        return self.initJob.use_plume
    use_plume = property(fset=_setuse_plume,fget=_getuse_plume)

    def _resolveJob(self):
        #resolve the contained job
        self.initJob._resolveJob()

    def assertCanStage(self):
        self.initJob.assertCanStage()
        _BaseJob.assertCanStage(self)

    def run(self):
        self.initJob.run()

    def stage(self, genInput=True):
        self.assertCanStage()
        self._resolveJob()
        self.initJob.stage(genInput=False)
        if (genInput):
            cmd = self._genInput()
            _check_call(cmd)

    def _genInput(self):
        cmd =['nc_regrid']
        cmd.extend([self.inputNcFile,
                    self.initJob.inputfile,
                    self._inputNcTimeIndex,
                    self.m,self.n,
                    0,0,0,0,0,0,0,0,
                    self.initJob.kinbcw, 0,0,0])             
        cmd = [_fortran_style('nc_regrid', c) for c in cmd]
        return cmd

class GivenInputJob(_AtomicJob):
    pass

class CheckpointJob(_BaseJob):
    pass

def _check_call(cmd):
    retcode = subprocess.call(cmd)
    if (retcode != 0):
        raise Exception("Error running:\n       %s" % ' '.join(cmd))
    
    

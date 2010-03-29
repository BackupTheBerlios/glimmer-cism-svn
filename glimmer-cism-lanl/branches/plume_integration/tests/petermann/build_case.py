

#  This script is used to build the following:
#         namelist file for plume model (if required)
#         netcdf input file for glimmer-cism
#         config file for glimmer-cism

# It works by taking values from an input file, which is itself
# a python module, and overwriting the default values with the
# provided values.

class plume_namelist(dict):

    def __init__(self, otherDict):

        # Default values are defined here first.
        
        self['mixlayer'] = 'F'
        self['in_glimmer'] = 'F'
        self['restart'] = 'F'
        self['frazil']  = 'F' 
        self['nonlin']  = 'T'  
        self['horturb'] = 'F' 
        self['entrain'] = 'T'  
        self['entype']  = '1'       
        self['basmelt'] = 'T'  
        self['rholinear'] = 'T'  
        self['thermobar'] = 'F' 
        self['intrace']  = 'F' 
        self['vardrag']  = 'F' 
        self['topedit']  = 'F' 
        self['tangle']  = 'F' 
        self['negfrz']  = 'F' 
        self['use_min_plume_thickness'] = 'T'
        self['tottim']  = '000.0d0'
        self['outtim']  = '00.1d0'
        self['labtim']  = '00.1d0'
        self['snottim'] = '00.01d0'
        self['lnottim'] = '00.1d0'
        self['dt1']     = '50.0d0'              
        self['m_grid'] = 46         
        self['n_grid'] = 86          
        self['hx'] = '200.d0'      
        self['hy'] = '200.d0'
        self['gldep'] = '1000.d0'
        self['ifdep'] = '800.d0'    
        self['wcdep'] = '1000.d0'   
        self['plume_min_thickness'] = '10.0d0'
        self['bsmoothit'] = '00'
        self['salttop'] = '34.500d0'
        self['saltbot'] = '34.500d0'
        self['temptop'] = '-1.000d0'
        self['tempbot'] = '-1.000d0'
        self['phi'] = '+0.d0'                   
        self['ah'] = '0000.d0'    
        self['kh'] = '0000.d0'
        self['cdb'] = '2.5d-3'    
        self['cl'] = '1.775d-2'   
        self['ef'] = '5.0d-1'     
        self['nus'] = '+1.d0'     
        self['nbar'] = '1.0d3'    
        self['nice'] = '10'        
        self['seedtype'] = '2' 
        self['cseedfix'] = '1.0d-7'
        self['cinffix']  = '0.0d0'

        #next we overwrite default values with any values provided 
        self.update(otherDict)

    def produce_namelist_contents(self):
        
        # return a string which can be used as the contents of
        # a plume namelist file
        
        lines = []
        for (k,v) in self.items():
            lines.append(', %s = %s\n' % (k,v))
        return '&plume_nml\n %s /' % ' '.join(lines)


class glimmer_cims_config(object):

    def __init__(self, specificVals):
        
        self.vals = { 'parameters' : 
                      { 'geothermal' : -42.0e-3,
                        'default_flwa' : 4.6e-18,
                        'flow_factor' : 1},
                      'Petermann shelf' :
                      {  'air_temperature' : -5.0,
                         'accumulation_rate' : 10.0,
                         'eustatic_sea_level' : 0.0 },
                      'options' :
                      { 'flow_law' : 1,
                        'evolution' : 3,
                        'termperature' : 0,
                        'vertical_integration' : 1,
                        'marine_margin' : 0,
                        'topo_is_relaxed' : 1,
                        'slip_coeff' : 1,
                        'periodic_ew' : 0,
                        'periodic_ns' : 0,
                        'hotstart' : 0,
                        'basal_water' : 2,
                        'which_bmlt' : 1}
                      }
                        
                       
                        
                      
            
def main(config_file):
    print config_file


if __name__ == '__main__':
    main(sys.argv[1:])

    
    

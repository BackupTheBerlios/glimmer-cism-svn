#!/usr/bin/python

# run a glimmer-cism-plume simulation and checkpoint
# then run periodically

import build_case
import subprocess
import sys
import math

def checkpoint_run(config_file):

    base_config = build_case.get_config(config_file)

    # check for required checkpoint-related variables
    for s in ['checkpoint_duration',
              'checkpoint_tstart',
              'checkpoint_tend',
              'checkpoint_initial_input_style']:
        if not s in base_config:
            raise Exception('Can not find %s in %s' %
                            (s, config_file))
    cp_dur = base_config['checkpoint_duration']
    ts = base_config['checkpoint_tstart']
    te = base_config['checkpoint_tend']

    #how are we going to generate the first netcdf input file?
    init_input_style = base_config['checkpoint_initial_input_style']

    # file names that will need to be modified for each checkpoint run
    plume_output_base = base_config['gc_vals']['plume']['plume_output_file']
    plume_nl_fname_base = base_config['gc_vals']['plume']['plume_nl_file']
    gc_config_fname_base= base_config['gc_config_filename']
    gc_output_fname_base = base_config['gc_vals']['CF output']['name']
    gc_input_fname_base = base_config['gc_vals']['CF input']['name']
    gc_cf_default_title_base = base_config['gc_vals']['CF default']['title']
    
    n_checkpoints = int(math.ceil(float(te - ts)/cp_dur))
    num_slices_per_cp = int(math.floor(cp_dur/ float(
                        base_config['gc_vals']['time']['dt']))) - 1
    

    for i in range(1,n_checkpoints+1):

        if (i == 1):
            base_config['input_style'] = init_input_style
            base_config['nc_gen_input_args'][1] = insert_cp_index(gc_input_fname_base,i)
        else:
            base_config['input_style'] = 'regrid'
            base_config['nc_regrid_args'][0] = \
                insert_cp_index(gc_output_fname_base,i-1) #output file from last run
            base_config['nc_regrid_args'][2] = -1 # this means read last time slice,            

        base_config['nc_regrid_args'][1] = \
                insert_cp_index(gc_input_fname_base,i)    #file for input to next run

        
        base_config['gc_vals']['time']['tstart'] = ts
        base_config['gc_vals']['time']['tend'] =  min(te, ts + cp_dur)
        base_config['gc_vals']['plume']['plume_output_file'] = \
                   insert_cp_index(plume_output_base, i)
        base_config['gc_config_filename'] = \
                   insert_cp_index(gc_config_fname_base, i)
        base_config['plume_nl_filename'] = \
                   insert_cp_index(plume_nl_fname_base, i)
        base_config['gc_vals']['plume']['plume_nl_filename'] = \
                   insert_cp_index(plume_nl_fname_base,i)
        base_config['gc_vals']['CF input']['name'] = \
                   insert_cp_index(gc_input_fname_base,i)
        base_config['gc_vals']['CF output']['name']= \
                  insert_cp_index(gc_output_fname_base,i)
        base_config['gc_vals']['CF output']['start'] = ts
        base_config['gc_vals']['CF output']['stop'] = min(te,ts + cp_dur)
        base_config['gc_vals']['CF default']['title'] = \
                   '%s.cp%s' % (gc_cf_default_title_base,i)

        build_case.build_case(base_config)

        cmd = ['shelf_driver', base_config['gc_config_filename']]
        print('Running %s' % ' '.join(cmd))
        retval = subprocess.call(cmd)
        if (retval != 0):
            raise Exception('Error running:\n %s' % ' '.join(cmd))
         
        ts = ts + cp_dur

def insert_cp_index(stringbase,ind):
    l = stringbase.split('.')
    last = l[-1]
    l[-1] = '0' * (3-len(str(ind))) + str(ind)
    l.append(last)
    return '.'.join(l)


def main():
    if (len(sys.argv) < 2):
        raise Exception('Need to call with config file as last argument')
    
    config_filename = sys.argv[-1] #last argument
    checkpoint_run(config_filename)

if __name__ == '__main__':
    main()

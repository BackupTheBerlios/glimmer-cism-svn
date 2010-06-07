#!/usr/bin/python

# run a glimmer-cism-plume simulation and checkpoint
# the run periodically

import build_case
import subprocess
import sys
import math


def checkpoint_run(config_file):

    base_config = build_case.get_config(config_file)
    
    print(base_config)

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

    init_input_style = base_config['checkpoint_initial_input_style']
    
    plume_output_base = base_config['gc_vals']['plume']['plume_output_file']
    plume_nl_fname_base = base_config['gc_vals']['plume']['plume_nl_filename']
    gc_config_fname_base= base_config['gc_config_filename']
    gc_output_fname_base = base_config['gc_vals']['CF output']['name']
    gc_input_fname_base = base_config['gc_vals']['CF input']['name']
    
    n_checkpoints = int(math.ceil(float(te - ts)/cp_dur))
    num_slices_per_cp = int(math.floor(cp_dur/ float(
                        base_config['gc_vals']['time']['dt']))) - 1
    
    # do initial run
    base_config['input_style'] = init_input_style
    base_config['gc_vals']['time']['tstart'] = ts
    base_config['gc_vals']['time']['tend'] = ts + cp_dur
    base_config['gc_vals']['CF output']['name'] = \
                               insert_cp_index(gc_output_fname_base,1)
    
    build_case.build_case(base_config)
    
    cmd = ['shelf_driver', gc_config_fname_base]
    retcode = subprocess.call(cmd)
    if (retcode != 0):
        raise Exception('Error running:\n %s' % ' '.join(cmd))

    for i in range(2,n_checkpoints+1):

        ts = ts + cp_dur

        base_config['input_style'] = 'regrid'
        base_config['nc_regrid_args'][0] = \
                insert_cp_index(gc_output_fname_base,i-1) #input file
        base_config['nc_regrid_args'][1] = \
                insert_cp_index(gc_input_fname_base,i)    #outpt file
        base_config['nc_regrid_args'][2] = \
                -1 # this means read last time slice,
                   # replacing num_slices_per_cp
        
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
        base_config['gc_vals']['CF output']['stop'] = ts + cp_dur
        
        build_case.build_case(base_config)

        cmd = ['shelf_driver', base_config['gc_config_filename']]
        retval = subprocess.call(cmd)
        if (retcode != 0):
            raise Exception('Error running:\n %s' % ' '.join(cmd))
    
def insert_cp_index(str,ind):
    l = str.split('.')
    last = l[-1]
    l[-1] = '%s' % ind
    l.append(last)
    return '.'.join(l)


def main():
    if (len(sys.argv) < 2):
        raise Exception('Need to call with config file as last argument')
    
    config_filename = sys.argv[-1] #last argument
    checkpoint_run(config_filename)

if __name__ == '__main__':
    main()

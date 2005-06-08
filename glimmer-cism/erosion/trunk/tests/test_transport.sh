#!/bin/bash

# test transport code
# Magnus Hagdorn, June 2005

DGRID="20 10"
NTIME="10 20"
GRID_FACTOR="1 4"

# create input velo fields
python trans_velo.py

# create configuration files
plots=""
for dg in $DGRID
do
  dgm=`echo 1000*$dg | bc`
  num=`echo 1500/$dg+1 | bc`
  for nt in $NTIME
  do
    for gf in $GRID_FACTOR
    do
      cfg="$dg"km-$nt-$gf
      echo $cfg
      sed "s/NUM/$num/; s/DELTA_GRID/$dgm/; s/DGRID/$dg/; s/NTIME/$nt/; s/GRID_FACTOR/$gf/" trans_velo.config.in > trans_velo.$cfg.config
      echo trans_velo.$cfg.config | ../src/test_transport
      plots="$plots trans_velo.$cfg.nc"
    done
  done
done

echo "0.      0.
1500000.        1500000." > prof

plotProfile.py -vseds1 -t12.9 -pprof $plots trans_velo.ps

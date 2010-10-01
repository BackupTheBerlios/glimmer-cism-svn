#!/bin/bash

######################################################################
# Test script version 1: 30th Sep 2010 
# This script runs simple_glide and compares the output to a golden 
# standard file
#
# Currently the tests which work are:
#   
# 1. EISMINT-1 moving margin experiments 1, 2 & 3 
# 2. ho-other (hump) tests using Payne Price and Pattyn/Bocek/Johonson solvers
# 3. ISMIP-HOM a, c using Payne Price and Pattyn/Bocek/Johonson solvers 
# at 160km resolution   
# 
# Usage: 
#
# Invoke the test from the working directory i.e the directory containing 
# configure, src, etc.. 
# 
# run a test by typing: make test expt_type=***** expt_num=****
#
# For example: 
#
# An EISMINT-1 test:    make test expt_type=EISMINT-1  expt_num=2  (1,2,3) 
#   
# A hump test:          make test expt_type=ho-other   expt_num=PP (PBJ) 
#
# An ISMIP-HOM test:    make test expt_type=ISMIP-HOM  expt_num=a  (a,c)
######################################################################

# the script requires the compiler type and directory path. This comes from the Makefile 
prefix=$1 
compiler=$2

# the script uses these command line arguments input by the user
expt_type=$3
expt_number=$4

#location of the netCDF comparison tool
DIFFTOOL=$prefix/tests/util/compare

#location of a log file to print details of whether a test has passed or failed. Useful for the nightly build   
logfile=$prefix/test.log


#----------- start of subroutines --------------------------------------

# run EISMINT-1 moving margin test and compare output to golden standard

run_eismint_1 (){
    cd ${prefix}/src/fortran

    EISMINT="e1-mm."$expt_number".config" 
   
    if [ -f ${prefix}/tests/SIA/EISMINT-1/config_files/${EISMINT} ] ;then
	
	echo ${prefix}/tests/SIA/EISMINT-1/config_files/${EISMINT} > ${EISMINT}.fname
	
	${prefix}/src/fortran/simple_glide < ${prefix}/src/fortran/${EISMINT}.fname
	
	fname=$(echo ${EISMINT}|sed 's/config/nc/g') >> $logfile
	    
	if $DIFFTOOL $fname $prefix/tests/SIA/EISMINT-1/golden_std/$fname -a 1e-8 ; then

#	if $DIFFTOOL $fname /data/ggsrs/test_suite/glimmer-cism-lanl/EISMINT-1/short/golden_std/$fname -a 1e-8 ; then
	    echo "EISMINT-1" $fname "**TEST OK**" "using" $compiler "compiler">> $logfile
	    echo "EISMINT-1" $fname "**TEST OK**" "using" $compiler "compiler"
	else
	    echo "EISMINT-1" $fname "**TEST FAILED**" "using" $compiler "compiler" >> $logfile
	    echo "EISMINT-1" $fname "**TEST FAILED**" "using" $compiler "compiler"
	fi

    else 
	echo "unable to find" ${EISMINT} 
	exit
    fi	

}


# This follows the instructions in the original /tests/HO/ISMIP-HOM/README file
# and makes use of the existing python script runISMIP-HOM.py 
# The resolution and choice of solver (i.e. Payne Price) is hard coded: perhaps there is a better way to do this? 


run_ismip_hom () {

    cd ${prefix}/tests/HO/ISMIP-HOM/config_files

    #copy the original configuration file provided as instructed in the README file
  
    ISMIP="ishom."$expt_number".PP.config" 

    cp  $ISMIP "ishom."$expt_number".config"

    #hard code the resolution
    resolution=160

    python runISMIPHOM.py --exp=$expt_number --size=$resolution --run=/${prefix}/src/fortran/simple_glide

    fname="ishom."$expt_number"."$resolution"km.out.nc"

    if $DIFFTOOL $prefix/tests/HO/ISMIP-HOM/golden_std/$fname $prefix/tests/HO/ISMIP-HOM/config_files/output/$fname -a 1e-8 ; then
	echo "ISMIP-HOM" $ISMIP "**TEST OK**" "using compiler" $compiler >> $logfile
	echo "ISMIP-HOM" $ISMIP "**TEST OK**" "using" $compiler "compiler"
    else
	echo "ISMIP-HOM" $ISMIP "**TEST FAILED**" "using compiler" $compiler >> $logfile
	echo "ISMIP-HOM" $ISMIP "**TEST FAILED**" "using" $compiler "compiler"
    fi
}


run_ho_other (){
    cd ${prefix}/src/fortran

    HUMP="hump."$expt_number".config" 


    if [ -f ${prefix}/tests/HO/ho-other/config_files/${HUMP} ] ;then
	
	echo ${prefix}/tests/HO/ho-other/config_files/${HUMP} > ${HUMP}.fname
	
	${prefix}/src/fortran/simple_glide < ${prefix}/src/fortran/${HUMP}.fname
	
	fname=$(echo ${HUMP}|sed 's/config/out.nc/g') >> $logfile
	    
	if $DIFFTOOL $fname $prefix/tests/HO/ho-other/golden_std/$fname -a 1e-8 ; then
	    echo "HUMP" $fname "**TEST OK**" "using" $compiler "compiler">> $logfile
	    echo "HUMP" $fname "**TEST OK**" "using" $compiler "compiler"
	else
	    echo "HUMP" $fname "**TEST FAILED**" "using" $compiler "compiler" >> $logfile
	    echo "HUMP" $fname "**TEST FAILED**" "using" $compiler "compiler"
	fi

    else 
	echo "unable to find" ${HUMP} 
	exit
    fi	

}
# ----------end of subroutines ---------------------------

# ----------start of main code ---------------------------

if [ -f $logfile ]; then
    rm $logfile
fi

case $expt_type in

    EISMINT-1 ) 
	case $expt_number in
	    [1-3]) run_eismint_1 ;;
	    * ) echo "error:select a number between 1-3"  ;;
	esac;;

   
    ISMIP-HOM ) 
	case $expt_number in
	    [a-c]) run_ismip_hom ;;
	    * ) echo "error:select a letter between a-c  ";;
	esac;;
      
    ho-other ) 
	case $expt_number in
	    PP | PBJ ) run_ho_other ;;
	    * ) echo "error:select expt_num=PP or PBJ  ";;
	esac;;

    
    * ) echo "------------------------------------------------------------------------"
	echo "Run a test by typing: make test expt_type=? expt_num=?"
	
	echo "These tests currently work:" 

	echo "An EISMINT-1 test:    make test expt_type=EISMINT-1  expt_num=2  (1,2,3)" 
   
	echo "A hump test:          make test expt_type=ho-other   expt_num=PP (PBJ)" 

	echo "An ISMIP-HOM test:    make test expt_type=ISMIP-HOM  expt_num=a  (a,c)";;   
esac


# ----------end of main code ---------------------------










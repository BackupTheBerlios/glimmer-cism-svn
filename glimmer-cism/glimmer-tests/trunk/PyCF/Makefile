# Magnus Hagdorn
# testing PyCF system
# This should also give an idea what the programs can do.
# We assume that the test suite has been run so these data files can be used.

all:	test1 test2 test3

# Random Topography
test1:
	echo "Test 1"
	echo "Generating random topography"
	ran_topo.py --title 'ran topo test' --institution 'UoE' --comment 'defaults' ran_topo.nc
	echo "Plotting random topography"
	plotCFvar.py --colourmap=None -vtopg ran_topo.nc ran_topo.ps
	echo ""

# Ice sheet statistics
test2:	test2a
test2a:
	echo "Test 2a"
	echo "Plotting EISMINT statistics"
	plotCFstats.py -p eismint.prof ../EISMINT-1/e1-*.nc eismint1-stat.ps

# EISmint plot
test3:
	echo "Test 3"
	echo "Plotting EISMINT info"
	plotEISMINT.py -p eismint.prof -t 100. ../EISMINT-2/e2.h.nc e2.h.info.ps

# plotting a single location
test4:
	echo "Test 4"
	echo "Plotting ice thickness at two locations"
	plotSpot.py -vthk -i15 15 -i 20 20 ../EISMINT-2/e2.h.nc e2.h.spot.ps
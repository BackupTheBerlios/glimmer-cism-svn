# Magnus Hagdorn
# testing PyCF system
# This should also give an idea what the programs can do.
# We assume that the test suite has been run so these data files can be used.

all:	test1 test2 test3 test4 test5 test6 test7 test8

# Random Topography
test1:
	echo "Test 1"
	echo "Generating random topography"
	ran_topo.py --title 'ran topo test' --institution 'UoE' --comment 'defaults' ran_topo.nc
	echo "Plotting random topography"
	plotCFvar.py --colourmap=None -vtopg ran_topo.nc ran_topo.ps

# Ice sheet statistics
test2:	test2a test2b
test2a:
	echo "Test 2a"
	echo "Plotting EISMINT statistics"
	plotCFstats.py -p eismint.prof ../EISMINT-1/e1-*.nc eismint1-stat.ps
test2b:
	echo "Test 2b"
	echo "Plotting EIS stats with reconstructed ice extent."
	plotCFstats.py -p fenscan.prof --not_projected -f fenscan.data ../eis/fenscan.nc fenscan.stat.ps

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

# plotting CF variables
test5:	test5a test5b
test5a:
	echo "Test 5a"
	echo "Plotting ice surfaces at 6 times"
	plotCFvar.py --width=7.5 -vis -cthk --land -t -100. -t -80. -t -60. -t -40. -t -20 -t 0 ../eis/fenscan.nc fenscan.is.ps
test5b:
	echo "Test 5b"
	echo "Plotting different variables"
	plotCFvar.py --width=7.5 -vacab -vcony -vslc -vtopg -t-17. ../eis/fenscan.nc fenscan.vars.ps

# plotting profiles
test6:	test6a test6b
test6a:
	echo "Test 6a"
	echo "Plotting some profiles"
	plotProfile.py --not_projected -p fenscan.prof -vis -vbtemp -vslc --pmt -t -20 ../eis/fenscan.nc fenscan.profs.ps
test6b:
	echo "Test 6b"
	echo "Plotting 2D temperature prof"
	plotProfile.py --not_projected -p fenscan.prof -vtemp  --pmt -t -20 ../eis/fenscan.2000a.nc fenscan.temp.ps

# plotting RSL curves
test7:
	echo "Test 7"
	echo "Plotting RSL curves"
	plotRSL.py ../eis/fenscan.nc fenscan.rsl.ps

# plotting EIS forcing
test8:
	echo "Test 8"
	echo "Plotting EIS forcing"
	plotEISvar.py --land -cthk -t -17. -vis --ela=../eis/mb.data --temp=../eis/s10.data --slc=../eis/specmap.data ../eis/fenscan.nc fenscan.forcing.ps

clean:
	rm -f *.ps ran_topo.nc

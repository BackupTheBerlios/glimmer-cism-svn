
runs:	isos.llfa.nc isos.llra.nc isos.elfa.nc isos.elra.nc
plots:	spots.ps profile.ps

.config.nc:
	glide_launch.py $<

spots.ps:	isos.llfa.nc isos.llra.nc isos.elfa.nc isos.elra.nc
	plotSpot.py -i 15 15 -i 15 22 -i 15 27 -vtopg $^ $@

profile.ps:	isos.llfa.nc isos.llra.nc isos.elfa.nc isos.elra.nc
	plotProfile.py -t 18  -v topg -p ../PyCF/eismint.prof $^ $@

.SUFFIXES: .config .nc .ps

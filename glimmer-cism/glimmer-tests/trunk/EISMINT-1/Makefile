
runs:	e1-fm.1.nc e1-fm.2.nc e1-fm.3.nc e1-mm.1.nc e1-mm.2.nc e1-mm.3.nc
plots:	e1-fm.1.ps e1-fm.2.ps e1-fm.3.ps e1-mm.1.ps e1-mm.2.ps e1-mm.3.ps

.config.nc:
	glide_launch.py $<

.nc.ps:
	plotEISMINT.py -t 200 --pmt  -p ../PyCF/eismint.prof $< $@

.SUFFIXES: .config .nc .ps

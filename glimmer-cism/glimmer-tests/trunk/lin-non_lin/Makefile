RUNS = l-5a.nc nl-5a.nc l-10a.nc nl-10a.nc l-20a.nc nl-20a.nc l-30a.nc nl-30a.nc

runs:	$(RUNS)
plots:	e2-comp.ps e2-comp-t.ps e2-btemp.ps e2-bvel.ps

e2-comp.ps:	$(RUNS)
	plotEISMINT2stats.py -e H $^ $@

e2-comp-t.ps:	$(RUNS)
	plotCFstats.py -e -m -p ../PyCF/eismint.prof $^ $@

e2-btemp.ps:	$(RUNS)
	plotCFvar.py --width 5 -T -1 -cthk -vtemp -l -1 $^ $@

e2-bvel.ps:	$(RUNS)
	plotCFvar.py --width 5 -T -1 -cthk -vvel -l -1 $^ $@

.config.nc:
	glide_launch.py $<
	bash move_data.sh picard_info.data $*.piccard

.SUFFIXES:	.config .nc .ps
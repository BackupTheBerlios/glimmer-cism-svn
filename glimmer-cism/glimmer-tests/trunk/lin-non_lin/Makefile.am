runs = l-5a.nc nl-5a.nc l-10a.nc nl-10a.nc l-20a.nc nl-20a.nc l-30a.nc nl-30a.nc

if BUILD_GRAPHICS
  graphics = e2-comp.ps e2-comp-t.ps e2-btemp.ps e2-bvel.ps
endif
BUILT_SOURCES = $(runs) $(graphics)

EXTRA_DIST = l-10a.config  l-20a.config  l-5a.config  nl-10a.config  nl-20a.config  nl-5a.config l-30a.config  nl-30a.config move_data.sh


e2-comp.ps:	$(runs)
	@GLIMMER_PREFIX@/bin/plotEISMINT2stats.py -e H $^ $@

e2-comp-t.ps:	$(runs)
	@GLIMMER_PREFIX@/bin/plotCFstats.py -e -m -p ../PyCF/eismint.prof $^ $@

e2-btemp.ps:	$(runs)
	@GLIMMER_PREFIX@/bin/plotCFvar.py --width 5 -T -1 -cthk -vtemp -l -1 $^ $@

e2-bvel.ps:	$(runs)
	@GLIMMER_PREFIX@/bin/plotCFvar.py --width 5 -T -1 -cthk -vvel -l -1 $^ $@

.config.nc:
	@GLIMMER_PREFIX@/bin/glide_launch.py $<
	bash move_data.sh picard_info.data $*.piccard

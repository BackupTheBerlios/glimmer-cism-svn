if BUILD_GRAPHICS
  graphics = spots.ps profile.ps
endif
BUILT_SOURCES = isos.llfa.nc isos.llra.nc isos.elfa.nc isos.elra.nc $(graphics)

EXTRA_DIST = isos.elfa.config  isos.elra.config  isos.llfa.config  isos.llra.config

MOSTLYCLEANFILES = $(graphics)
CLEANFILES = $(BUILT_SOURCES)

spots.ps:	isos.llfa.nc isos.llra.nc isos.elfa.nc isos.elra.nc
	@GLIMMER_PREFIX@/bin/plotSpot.py -i 15 15 -i 15 22 -i 15 27 -vtopg $^ $@

profile.ps:	isos.llfa.nc isos.llra.nc isos.elfa.nc isos.elra.nc
	@GLIMMER_PREFIX@/bin/plotProfile.py -t 18  -v topg -p ../PyCF/eismint.prof $^ $@

include $(top_srcdir)/make.rules

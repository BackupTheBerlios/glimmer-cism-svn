graphics = topo.eps is.eps btemp.eps bvel.eps bvelp.eps extent.eps prof_is.eps prof_vars.eps prof_temp.eps
EXTRA_DIST = $(graphics)
BUILT_SOURCES = $(graphics)
CLEANFILES = $(BUILT_SOURCES)


topo.ps:	@GLIMMER_TESTS@/eis/fenscan.nc
		plotCFvar.py -vtopg $< $@

is.ps:		@GLIMMER_TESTS@/eis/fenscan.nc
		plotCFvar.py -vis -cthk -iis -t -17 --land $< $@

btemp.ps:	@GLIMMER_TESTS@/eis/fenscan.nc
		plotCFvar.py -vbtemp -cthk -gbvel --pmt -t -17 --land $< $@

bvel.ps:	@GLIMMER_TESTS@/eis/fenscan.nc
		plotCFvar.py -vbvel --urg=20 60 -t -17  $< $@

bvelp.ps:	@GLIMMER_TESTS@/eis/fenscan.nc
		plotCFvar.py -vbvel -cthk -t -17 --land -p @GLIMMER_TESTS@/PyCF/fenscan.prof --not_p $< $@

extent.ps:	@GLIMMER_TESTS@/eis/fenscan.nc @GLIMMER_TESTS@/eis/fenscan-gthf.nc 
		plot_extent.py -t -17 $^ $@

prof_is.ps:	@GLIMMER_TESTS@/eis/fenscan.nc @GLIMMER_TESTS@/eis/fenscan-gthf.nc
		plotProfile.py -vis  -t -17 -p @GLIMMER_TESTS@/PyCF/fenscan.prof --not_p $^ $@

prof_vars.ps:	@GLIMMER_TESTS@/eis/fenscan.nc @GLIMMER_TESTS@/eis/fenscan-gthf.nc
		plotProfile.py -vis -vbvel -vbmlt -t -17 -p @GLIMMER_TESTS@/PyCF/fenscan.prof --not_p $^ $@

prof_temp.ps:	@GLIMMER_TESTS@/eis/fenscan.2000a.nc
		plotProfile.py -vtemp  -t -18 -p @GLIMMER_TESTS@/PyCF/fenscan.prof --not_p $^ $@

include $(top_srcdir)/doc/make_doc.rules

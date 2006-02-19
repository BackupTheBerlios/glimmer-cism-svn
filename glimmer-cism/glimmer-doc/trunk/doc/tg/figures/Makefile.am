if HAVE_RSL_DB
  rsl_graphics = rsl_loc.eps rsl_curves.eps rsl_data.eps rsl_res.eps rsl_hist.eps
endif

graphics = topo.eps is.eps thk.eps btemp.eps bvel.eps bvelp.eps extent.eps prof_is.eps prof_vars.eps prof_temp.eps prof_temp2.eps \
           stats.eps spot_temp.eps spot_thk.eps time_dist1.eps time_dist2.eps streams1.eps streams2.eps diff.eps soft.eps\
           $(rsl_graphics)

EXTRA_DIST = $(graphics)
BUILT_SOURCES = $(graphics)

topo.ps:	@GLIMMER_TESTS@/eis/fenscan.nc
		plotCFvar.py -vtopg $< $@

is.ps:		@GLIMMER_TESTS@/eis/fenscan.nc
		plotCFvar.py -vis -cthk -iis -t -17 --land $< $@

thk.ps:		@GLIMMER_TESTS@/eis/fenscan.nc
		plotCFvar.py -vthk -cthk -t-17 --land --contours=0/4000/500/1000 $< $@

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

prof_temp2.ps:	@GLIMMER_TESTS@/eis/fenscan.2000a.nc @GLIMMER_TESTS@/eis/fenscan-gthf.2000a.nc
		plot3DProfiles.py -vtemp  -t -18 -p @GLIMMER_TESTS@/PyCF/fenscan.prof --not_p $^ $@

stats.ps:	@GLIMMER_TESTS@/eis/fenscan.nc @GLIMMER_TESTS@/eis/fenscan-gthf.nc
		plotCFstats.py -mi -p @GLIMMER_TESTS@/PyCF/fenscan.prof --not_p -f @GLIMMER_TESTS@/PyCF/fenscan.data $^ $@

diff.ps:	@GLIMMER_TESTS@/eis/fenscan.nc @GLIMMER_TESTS@/eis/fenscan-gthf.nc
		plotCFdiff.py $^ $@

spot_temp.ps:	@GLIMMER_TESTS@/eis/fenscan.2000a.nc
		plotSpot.py --ij 57 42 --ij 77 76 -vtemp  -t-18 $< $@

spot_thk.ps:	@GLIMMER_TESTS@/eis/fenscan.nc @GLIMMER_TESTS@/eis/fenscan-gthf.nc
		plotSpot.py --ij 77 76 -vthk $^ $@

time_dist1.ps:	@GLIMMER_TESTS@/eis/fenscan.nc
		plotProfileTS.py -t -30 0 -vbtemp --pmt -cthk -p @GLIMMER_TESTS@/PyCF/fenscan.prof --not_p $^ $@

time_dist2.ps:	@GLIMMER_TESTS@/eis/fenscan.nc
		plotProfileTS.py -e stages --profvar is -vbtemp --pmt -cthk -p @GLIMMER_TESTS@/PyCF/fenscan.prof --not_p $^ $@

rsl_loc.ps:	@GLIMMER_TESTS@/eis/fenscan.nc
		plotRSLloc.py $< $@

rsl_curves.ps:	@GLIMMER_TESTS@/eis/fenscan.nc @GLIMMER_TESTS@/eis/fenscan-gthf.nc
		plotRSL.py $^ $@

rsl_data.ps:	@GLIMMER_TESTS@/eis/fenscan.nc  @GLIMMER_TESTS@/eis/fenscan-gthf.nc
		plotRSLdata.py --id=103 $^ $@

rsl_res.ps:	@GLIMMER_TESTS@/eis/fenscan.nc
		plotRSLres.py $< $@

rsl_hist.ps:	@GLIMMER_TESTS@/eis/fenscan.nc
		plotRSLhist.py $< $@

streams1.ps:	@GLIMMER_TESTS@/eis/fenscan.nc
		plotStreams.py --legend -t-10 --deltat=10000. $< $@

streams2.ps:	@GLIMMER_TESTS@/eis/fenscan.nc
		plotStreams.py -t-10 --deltat=10000. --velocity $< $@

soft1.nc:	@GLIMMER_TESTS@/eis/fenscan.nc
		construct_field.py --title="1 patch" -vsoft -o 1.e-2 -i 5.e-2 -p 4.5 62.5 -p 5 59 -p 7.5 58 -p 10 59 -p 11 58 -p 8 57.5 -p 5 58 -p 2.5 62 $< $@

soft2.nc:	soft1.nc
		construct_field.py  --title="2 patches" -vsoft -i 5.e-2 -p 22.5 58 -p 26 55 -p 28 55.5 -p 23 58.8 $< $@

soft.ps:	soft1.nc soft2.nc
		plotCFvar.py --width=9 -vsoft $^ $@

include $(top_srcdir)/doc/make_doc.rules

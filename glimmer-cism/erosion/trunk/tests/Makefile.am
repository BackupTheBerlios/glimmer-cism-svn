TEST_CONFIGS = erosion0.config erosion1.config erosion2.config erosion3.config erosion4.config \
               er-transport0.config er-transport1.config er-transport2.config er-transport3.config er-transport4.config fenscan.config

EXTRA_DIST = trans_velo.config.in trans_velo.py $(TEST_CONFIGS)

if RUN_TESTS
  RUNS = $(subst config,nc,$(TEST_CONFIGS))
endif

BUILT_SOURCES = $(RUNS)
CLEANFILES = $(BUILT_SOURCES)

fenscan.nc:	fenscan.config
		glide_launch.py -m ../src/eis_erosion $<

.config.nc::
		glide_launch.py -m ../src/simple_erosion $<

include $(top_srcdir)/doc/make_doc.rules

TEST_CONFIGS = erosion0.config erosion1.config erosion2.config erosion3.config erosion4.config \
               er-transport0.config er-transport1.config er-transport2.config er-transport3.config er-transport4.config

EXTRA_DIST = topo.nc $(TEST_CONFIGS)

data = $(subst config,nc,$(TEST_CONFIGS))

CLEANFILES = $(data) $(graphics) *.jid $(BUILT_SOURCES)

MODEL_BINARY=$(top_builddir)/tests/EISMINT/src/simple_erosion

include $(top_srcdir)/extra_rules.am

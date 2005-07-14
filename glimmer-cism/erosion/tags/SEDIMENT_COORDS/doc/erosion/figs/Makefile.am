graphics = ice_velo.pstex ice_velo.pstex_t erosion_layers.pstex erosion_layers.pstex_t plot_basal.eps
EXTRA_DIST = ice_velo.fig erosion_layers.fig $(graphics)
BUILT_SOURCES = $(graphics)
CLEANFILES = $(BUILT_SOURCES)

plot_basal.ps:	plot_basal.py
	$(PYTHON) plot_basal.py

include $(top_srcdir)/doc/make_doc.rules

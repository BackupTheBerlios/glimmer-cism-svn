if BUILD_GRAPHICS
  graphics = 
endif
BUILT_SOURCES = fenscan.nc $(graphics)

EXTRA_DIST = fenscan.config europe.20km.nc mb.data s10.data specmap.data

MOSTLYCLEANFILES = $(graphics)
CLEANFILES = $(BUILT_SOURCES)

include $(top_srcdir)/make.rules


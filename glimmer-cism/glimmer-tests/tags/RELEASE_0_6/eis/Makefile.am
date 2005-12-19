if BUILD_GRAPHICS
  graphics = 
endif
BUILT_SOURCES = fenscan.nc fenscan-gthf.nc $(graphics)

EXTRA_DIST = fenscan.config fenscan-gthf.config europe.20km.nc mb2.data temp-exp.model specmap.data

MOSTLYCLEANFILES = $(graphics)
CLEANFILES = $(BUILT_SOURCES)

include $(top_srcdir)/make.rules


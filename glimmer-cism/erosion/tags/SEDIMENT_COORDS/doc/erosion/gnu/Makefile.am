graphics = stress-strain.pslatex stress-strain-fit.pslatex

EXTRA_DIST = stress-strain.gp fit.gp stress-strain-fit.gp stress-strain.data $(graphics)
 
BUILT_SOURCES = $(graphics)
CLEANFILES = $(BUILT_SOURCES)

stress-strain.gp:	fit.gp	stress-strain.data
stress-strain-fit.gp:	fit.gp	stress-strain.data

include $(top_srcdir)/doc/make_doc.rules

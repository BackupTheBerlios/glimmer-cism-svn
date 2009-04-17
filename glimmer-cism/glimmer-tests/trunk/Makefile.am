ACLOCAL_AMFLAGS = -I $(GLIMMER_PREFIX)/share/aclocal

SUBDIRS = EISMINT-1 EISMINT-2 lin-non_lin eis isos PyCF benchmarks

graphics::
		list='$(SUBDIRS)'; for dir in $$list; do \
		  (cd $$dir && $(MAKE) $@); \
		done

data::
		list='$(SUBDIRS)'; for dir in $$list; do \
		  (cd $$dir && $(MAKE) $@); \
		done


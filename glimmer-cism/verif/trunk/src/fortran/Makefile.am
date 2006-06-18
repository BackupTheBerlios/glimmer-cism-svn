bin_PROGRAMS = verif_glide
lib_LIBRARIES = libverif.a

variable_DATA = verif_vars.def
variabledir = $(datadir)/glimmer
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak
BUILT_SOURCES = verif_io.F90
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak

libverif_a_SOURCES = verif.F90 verifBC.F90 verifD.F90 verif_io.F90
verif_glide_SOURCES = verif_glide.F90
verif_glide_LDADD = -lglide -lglimmer -lslap libverif.a

#some special rules for automatically generating ncdf src code
verif_io.F90:	verif_vars.def @GLIMMER_PREFIX@/share/glimmer/ncdf_template.F90.in
		$(PYTHON) @GLIMMER_PREFIX@/bin/generate_ncvars.py $^
f90_dep.mak:	$(verif_glide_SOURCES) $(libverif_a_SOURCES)
		$(PYTHON) @GLIMMER_PREFIX@/bin/f90_dependencies.py -o $@ $^

-include f90_dep.mak

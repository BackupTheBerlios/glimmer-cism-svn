lib_LIBRARIES = liberosion.a
bin_PROGRAMS  = simple_erosion eis_erosion
check_PROGRAMS = 

variable_DATA = erosion_vars.def
variabledir = $(datadir)/glimmer
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS = *.mod
EXTRA_DIST = f90_dep.mak
BUILT_SOURCES = erosion_io.f90
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak

liberosion_a_SOURCES = erosion.f90 erosion_types.f90 erosion_setup.f90 erosion_io.f90 erosion_version.f90

simple_erosion_SOURCES = simple_erosion.f90
simple_erosion_LDADD = liberosion.a -lclimate -lglide -lglimmer -lslap

eis_erosion_SOURCES = eis_erosion.f90
eis_erosion_LDADD = liberosion.a -lclimate -lglide -lglimmer -lslap

#some special rules for automatically generating ncdf src code
erosion_io.f90:	erosion_vars.def @GLIMMER_PREFIX@/share/glimmer/ncdf_template.f90.in
		$(PYTHON) @GLIMMER_PREFIX@/bin/generate_ncvars.py $^

f90_dep.mak:	*.f90 $(BUILT_SOURCES)
		$(PYTHON) @GLIMMER_PREFIX@/bin/f90_dependencies.py -o $@ $^

-include f90_dep.mak

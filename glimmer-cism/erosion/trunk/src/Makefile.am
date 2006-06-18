lib_LIBRARIES = liberosion.a
bin_PROGRAMS  = simple_erosion eis_erosion
check_PROGRAMS = test_rk4 test_geometry test_advect test_interpolate2d

variable_DATA = erosion_vars.def
variabledir = $(datadir)/glimmer
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak
BUILT_SOURCES = erosion_io.F90
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak

liberosion_a_SOURCES = erosion.F90 erosion_types.F90 erosion_setup.F90 erosion_io.F90 erosion_version.F90 \
                       erosion_advect.F90 geometry.F90 runge_kutta.F90  erosion_transport.F90 \
                       glimmer_interpolate2d.F90 erosion_sediment.F90 advect_2ndmo.F90  erosion_transport_type.F90  

simple_erosion_SOURCES = simple_erosion.F90
simple_erosion_LDADD = liberosion.a -lclimate -lglide -lglimmer -lslap

eis_erosion_SOURCES = eis_erosion.F90
eis_erosion_LDADD = liberosion.a -lclimate -lglide -lglimmer -lslap

test_interpolate2d_SOURCES = test_interpolate2d.F90 geometry.F90 glimmer_interpolate2d.F90
test_interpolate2d_LDADD = -lglimmer

test_geometry_SOURCES = test_geometry.F90 geometry.F90

test_rk4_SOURCES = test_rk4.F90 runge_kutta.F90

test_advect_SOURCES = test_advect.F90 erosion_advect.F90 runge_kutta.F90 glimmer_interpolate2d.F90
test_advect_LDADD = -lglimmer

#some special rules for automatically generating ncdf src code
erosion_io.F90:	erosion_vars.def @GLIMMER_PREFIX@/share/glimmer/ncdf_template.F90.in
		$(PYTHON) @GLIMMER_PREFIX@/bin/generate_ncvars.py $^

f90_dep.mak:	$(liberosion_a_SOURCES) $(simple_erosion_SOURCES) $(eis_erosion_SOURCES) $(test_transport_SOURCES) \
                $(test_interpolate2d_SOURCES) $(test_geometry_SOURCES) $(test_rk4_SOURCES) $(test_advect_SOURCES)
		$(PYTHON) @GLIMMER_PREFIX@/bin/f90_dependencies.py -o $@ $^

-include f90_dep.mak

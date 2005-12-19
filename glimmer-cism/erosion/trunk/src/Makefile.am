lib_LIBRARIES = liberosion.a
bin_PROGRAMS  = simple_erosion eis_erosion
check_PROGRAMS = test_rk4 test_geometry test_advect test_transport test_interpolate2d

variable_DATA = erosion_vars.def
variabledir = $(datadir)/glimmer
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak
BUILT_SOURCES = erosion_io.f90
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak

liberosion_a_SOURCES = erosion.f90 erosion_types.f90 erosion_setup.f90 erosion_io.f90 erosion_version.f90 \
                       erosion_advect.f90 geometry.f90 runge_kutta.f90  erosion_transport.f90 \
                       glimmer_interpolate2d.f90 erosion_sediment.f90 advect_2ndmo.f90  erosion_transport_type.f90  

simple_erosion_SOURCES = simple_erosion.f90
simple_erosion_LDADD = liberosion.a -lclimate -lglide -lglimmer -lslap

eis_erosion_SOURCES = eis_erosion.f90
eis_erosion_LDADD = liberosion.a -lclimate -lglide -lglimmer -lslap

test_transport_SOURCES = test_transport.f90
test_transport_LDADD = liberosion.a -lclimate -lglide -lglimmer -lslap

test_interpolate2d_SOURCES = test_interpolate2d.f90 geometry.f90 glimmer_interpolate2d.f90
test_interpolate2d_LDADD = -lglimmer

test_geometry_SOURCES = test_geometry.f90 geometry.f90

test_rk4_SOURCES = test_rk4.f90 runge_kutta.f90

test_advect_SOURCES = test_advect.f90 erosion_advect.f90 runge_kutta.f90 glimmer_interpolate2d.f90
test_advect_LDADD = -lglimmer

#some special rules for automatically generating ncdf src code
erosion_io.f90:	erosion_vars.def @GLIMMER_PREFIX@/share/glimmer/ncdf_template.f90.in
		$(PYTHON) @GLIMMER_PREFIX@/bin/generate_ncvars.py $^

f90_dep.mak:	$(liberosion_a_SOURCES) $(simple_erosion_SOURCES) $(eis_erosion_SOURCES) $(test_transport_SOURCES) \
                $(test_interpolate2d_SOURCES) $(test_geometry_SOURCES) $(test_rk4_SOURCES) $(test_advect_SOURCES)
		$(PYTHON) @GLIMMER_PREFIX@/bin/f90_dependencies.py -o $@ $^

-include f90_dep.mak

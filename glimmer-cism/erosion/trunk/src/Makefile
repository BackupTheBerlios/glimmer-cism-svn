# Makefile for building the example climate driver
# Magnus Hagdorn, November 2004

# setup compilation flags: tell the compiler where to look for f90 module files
GLIMMER_FFLAGS = -I$(GLIMMER_PREFIX)/include -I$(GNU)/include
# setup linker flags: tell the linker where to look for libraries and which libraries to link
GLIMMER_LIBRARIES = -L$(GNU)/lib -L$(GLIMMER_PREFIX)/lib -lglide -lglimmer -lslap -lnetcdf
# extra compiler flags, e.g. for setting debug or profiler flags
FFLAGS = -C -std95 -traceback -g
# extra linker flags
LD_FLAGS =-C -std95 -traceback -g

AUTOGEN = erosion_io.f90
EROSION = erosion.o erosion_types.o erosion_setup.o erosion_io.o
PROGS = simple_erosion eis_erosion
TESTS = test_rk4 test_geometry test_advect test_integrate2d

all:		liberosion.a $(PROGS)
test:		$(TESTS)

liberosion.a:	$(EROSION)
		ar rcv $@ $^
		if [ "$(RANLIB)" != "" ] ; then $(RANLIB) $@ ; fi

eis_erosion:	eis_erosion.o liberosion.a
		$(FC) -o $@ $^ $(LD_FLAGS) -lclimate $(GLIMMER_LIBRARIES)

simple_erosion:	simple_erosion.o liberosion.a
		$(FC) -o $@ $^ $(LD_FLAGS) -lclimate $(GLIMMER_LIBRARIES) 

test_integrate2d:	test_integrate2d.o erosion_integrate2d.o sparse.o geometry.o glimmer_coordinates.o
		$(FC) -o $@ $^ $(LD_FLAGS)

test_geometry:	test_geometry.o geometry.o
		$(FC) -o $@ $^ $(LD_FLAGS)

test_rk4:	test_rk4.o runge_kutta.o
		$(FC) -o $@ $^ $(LD_FLAGS)

test_advect:	test_advect.o erosion_advect.o runge_kutta.o
		$(FC) -o $@ $^ $(LD_FLAGS) -L$(GNU)/lib -lnetcdf

# implicit rule, telling the make system how to build an object file from a .f90 file
.f90.o:
		$(FC) -fpp $(FFLAGS) $(GLIMMER_FFLAGS) -c $< -o $@

erosion_io.f90:	erosion_vars.def
		generate_ncvars.py $< $(GLIMMER_PREFIX)/share/glimmer/ncdf_template.f90.in


# clean up
clean : 
	rm -f *.o *.mod f90_dep.mak $(AUTOGEN) $(PROGS) $(LIBS)
	for dir in $(SUBDIRS); do \
	  $(MAKE) -C $$dir clean; \
	done

f90_dep.mak: *.f90 $(AUTOGEN)
	f90_dependencies.py -o $@ $^

include f90_dep.mak

.SUFFIXES:	.f .f90 .o

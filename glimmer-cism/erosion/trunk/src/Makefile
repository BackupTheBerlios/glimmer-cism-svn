# Makefile for building the example climate driver
# Magnus Hagdorn, November 2004

# setup compilation flags: tell the compiler where to look for f90 module files
GLIMMER_FFLAGS = -I$(GLIMMER_PREFIX)/mod -I$(GNU)/include
# setup linker flags: tell the linker where to look for libraries and which libraries to link
GLIMMER_LIBRARIES = -L$(GNU)/lib -L$(GLIMMER_PREFIX)/lib -lglide -lglimmer -lslap -lnetcdf
# extra compiler flags, e.g. for setting debug or profiler flags
FFLAGS = 
# extra linker flags
LD_FLAGS =

AUTOGEN = erosion_io.f90
EROSION = erosion.o erosion_types.o erosion_setup.o erosion_io.o

liberosion.a:	$(EROSION)
		ar rcv $@ $^
		if [ "$(RANLIB)" != "" ] ; then $(RANLIB) $@ ; fi


simple_erosion:	simple_erosion.o liberosion.a
		f95 -o $@ $^ $(LD_FLAGS) $(GLIMMER_LIBRARIES) -lclimate 

# implicit rule, telling the make system how to build an object file from a .f90 file
.f90.o:
		f95 -fpp $(FFLAGS) $(GLIMMER_FFLAGS) -c $< -o $@

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

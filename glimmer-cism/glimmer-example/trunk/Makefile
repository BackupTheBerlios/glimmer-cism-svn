# Makefile for building the example climate driver
# Magnus Hagdorn, November 2004

# setup compilation flags: tell the compiler where to look for f90 module files
GLIMMER_FFLAGS = `glimmer-config --fcflags`
# setup linker flags: tell the linker where to look for libraries and which libraries to link
GLIMMER_LIBRARIES = `glimmer-config --libs` -lglide -lglimmer -lslap -lnetcdf

# set variables FCFLAGS and LDFLAGS for extra libraries, e.g. netCDF

# how to build the example executable
example:	example.o
		$(FC) -o $@ $< $(LDFLAGS) $(GLIMMER_LIBRARIES)

# implicit rule, telling the make system how to build an object file from a .f90 file
.f90.o:
		$(FC) -fpp $(FCFLAGS) $(GLIMMER_FFLAGS) -c $< -o $@


.SUFFIXES:	.f90 .o

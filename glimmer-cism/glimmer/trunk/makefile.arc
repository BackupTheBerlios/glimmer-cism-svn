#
# Architecture specific defines for genie compilation - this version for
# Linux 
#
MACHINE=Linux
#
# C compiler to use
#
CC=gcc
#
# Fortran compiler to use
#
F77=ifc
#
# Compiler to use for final linkage stage
#
F77_LD=ifc
#
# Compilation flags for C compiler
#
CCFLAGS= -DIFC -g -DLINUX
# 
# compilation flags for fortran compiler
#
F77FLAGS=-fpp -w
#F77FLAGSR8=-r8
#
# final stage loader flags
#
LD_FLAGS=
#
# locations of wishx and xqplot.tcl
# they are in fact just empty files
# they are on cvs in the genie-main/inputdata directory.
# replace /home/ggdjl/genie/ with your genie code directory.
LOCFLAGS=-DWISHX=\"/home/ggdjl/genie/genie-main/inputdata/wishx\" -DXQPLOTTCL=\"/home/ggdjl/genie/genie-main/inputdata/xqplot.tcl\" -DXFONT=\"/home/ggdjl/genie/genie-main/inputdata/fonts\"
#
# How to get ordered list of objects to put in a library
#
ORDER_OBJECTS=*.o
#
# ranlib command, leave blank for no ranlib
#
RANLIB=ranlib
#
# *******GRAPHICS*******
# Libraries needed for final X (or other graphics) linkage
# and compile-time options for graphics. 
#(1) is usual for graphics
#(2) is for no graphics
#
#(1)
#XLIB=-L/usr/X11R6/lib -lX11
#DOPTS=-Dlgraph
#(2)
XLIB=
DOPTS=
#
# Libraries needed for final fortran linkage
#
FORTRANLIBS=
#
# netcdf library
# This may have to be installed from:
# http://www.unidata.ucar.edu/packages/netcdf/index.html
#
NETCDF=-lnetcdf_ifc

.c.o:
	$(CC) -c $(CCFLAGS) $*.c

.f.o:
	$(F77) -c $(F77FLAGS) $*.f

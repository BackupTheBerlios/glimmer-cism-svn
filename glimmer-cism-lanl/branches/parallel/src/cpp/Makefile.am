if TRILINOS

lib_LIBRARIES = libcalltrilinos.a

variabledir = $(pkgdatadir)
dist_variable_DATA = $(variable_DATA)

#nodist_include_HEADERS := $(wildcard *.mod)
BUILT_SOURCES = 
CLEANFILES = $(BUILT_SOURCES) 

libcalltrilinos_a_SOURCES = solveTriadMatrixWithTrilinos.cpp \
                            trilinosLinearSolver.cpp \
                            matrixInterface.cpp 
libcalltrilinos_a_DEPENDENCIES = matrixInterface.hpp

ALL_SRC= $(libcalltrilinoa_a_SOURCES)

AM_FCFLAGS = $(CXXFLAGS) $(NETCDF_CPPFLAGS) $(HDF5_CPPFLAGS)
AM_CXXFLAGS = $(CXXFLAGS)

LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) 

AM_CFLAGS=$(CFORTRAN_CFLAGS)

endif TRILINOS

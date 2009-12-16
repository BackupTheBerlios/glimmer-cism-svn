if TRILINOS

lib_LIBRARIES = libcalltrilinos.a
#bin_PROGRAMS  = helloEpetra

variabledir = $(pkgdatadir)
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS := $(wildcard *.mod)
BUILT_SOURCES = 
CLEANFILES = $(BUILT_SOURCES) *.mod

libcalltrilinos_a_SOURCES = call_trilinos.cpp

#helloEpetra_SOURCES = helloCXX.cpp interface_modCXX.F90 helloF.F90 
#helloEpetra_LDADD = $(BLAS_LIBS) $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

#ALL_SRC=$(helloEpetra_SOURCES) $(libcalltrilinoa_a_SOURCES)
ALL_SRC= $(libcalltrilinoa_a_SOURCES)

AM_FCFLAGS = $(CXXFLAGS) $(NETCDF_CPPFLAGS) $(HDF5_CPPFLAGS)
AM_CXXFLAGS = $(CXXFLAGS)

LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) 

AM_CFLAGS=$(CFORTRAN_CFLAGS)

CXXLD = gfortran
#call_trilinos.cpp: call_trilinos.cpp
#	$(NOX_CXX_COMPILER) $(NOX_CXX_FLAGS) $(CXXFLAGS) -c $^
endif TRILINOS

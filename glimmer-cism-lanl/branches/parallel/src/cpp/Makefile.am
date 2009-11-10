if TRILINOS

bin_PROGRAMS  = helloEpetra

variabledir = $(pkgdatadir)
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS := $(wildcard *.mod)
BUILT_SOURCES = 
CLEANFILES = $(BUILT_SOURCES) *.mod

helloEpetra_SOURCES = helloCXX.cpp interface_modCXX.F90 helloF.F90
helloEpetra_LDADD = $(BLAS_LIBS) $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

ALL_SRC=$(helloEpetra_SOURCES)

AM_FCFLAGS = $(CXXFLAGS) $(NETCDF_CPPFLAGS) $(HDF5_CPPFLAGS)
AM_CXXFLAGS = $(CXXFLAGS)

LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) 

AM_CFLAGS=$(CFORTRAN_CFLAGS)

CXXLD = gfortran

endif TRILINOS

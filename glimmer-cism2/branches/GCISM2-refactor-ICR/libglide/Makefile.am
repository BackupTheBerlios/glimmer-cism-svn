lib_LTLIBRARIES = libglide.la
bin_PROGRAMS  =  
check_PROGRAMS = test_setup

variable_DATA = glide_vars.def lithot_vars.def time_vars.def
variabledir = $(pkgdatadir)
dist_variable_DATA = $(variable_DATA)

noinst_LTLIBRARIES = libisostasy.la liblithot.la libvelo.la

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak $(RESTART_INC_INS)
BUILT_SOURCES = glide_io.F90  lithot_io.F90
CLEANFILES = *.$(F90_MODULE_EXTENSION) $(BUILT_SOURCES) f90_dep.mak

libglide_la_SOURCES = glide.F90 glide_setup.F90 glide_types.F90 glide_temp.F90 \
                      glide_thck.F90  \
                      glide_stop.F90 glide_io.F90 lithot_io.F90 glide_nc_custom.F90 \
                      glide_profile.F90 glimmer_routing.F90 \
                      glide_glenflow.F90 \
                      glide_bwat.F90 glide_tempFullSoln.F90 \
                      glide_thckADI.F90 glide_thckCommon.F90
libglide_la_LIBADD = libisostasy.la liblithot.la libvelo.la

libisostasy_la_SOURCES = isostasy.F90 isostasy_el.F90 isostasy_setup.F90 isostasy_types.F90 
liblithot_la_SOURCES = lithot.F90 lithot_setup.F90 lithot_types.F90 lithot1d.F90 
libvelo_la_SOURCES = velo_types.F90 glide_velo.F90

test_setup_SOURCES = test_setup.F90
test_setup_LDADD = libglide.la ../libglimmer/libglimmer.la ../libglimmer/libglimmer-IO.la $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS)

AM_FCFLAGS = $(F90_MODFLAG)../libglimmer $(F90_MODFLAG)../libglimmer-solve -I$(top_srcdir)/libglimmer

#AM_FCFLAGS = $(NETCDF_FCFLAGS)
#LDADD = $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS)
#AM_CFLAGS=$(CFORTRAN_CFLAGS)

#some special rules for automatically generating ncdf src code
%_io.F90:	$(srcdir)/%_vars.def $(top_srcdir)/libglimmer/ncdf_template.F90.in
		$(PYTHON) $(top_srcdir)/scripts/generate_ncvars.py $^


f90_dep.mak:	$(srcdir)/*.F90 $(BUILT_SOURCES)
		$(PYTHON) $(top_srcdir)/scripts/f90_dependencies.py -l -o $@ $^

-include f90_dep.mak

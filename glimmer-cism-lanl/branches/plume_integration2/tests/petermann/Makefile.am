#by default build nc_gen_input and nc_regrid tools

bin_PROGRAMS = nc_gen_input nc_regrid nc_gen_plume nc_find_last_nonnan_time

nc_gen_input_SOURCES = nc_gen_input.F90
nc_gen_input_LDADD = $(NETCDF_LIBS)

nc_regrid_SOURCES = nc_regrid.F90
nc_regrid_LDADD = $(NETCDF_LIBS)

nc_gen_plume_SOURCES = nc_gen_plume.F90
nc_gen_plume_LDADD = $(NETCDF_LIBS)

nc_find_last_nonnan_time_SOURCES = nc_find_last_nonnan_time.F90
nc_find_last_nonnan_time_LDAA = $(NETCDF_LIBS)

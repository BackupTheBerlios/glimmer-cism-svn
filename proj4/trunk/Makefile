libfproj4.a:	fort-proj.o proj4.o
		$(AR) rc $@ $^

test-proj:	test-proj.o libfproj4.a 
		$(F90) $(LDFLAGS) -o $@ $^ -L$(PROJ_PREFIX)/lib -lproj

fort-proj.o:	fort-proj.c
		$(CC) $(CFLAGS) -I$(PROJ_PREFIX)/include -c $< -o $@

proj4.o:	proj4.f90 proj4.inc
		$(F90) $(FFLAGS) -c $< -o $@

test-proj.o:	libfproj4.a

clean:
		rm -f *.o libfproj4.a test-proj

ifdef GLIMMER_PREFIX
install:	libfproj4.a
		install -d $(GLIMMER_PREFIX)/lib
		cp $^ $(GLIMMER_PREFIX)/lib
		install -d $(GLIMMER_PREFIX)/mod
		cp *.mod $(GLIMMER_PREFIX)/mod
		install -d $(GLIMMER_PREFIX)/include
		cp *.inc $(GLIMMER_PREFIX)/include
endif

f90_proj4.tar.gz:	ChangeLog COPYING README cfortran.doc cfortran.h fort-proj.c proj4.f90 proj4.inc test-proj.f90
		tar -cvf - $^ | gzip > $@

.f90.o:
		$(F90) $(FFLAGS) -c $< -o $@

.SUFFIXES:	.c .f90 .o
TRILINOS_PREFIX=/home/hdnong/work/Trilinos_DB/200910_SERIAL
include $(TRILINOS_PREFIX)/include/Makefile.export.NOX
INC_PATH=-I$(TRILINOS_PREFIX)/include
LIB_PATH=-L$(TRILINOS_PREFIX)/lib
FORTRAN_COMPILER=gfortran
FORTRAN_LIBRARIES=-lgfortranbegin -lgfortran

test1 = SI_Test
test2 = simple

.cpp.o:
	$(NOX_CXX_COMPILER) $(NOX_CXX_FLAGS) $(INC_PATH) -I/usr/lib64/openmpi/1.3.2-gcc/include -c $<

sparse_type.o:
	$(FORTRAN_COMPILER) -c sparse_type.F90

$(test2).o: sparse_type.o
	$(FORTRAN_COMPILER) -c $(test2).F90

$(test1).exe: $(test1).o Simple_Interface.o
	$(NOX_CXX_COMPILER) $(NOX_CXX_FLAGS) -o $(test1).exe $^ $(LIB_PATH) $(NOX_LIBRARIES) $(NOX_TPL_LIBRARIES)

$(test2).exe: $(test2).o Simple_Interface.o SI_Drive.o sparse_type.o
	$(NOX_CXX_COMPILER) $(NOX_CXX_FLAGS) -o $(test2).exe $^ $(FORTRAN_LIBRARIES) $(LIB_PATH) $(NOX_LIBRARIES) $(NOX_TPL_LIBRARIES)

clean :
	rm -f *.o *.mod *.exe *~

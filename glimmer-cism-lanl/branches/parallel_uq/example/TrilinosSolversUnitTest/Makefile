TRILINOS_PREFIX=/home/hdnong/work/Trilinos_DB/200910_SERIAL
include $(TRILINOS_PREFIX)/include/Makefile.export.NOX
INC_PATH=-I$(TRILINOS_PREFIX)/include
LIB_PATH=-L$(TRILINOS_PREFIX)/lib
FORTRAN_COMPILER=gfortran
FORTRAN_LIBRARIES=-lgfortranbegin -lgfortran

test = TrilinosSolversUnitTest

.cpp.o:
	$(NOX_CXX_COMPILER) $(NOX_CXX_FLAGS) $(INC_PATH) -I/usr/lib64/openmpi/1.3.2-gcc/include -c $<

sparse_type.o:
	$(FORTRAN_COMPILER) -c sparse_type.F90

sparse_trilinos.o: sparse_type.o
	$(FORTRAN_COMPILER) -c sparse_trilinos.F90

$(test).o: sparse_type.o sparse_trilinos.o
	$(FORTRAN_COMPILER) -c $(test).F90

$(test).exe: $(test).o call_trilinos.o sparse_trilinos.o sparse_type.o
	$(NOX_CXX_COMPILER) $(NOX_CXX_FLAGS) -o $(test).exe $^ $(FORTRAN_LIBRARIES) $(LIB_PATH) $(NOX_LIBRARIES) $(NOX_TPL_LIBRARIES)

clean :
	rm -f *.o *.mod *.exe *~

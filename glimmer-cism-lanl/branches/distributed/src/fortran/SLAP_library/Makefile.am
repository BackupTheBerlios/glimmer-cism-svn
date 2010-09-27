lib_LIBRARIES = libslap.a
F77=ftn -O2 
FC=ftn -O2
F90=ftn -O2
if EXT_BLAS
 BLAS =
else
 BLAS = blas.f
endif
libslap_a_SOURCES = $(BLAS) dcgs.f dlaputil.f dmset.f dmvops.f mach.f xersla.f \
                            dslugm.f dgmres.f dpigmr.f dorth.f dheqr.f isdgmr.f drlcal.f \
                            dhels.f dxlcal.f
EXTRA_libslap_a_SOURCES = blas.f

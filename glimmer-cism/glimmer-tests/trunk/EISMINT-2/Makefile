
runs:	e2.a.nc e2.b.nc e2.c.nc e2.d.nc e2.f.nc e2.g.nc e2.h.nc e2.i.nc e2.j.nc e2.k.nc e2.l.nc
plots:	e2.a.ps e2.b.ps e2.c.ps e2.d.ps e2.f.ps e2.g.ps e2.h.ps e2.i.ps e2.j.ps e2.k.ps e2.l.ps e2.a-comp.ps e2.b-comp.ps e2.c-comp.ps e2.d-comp.ps e2.g-comp.ps e2.h-comp.ps

e2.b.nc:	e2.a.nc
e2.c.nc:	e2.a.nc
e2.d.nc:	e2.a.nc
e2.j.nc:	e2.i.nc
e2.l.nc:	e2.k.nc

e2.a-comp.ps:	e2.a.nc
	plotEISMINT2stats.py -e A $^ $@
e2.b-comp.ps:	e2.b.nc e2.a.nc
	plotEISMINT2stats.py -e B $^ $@
e2.c-comp.ps:	e2.c.nc e2.a.nc
	plotEISMINT2stats.py -e C $^ $@
e2.d-comp.ps:	e2.d.nc e2.a.nc
	plotEISMINT2stats.py -e D $^ $@
e2.g-comp.ps:	e2.g.nc
	plotEISMINT2stats.py -e G $^ $@
e2.h-comp.ps:	e2.h.nc
	plotEISMINT2stats.py -e H $^ $@

.config.nc:
	glide_launch.py $<

.nc.ps:
	plotEISMINT.py -T -1  -p ../PyCF/eismint.prof $< $@

.SUFFIXES: .config .nc .ps


runs:	e2.a.nc e2.b.nc e2.c.nc e2.d.nc e2.f.nc e2.g.nc e2.h.nc e2.i.config e2.j.nc e2.k.nc e2.l.nc

e2.b.nc:	e2.a.nc
e2.c.nc:	e2.a.nc
e2.d.nc:	e2.a.nc
e2.j.nc:	e2.i.nc
e2.l.nc:	e2.k.nc

.config.nc:
	glide_launch.py $<

.SUFFIXES: .config .nc

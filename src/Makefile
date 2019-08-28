FC = mpif90
FFLAGS = -O2 -ffast-math -w

FILES = global.f90 comm.f90 rndmod.f90 sysdef.f90 config.f90 nnlist.f90 \
	measure.f90 forcem.f90 move.f90 main.f90

all: clean

mdsim: $(FILES)
	$(FC) $(FFLAGS) -o $@ $(FILES)

clean: mdsim
	rm *.mod

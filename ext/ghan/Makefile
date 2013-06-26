FORT=gfortran
FFLAGS=-ffixed-form

all: subgrid_mod.o ghan_interface.o
subgrid_mod.o:
	$(FORT) -c subgrid_mod.f90 -o subgrid_mod.o -ffixed-form
ghan_interface.o: subgrid_mod.o
	$(FORT) -c ghan_interface.f90 -o ghan_interface.o
	f2py -c subgrid_mod.o -m ghan_interface ghan_interface.f90
clean:
	rm -f *.o *.so *.mod fort.100
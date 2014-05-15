FORT=gfortran
FFLAGS=-ffixed-form
CC=gcc
CCFLAGS=

all: subgrid_mod.o ghan_interface.o
	#python ghan_test.py
subgrid_mod.o:
	$(FORT) -c subgrid_mod.f90 -o subgrid_mod.o $(FFLAGS)
ghan_interface.o: subgrid_mod.o
	$(FORT) -c ghan_interface.f90 -o ghan_interface.o
	CC="$(CC) $(CCFLAGS)" f2py -c subgrid_mod.o -m ghan_interface ghan_interface.f90
clean:
	rm -f *.pyc *.o *.so *.mod fort.100
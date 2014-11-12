all:
	gfortran -c constants.f90 
	gfortran -c util.f90
	gfortran -c thermo.f90 
	gfortran -c aerosol.f90
	gfortran -c dvode_f90_m.f90 
	gfortran -c model.f90 
	gfortran -c main.f90 
	gfortran *.o -o main 
clean:
	rm -rf *.mod *.o
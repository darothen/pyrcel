# FC = the compiler to use 
FC=gfortran

# Compiler options 
FFLAGS=-O 

# List libraries used by the program here 
LIBS=


# Suffix-rules:  Begin by throwing away all old suffix- 
# rules, and then create new ones for compiling  
# *.f90-files. 
.SUFFIXES: 
.SUFFIXES: .f90 .o

.f90.o: 
	$(FC) -c $(FFLAGS) $<

# Include the dependency-list created by makedepf90 below 
include .depend

# target 'clean' for deleting object- *.mod- and other  
# unwanted files 
clean: 
	cp dvode_f90_m.mod dvode_f90_m.mod.bak
	cp dvode_f90_m.o dvode_f90_m.o.bak
	rm -f *.o *.mod core *.out
	cp dvode_f90_m.o.bak dvode_f90_m.o
	cp dvode_f90_m.mod.bak dvode_f90_m.mod
	rm *.bak

# Create a dependency list using makedepf90.  All files  
# that needs to be compiled to build the program,  
# i.e all source files except include files, should  
# be given on the command line to makedepf90.   
# 
# The argument to the '-o' option will be the name of the 
# resulting program when running 'make', in this case  
# 'foobar' 
depend .depend: 
	makedepf90 -o main *.f90 > .depend
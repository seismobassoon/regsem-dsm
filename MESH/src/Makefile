# Makefile

FC = ifort   # Fortran 90 compiler
OPT = -O3 #-ax -ipo -ftz
MMOD = -module
MODULEDIR = ./Elliptic/
F90FLAGS = $(MMOD) $(MODULEDIR)

OBJ = $(patsubst %.f90,%.o, $(wildcard *.f90) )
OBJ_MOD = $(patsubst %.f90,%.o, $(wildcard $(MODULEDIR)*.f90))

all : Mesher.x

Elliptic/nrutil.o :	Elliptic/nrtype.o
Elliptic/module_spline.o :	Elliptic/nrutil.o Elliptic/nrtype.o
Elliptic/module_ellipticity.o :	Elliptic/module_spline.o Elliptic/funaro.o

Mesher.x : make_module $(OBJ) 
	$(FC) $(F90FLAGS) $(OPT) -o ../bin/$@ $(OBJ) $(OBJ_MOD) -L. Metis/libmetis.a

make_module :   $(OBJ_MOD)

Elliptic/%.o : Elliptic/%.f90
	$(FC) $(MMOD) Elliptic $(OPT) -c $<  -o $@

%.o : %.f90 
	$(FC) $(F90FLAGS) $(OPT) -c $<  -o $@ 


clean: clean_modules cleanx
	rm -f *.o

clean_modules:
	rm -f  $(MODULEDIR)*.o $(MODULEDIR)*.mod

cleanx:
	rm -f ../bin/*.x

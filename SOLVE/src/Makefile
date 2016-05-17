# Makefile

FC = 		mpif90   # MPI - Fortran 90 compiler
CC = 		gcc   # C compiler
OPT = 		-O3 #-axS -ipo -ftz
MMOD = 		-module
MODULEDIR = 	./Modules/
F90FLAGS =	$(MMOD) $(MODULEDIR)

# patsubst = path substitution, wildcard = tout
# si le repertoire contient Fn1.c Fn2.c Fn3.f90 et Fn4.f90 alors OBJ = "Fn1.o Fn2.o Fn3.o Fn4.o"
OBJ = 		$(patsubst %.c,%.o, $(wildcard *.c))  $(patsubst %.f90,%.o, $(wildcard *.f90) ) 
OBJ2 = 		$(patsubst %.f,%.o, $(wildcard BLAS/*.f)) 
OBJ_MOD = 	$(patsubst %.c,%.o, $(wildcard $(MODULEDIR)*.c)) $(patsubst %.f90,%.o, $(wildcard $(MODULEDIR)*.f90))

all : spec.x

Modules/Domain.o :	$(MODULEDIR)Element.o $(MODULEDIR)Face.o $(MODULEDIR)Source.o $(MODULEDIR)Receiver.o \
                        $(MODULEDIR)Vertex.o $(MODULEDIR)TimeParam.o $(MODULEDIR)Subdomain.o $(MODULEDIR)Edge.o \
                        $(MODULEDIR)Adjoint.o $(MODULEDIR)Comm.o
Modules/Element.o :	$(MODULEDIR)Simu.o
Modules/Face.o :	$(MODULEDIR)Simu.o
Modules/Edge.o :	$(MODULEDIR)Simu.o
Modules/Vertex.o :	$(MODULEDIR)Simu.o
Modules/Source.o :	$(MODULEDIR)angles.o
Modules/nrutil.o :	$(MODULEDIR)nrtype.o
Modules/module_spline.o :	$(MODULEDIR)nrutil.o $(MODULEDIR)nrtype.o
Modules/init_cond.o :	$(MODULEDIR)module_spline.o
Modules/module_ellipticity.o :	$(MODULEDIR)module_spline.o $(MODULEDIR)funaro.o
Modules/module_A3d.o :	$(MODULEDIR)def_gparam.o $(MODULEDIR)earth_modele.o \
                        $(MODULEDIR)module_spline.o $(MODULEDIR)spl_A3d.o


spec.x : make_module $(OBJ) make_BLAS
	$(FC) $(F90FLAGS) $(OPT) -o ../bin/$@ $(OBJ) $(OBJ_MOD) $(OBJ2)

make_module : 	$(OBJ_MOD)
make_BLAS :	$(OBJ2)

BLAS/%.o : BLAS/%.f
	$(FC) $(OPT) -c $< -o $@

Modules/%.o : Modules/%.f90 
	$(FC) $(MMOD) Modules $(OPT) -c $< -o $@

Modules/%.o : Modules/%.c
	$(CC) -O2 -c $< -o $@

%.o : %.f90 
	$(FC) $(F90FLAGS) $(OPT) -c $< -o $@

%.o : %.c
	$(CC) -O2 -c $<


clean: clean_modules cleanx clean_BLAS
	rm -f *.o

clean_modules:
	rm -rf  $(MODULEDIR)/*.o $(MODULEDIR)/*.mod

cleanx:
	rm -f ../bin/*.x

clean_BLAS:
	rm -f BLAS/*.o

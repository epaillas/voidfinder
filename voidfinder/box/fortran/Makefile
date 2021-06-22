# Compiler choice
FC = gfortran
MPIFC = gfortran

MLIBS =

#################

EXES = grow_spheres.exe 

all: $(EXES)

grow_spheres.exe: grow_spheres.f90
	$(FC) -fopenmp procedures.f90 grow_spheres.f90 -o grow_spheres.exe $(MLIBS)


install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod

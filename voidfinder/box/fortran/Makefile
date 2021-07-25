# Compiler choice
FC = gfortran
MPIFC = gfortran

MLIBS =

#################

EXES = grow_spheres.exe recentring.exe overlapping.exe

all: $(EXES)

grow_spheres.exe: grow_spheres.f90
	$(FC) -fopenmp procedures.f90 grow_spheres.f90 -o grow_spheres.exe $(MLIBS)

recentring.exe: recentring.f90
	$(FC) -fopenmp procedures.f90 recentring.f90 -o recentring.exe $(MLIBS)

overlapping.exe: overlapping.f90
	$(FC) overlapping.f90 -o overlapping.exe $(MLIBS)

install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod

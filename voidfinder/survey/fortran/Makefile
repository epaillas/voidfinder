# Compiler choice
FC = gfortran
MPIFC = gfortran
IFC = ifort

MLIBS =

#################

EXES = tpcf.exe tpccf_monopole.exe tpccf_monopole_2d.exe tpccf_rmu.exe

all: $(EXES)

tpcf.exe: tpcf.f90
	$(FC) -fopenmp procedures.f90 tpcf.f90 -o tpcf.exe $(MLIBS)

tpccf_monopole.exe: tpccf_monopole.f90
	$(FC) -fopenmp procedures.f90 tpccf_monopole.f90 -o tpccf_monopole.exe $(MLIBS)

tpccf_monopole_2d.exe: tpccf_monopole_2d.f90
	$(FC) -fopenmp procedures.f90 tpccf_monopole_2d.f90 -o tpccf_monopole_2d.exe $(MLIBS)

tpccf_rmu.exe: tpccf_rmu.f90
	$(FC) -fopenmp procedures.f90 tpccf_rmu.f90 -o tpccf_rmu.exe $(MLIBS)


install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod

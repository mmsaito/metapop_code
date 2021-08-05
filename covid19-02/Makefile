#FC  = mpif90
FC = gfortran
CC  = gcc
#LD  = mpif90

#MPIFC = mpiifort
#MPICC = mpiicc
#MPILD = mpiifort

#CC = gcc
#FC = $(MPIFC)
#LD = $(MPILD)

#OPT = -c -Wall -O2 -ftrace=full -CB -g
#OPT = -cpp -c -CB -O2 -nogen-interfaces -warn -traceback -g
#OPT = -cpp -c -CB -O3 -nogen-interfaces -traceback -g 
OPT = -cpp -c  -O3  -g  -I./ -I../msmpi-x86-patched -fno-range-check
MPILIB = c:/windows/system32/msmpi.dll

all: metapop_da  metapop_pmcmc

metapop_pmcmc.o: ../metapop_pmcmc.f90 modify_par.f90

c_misc.o: ../c_misc.c
	$(CC) $(OPT) -c $<

OBJ_EXE = c_misc.o misc.o rand.o model.o pf.o rubella_03_lik2_step7.o
rubella_03_lik2_step7: $(OBJ_EXE)
	$(FC) $(OBJ_EXE) -o $@

OBJ_EXEb = c_misc.o misc.o rand.o model.o pf.o metapop_da.o
metapop_da: $(MPILIB) $(OBJ_EXEb)
	$(FC) $(OBJ_EXEb) $(MPILIB) -o $@

OBJ_EXE2 = c_misc.o misc.o rand.o model.o pf.o rubella_04_lik2_step7.o
rubella_04_lik2_step7: $(OBJ_EXE2)
	$(FC) $(OBJ_EXE2) -o $@

# single jump of beta
OBJ_EXE3 = c_misc.o misc.o rand.o model.o pf.o metapop_pmcmc.o
metapop_pmcmc: $(OBJ_EXE3)
	$(FC) $(OBJ_EXE3) $(MPILIB) -o $@

# multi jump of beta
OBJ_EXE4 = c_misc.o misc.o rand.o model.o pf.o rubella_06_lik2_step7.o
rubella_06_lik2_step7: $(OBJ_EXE4)
	$(FC) $(OBJ_EXE4) -o $@

OBJ_EXE999 = c_misc.o misc.o rand.o gomi.o
a.out: $(OBJ_EXE999)
	$(FC) $(OBJ_EXE999) -o $@

.SUFFIXES:
.SUFFIXES: .f .f90 .c .obj .o .exe

.f90.o:
	$(FC) $(OPT) $< -o $@

#.c.o:
#	$(MPICC) $(OPT) $< -o $@

clean:
	rm *.o *.mod *.exe


OBJ_EXE999 = c_misc.o misc.o rand.o mpi_sort_test2.o
a.out: $(OBJ_EXE999)
	$(FC) $(OBJ_EXE999) -o $@

.SUFFIXES:
.SUFFIXES: .f .f90 .c .obj .o .exe

vpath %.c ../
vpath %.f90 ../

.f90.o:
	$(FC) $(OPT) $< -o $@

#.c.o:
#	$(MPICC) $(OPT) $< -o $@

pwd:
	gfortran -E -DX=\"`pwd`/const\"  hoge.F90

clean:
	rm *.o *.mod *.exe metapop_da  metapop_pmcmc

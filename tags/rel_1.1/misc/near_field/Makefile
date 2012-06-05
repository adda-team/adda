# Makefile for Near Field tools

# ---- Intel Fortran Compiler Linux
#FC     = ifort
#FFLAGS = -warn all -O3
#MPIFC  = mpif90
#LIBS   =
#MPILIBS= -L/opt/hpmpi/lib/linux_64 -lmpi

# ----- HP-UX
#FC      = f90
#MPIFC   = mpif90
#FFLAGS  = +O3
#LIBS    =
#MPILIBS =  -L/opt/mpi/lib/hpux64/ -lmpio

# ---- GNU Fortran Compiler
FC     = gfortran
FFLAGS = -O3 -ffast-math -funroll-loops
MPIFC  = $(FC)
LIBS   =
MPILIBS= -lfmpich2g
# You might need something like the following on Windows (with MinGW)
MPILIBS+=  -L"C:/Program Files/MPICH2/lib" -I"C:/Program Files/MPICH2/include"

ALLPROGS = gendip genplane genline nearfield checknfsymm


all: $(ALLPROGS)
	@echo 'everything done!'

xfuty.o : xfuty.f Makefile
	$(FC)  -c $(FFLAGS)  xfuty.f

nfsubs.o : nfsubs.f Makefile
	$(FC)  -c $(FFLAGS)  nfsubs.f

gendip : gendip.f   Makefile
	$(FC)   $(FFLAGS) -o $@  gendip.f

genline : genline.f Makefile
	$(FC)   $(FFLAGS) -o $@  genline.f

genplane : genplane.f Makefile
	$(FC)   $(FFLAGS) -o $@  genplane.f

nearfield : nearfield.f  xfuty.o nfsubs.o Makefile
	$(MPIFC)  $(FFLAGS) -o $@  nearfield.f $(MPILIBS) xfuty.o nfsubs.o

checknfsymm : checknfsymm.f  nfsubs.o Makefile
	$(FC)  $(FFLAGS) -o $@  checknfsymm.f nfsubs.o

clean:
	rm -f *.o $(ALLPROGS) $(addsuffix .exe,$(ALLPROGS))

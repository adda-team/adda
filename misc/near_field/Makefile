# Makefile for Near Field tools

# ---- Intel Fortran Compiler Linux
#FC     = ifort
#FFLAGS = -warn all -O3
#MPIFC  = mpif90
#MPILIBS= -L/opt/hpmpi/lib/linux_64 -lmpi

# ----- HP-UX
#FC      = f90
#MPIFC   = mpif90
#FFLAGS  = +O3
#MPILIBS =  -L/opt/mpi/lib/hpux64/ -lmpio

# ---- GNU Fortran Compiler
FC     = gfortran
# -ffast-math causes IEEE_INVALID_FLAG in checknfsymm, this is probably OK but safer is to skip this optimization
FFLAGS = -O3 -funroll-loops $(EXTRA_FLAGS)
# depending on your MPI implementation you may need one of the following (either special compiler wrapper, which takes
# care of libraries), or standard compiler with explicit specification of MPI Fortran libraries
MPIFC = mpifort
#MPIFC  = $(FC)
#MPILIBS= -lfmpich2g
# You might need something like the following on Windows (with MinGW)
#MPILIBS += -I"C:/Program Files/MPICH2/include"
#MPILIBS += -L"C:/Program Files/MPICH2/lib" -I"C:/Program Files/MPICH2/include"

ALLPROGS = gendip genplane genline nearfield checknfsymm

srcdir = .
vpath %.f $(srcdir)/
vpath Makefile $(srcdir)/

#=======================================================================================================================

all: $(ALLPROGS)
	@echo 'everything done!'

xfuty.o : xfuty.f Makefile
	$(FC)  -c $(FFLAGS) $<

nfsubs.o : nfsubs.f Makefile
	$(FC)  -c $(FFLAGS) $<

gendip : gendip.f   Makefile
	$(FC)   $(FFLAGS) -o $@ $<

genline : genline.f Makefile
	$(FC)   $(FFLAGS) -o $@ $<

genplane : genplane.f Makefile
	$(FC)   $(FFLAGS) -o $@ $<

nearfield : nearfield.f  xfuty.o nfsubs.o Makefile
	$(MPIFC)  $(FFLAGS) -o $@ $< $(MPILIBS) xfuty.o nfsubs.o

checknfsymm : checknfsymm.f  nfsubs.o Makefile
	$(FC)  $(FFLAGS) -o $@ $< nfsubs.o

clean:
	rm -f *.o $(ALLPROGS) $(addsuffix .exe,$(ALLPROGS))

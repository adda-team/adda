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
FFLAGS = -O3 -ffast-math -funroll-loops $(EXTRA_FLAGS)
MPIFC  = $(FC)
LIBS   =
MPILIBS= -lfmpich2g
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

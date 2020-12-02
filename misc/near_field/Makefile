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

# Microsoft MPI on Windows is detected automatically. If you have it installed, but want to use another 
# implementation, delete the last symbol on the following line (then this case will not be satisfied).
# Alternatively, if this logic is not supported by your old make, remove the case completely.
ifdef MSMPI_INC
  MPIFC    = $(FC)
  MPILIBS  = -lmsmpifec -lmsmpi -I'$(MSMPI_INC)' -I'$(MSMPI_INC)\x64' -L'$(MSMPI_LIB64)' -fno-range-check
  # the following is due to the bug in MS_MPI - https://github.com/microsoft/Microsoft-MPI/issues/7
  # this workaround was proposed by Maik Riechert; still warnings remain
  CFG_STUB = cfg_stub.o
  CC       = gcc
else
  # depending on your MPI implementation you may need one of the following (either special compiler wrapper, which 
  # takes care of libraries), or standard compiler with explicit specification of MPI Fortran libraries and maybe
  # pathes to libraries and headers
  MPIFC = mpifort
  #MPIFC  = $(FC)
  #MPILIBS= -lfmpich2g
  #MPILIBS += -L"C:/Program Files/MPICH2/lib"
  #MPILIBS += -I"C:/Program Files/MPICH2/include"  
  # Depending on a particular MPI installation one may need to manually specify paths to MPI headers and libraries.
endif

ALLPROGS = gendip genplane genline nearfield checknfsymm

srcdir = .
vpath %.f $(srcdir)/
vpath Makefile $(srcdir)/

#=======================================================================================================================

all: $(ALLPROGS)
	@echo 'everything done!'

xfuty.o : xfuty.f Makefile
	$(FC) -c $(FFLAGS) $<

nfsubs.o : nfsubs.f Makefile
	$(FC) -c $(FFLAGS) $<

gendip : gendip.f   Makefile
	$(FC) $(FFLAGS) -o $@ $<

genline : genline.f Makefile
	$(FC) $(FFLAGS) -o $@ $<

genplane : genplane.f Makefile
	$(FC) $(FFLAGS) -o $@ $<

cfg_stub.o: cfg_stub.c Makefile
	$(CC) -c $<

nearfield : nearfield.f xfuty.o nfsubs.o Makefile $(CFG_STUB)
	$(MPIFC) $(FFLAGS) -o $@ $< $(MPILIBS) xfuty.o nfsubs.o $(CFG_STUB)

checknfsymm : checknfsymm.f nfsubs.o Makefile
	$(FC) $(FFLAGS) -o $@ $< nfsubs.o

clean:
	rm -f *.o $(ALLPROGS) $(addsuffix .exe,$(ALLPROGS))

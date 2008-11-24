# Main makefile for ADDA package
# Actual compiling goes in make_seq and make_mpi for
# sequential and parallel (mpi) version respectively
# AUTH: Maxim Yurkin

#==========================================================
# Optional flags; can be modified by user
# information can be found in corresponding headers  

# for faster recompiling turn on in debug.h
#CFLAGS += -DDEBUG             

# info in fft.h; change here 
#CFLAGS += -DFFT_TEMPERTON

# for faster recompiling turn on in timing.h
#CFLAGS += -DPRECISE_TIMING

# for faster recompiling turn on in io.h
#CFLAGS += -DNOT_USE_LOCK

# Choose either gcc or icc;
# other should be implemented in sub-makefiles
COMPILER = gnu
#COMPILER = intel
#COMPILER = other

#==========================================================
# Unconditional variables 

CSOURCE  = ADDAmain.c  calculator.c make_particle.c CalculateE.c \
           GenerateB.c debug.c timing.c linalg.c \
           iterative.c memory.c comm.c matvec.c param.c \
           fft.c crosssec.c Romberg.c io.c prec_time.c vars.c
# One file is added below
FSOURCE  =  

COBJECTS = $(CSOURCE:.c=.o)
CDEPEND  = $(CSOURCE:.c=.d)
FOBJECTS = $(FSOURCE:.f=.o)
PROG     = adda
SHELL    = /bin/sh
LDLIBS   = -lm
DEPFLAG  = -MM
MFILES   = Makefile
# Fortran sources generate a lot of warnings
FWARN    = -w

#==========================================================
# Conditional variables, that depend on the
# values of optional flags

# adds link to fftw3 library when above is commented 
ifeq ($(findstring $(CFLAGS),-DFFT_TEMPERTON),)
  LDLIBS += -lfftw3
else
  FSOURCE += cfft99D.f
endif


#==========================================================
# main target part

.EXPORT_ALL_VARIABLES:
.PHONY: seq mpi clean

seq: 
	make -f make_seq

mpi:
	make -f make_mpi

clean:
	rm -f *.o
	rm -f *.d
	rm -f $(PROG)

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

# Choose compiler (gcc, icc, or other);
# other should be implemented below and in sub-makefiles;
# gcc tested for version 3.2.3 and higher;
# higher versions of intel compilers should work with intel9.0
COMPILER = gnu
#COMPILER = intel9.0
#COMPILER = intel8.1
#COMPILER = other

#If uncommented, all warning are suppressed
RELEASE = on

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
DEPFLAG  = -M
MFILES   = Makefile
# Fortran sources generate a lot of warnings
FWARN    = -w

#==========================================================
# Conditional variables, that depend on the
# values of optional flags

# adds link to fftw3 library when above is commented 
ifeq ($(findstring -DFFT_TEMPERTON,$(CFLAGS)),)
  LDLIBS += -lfftw3
else
  FSOURCE += cfft99D.f
endif
# if 'release' turn off warningns
ifdef RELEASE
  CWARN = -w
  LWARN = -w
endif 
# compiler warnings and optimization flags
FOPT    = $(COPT)
ifeq ($(COMPILER),gnu)
# You can add option '-march=...' to COPT and FOPT in this section
  COPT  = -O3 -ffast-math -funroll-loops
  FOPT  = -O
  ifndef RELEASE
    CWARN = -pedantic -Wall -W -Winline -Wpointer-arith -Wcast-qual \
            -Wwrite-strings -Wstrict-prototypes -Wno-uninitialized \
            -Wno-unknown-pragmas -Wno-comment -Wno-unused-parameter
  endif
endif
ifeq ($(COMPILER),intel9.0)
  COPT  = -fast -fp-model fast
  FOPT  = -fast
  ifndef RELEASE
    CWARN = -Wall -Wno-uninitialized -Wno-comment \
            -wd279,869,981,1418,1419,1572
  endif
endif  
ifeq ($(COMPILER),intel8.1)
  COPT  = -fast
  ifndef RELEASE
    CWARN = -Wall -wd9,279,810,869,981,1418,1419,1572
  endif
endif  
ifeq ($(COMPILER),other)
endif
# Finalize option flags (almost)
CFLAGS += $(COPT) $(CWARN)
FFLAGS += $(FOPT) $(FWARN)
LFLAGS += $(COPT) $(LWARN) 

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

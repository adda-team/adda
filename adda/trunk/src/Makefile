# Main makefile for ADDA package
# Actual compiling goes in make_seq and make_mpi for
# sequential and parallel (mpi) version respectively
# AUTH: Maxim Yurkin

#===============================================================================
# !!! Start of control section. Flags and options here are designed to be
# modified by user to choose particular options for compilation. However, the
# default set of options may work out-of-box on some systems. If compiling with
# MPI support please also look at control section in file 'make_mpi'.
#===============================================================================

# --Optional flags-- Uncomment those that you find appropriate;
# information can be found in corresponding headers

# for faster recompiling turn on in debug.h
#CFLAGS += -DDEBUG

# info in fft.h; change here
#CFLAGS += -DFFT_TEMPERTON

# for faster recompiling turn on in prec_timing.h
#CFLAGS += -DPRECISE_TIMING

# for faster recompiling turn on in io.h
#CFLAGS += -DNOT_USE_LOCK
#CFLAGS += -DONLY_LOCKFILE

# --Compilers-- Choose (uncomment) one of the lines below

# tested for gcc 3.2.3 and higher;
COMPILER = gnu

# uses -fast option, which is general, but includes static linking, which may
# cause problems in some cases
#COMPILER = intel9.x

# the same as intel9.x but do not use static linking; it may also not use all
# possible optimizations, depending on a particular processor
#COMPILER = intel9.x_ns

# this is almost the same as intel9.x with few changes for old version, obsolete
#COMPILER = intel8.1

# tested on Compaq C V6.5-303 (dtk)
#COMPILER = compaq

# All options should be implemented below and in sub-makefiles
#COMPILER = other

# --Warnings-- If uncommented, all warning are suppressed
RELEASE = on

# --FFTW3 paths--
# Specify path to headers and libraries of FFTW3. Some systems
# do not need them at all, some specify special global variables
# (first 2 lines), on some - FFTW3 is installed under user account
# (next 2 lines). Under Windows it may be required to specify manually
# relative or absolute paths (last 2 lines).
#FFTW3_INC_PATH = "$(FFTW_INC)"
#FFTW3_LIB_PATH = "$(FFTW_LIB)"
#FFTW3_INC_PATH = "$(HOME)/include"
#FFTW3_LIB_PATH = "$(HOME)/lib"
#FFTW3_INC_PATH = "./../lib"
#FFTW3_LIB_PATH = "./../lib"

#===============================================================================
# !!! End of control section. Everything below is not designed to be modified by
# user. However, advanced users may wish to modify some compilers  flags below,
# especially when using 'other' compiler.

# Unconditional variables
#===============================================================================


CSOURCE  = ADDAmain.c  calculator.c make_particle.c CalculateE.c \
           GenerateB.c debug.c timing.c linalg.c \
           iterative.c memory.c comm.c matvec.c param.c \
           fft.c crosssec.c Romberg.c io.c prec_time.c vars.c \
           mt19937ar.c sinint.c
# One file is added below
FSOURCE  =

COBJECTS = $(CSOURCE:.c=.o)
CDEPEND  = $(CSOURCE:.c=.d)
FOBJECTS = $(FSOURCE:.f=.o)
PROGSEQ  = adda
PROGMPI  = adda_mpi
LASTSEQ  = .lastseq
LASTMPI  = .lastmpi
SHELL    = /bin/sh
LDLIBS   = -lm
DEPFLAG  = -M
MFILES   = Makefile
# Fortran sources generate a lot of warnings
FWARN    = -w

#===============================================================================
# Conditional variables that depend on the values of optional flags
#===============================================================================

# adds link to fftw3 library and modifies paths when above is commented
ifeq ($(findstring -DFFT_TEMPERTON,$(CFLAGS)),)
  LDLIBS += -lfftw3
  ifneq ($(FFTW3_INC_PATH),)
    CFLAGS += -I$(FFTW3_INC_PATH)
  endif
  ifneq ($(FFTW3_LIB_PATH),)
    LFLAGS += -L$(FFTW3_LIB_PATH)
  endif
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
  CC    = gcc
  CF    = g77
  COPT  = -std=c99 -O3 -ffast-math -funroll-loops
  FOPT  = -O
  ifndef RELEASE
    CWARN = -pedantic -Wall -W -Winline -Wpointer-arith -Wcast-qual \
            -Wwrite-strings -Wstrict-prototypes -Wno-uninitialized \
            -Wno-unknown-pragmas -Wno-unused-parameter -Wno-comment
  endif
endif
ifeq ($(COMPILER),intel9.x)
  CC    = icc
  CF    = ifort
  COPT  = -fast -fp-model fast
  FOPT  = -fast
  ifndef RELEASE
    CWARN = -Wall -Wno-uninitialized -Wno-comment \
            -wd279,810,869,981,1418,1419,1572
  endif
endif
ifeq ($(COMPILER),intel9.x_ns)
  CC    = icc
  CF    = ifort
  COPT  = -fp-model fast -O3 -ipo -no-prec-div -xP
  FOPT  = -O3 -ipo -no-prec-div -xP
  ifndef RELEASE
    CWARN = -Wall -Wno-uninitialized -Wno-comment \
            -wd279,810,869,981,1418,1419,1572
  endif
endif
ifeq ($(COMPILER),intel8.1)
  CC    = icc
  CF    = ifort
  COPT  = -fast
  ifndef RELEASE
    CWARN = -Wall -wd9,279,810,869,981,1418,1419,1572
  endif
endif
ifeq ($(COMPILER),compaq)
  # You can add option '-arch host' to COPT and FOPT in this section
  CC    = cc
  CF    = f77
  COPT  = -fast
  ifndef RELEASE
    CWARN = -w0 -msg_disable nestedcomment,unknownpragma,unreachcode
  endif
endif
ifeq ($(COMPILER),other)
endif
# Finalize option flags (almost)
CFLAGS += $(COPT) $(CWARN)
FFLAGS += $(FOPT) $(FWARN)
LFLAGS += $(COPT) $(LWARN)

#===============================================================================
# Main action part
#===============================================================================

.EXPORT_ALL_VARIABLES:
.PHONY: seq mpi clean

seq:
	rm -f $(LASTMPI)
	$(MAKE) -f make_seq

mpi:
	rm -f $(LASTSEQ)
	$(MAKE) -f make_mpi

clean:
	rm -f *.o *.d $(PROGSEQ) $(PROGMPI) $(LASTSEQ) $(LASTMPI)

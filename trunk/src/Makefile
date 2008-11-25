# Main makefile for ADDA package
# Actual compiling goes in make_seq and make_mpi for
# sequential and parallel (mpi) version respectively
# $Author$
# $Date::                            $
#
# Copyright (C) 2006-2008 University of Amsterdam
# This file is part of ADDA.
#
# ADDA is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# ADDA is distributed in the hope that it will be useful,but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# ADDA. If not, see <http://www.gnu.org/licenses/>.

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

# --Compilers-- Choose one of the following
# gnu - tested for gcc 3.2.3 - 3.4.6
# intel - tested on icc 9.0 and 9.1
# compaq - tested on Compaq C V6.5-303 (dtk)
# ibm - tested on xlc 8.0
# other -
COMPILER = gnu

# Compilation mode (bits): 32 or 64
# Although all options in this Makefile are adjusted automatically, different
# modes may require different configurations (environments) of the system,
# especially if MPI mode is used.
BITS = 32

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
#FTW3_LIB_PATH = "./../lib"

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
DFFLAG   = -MF
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
# compiler warnings and optimization flags;
# when release, warning are turned off below
FOPT    = $(COPT)
ifeq ($(COMPILER),gnu)
  # You can add option '-march=...' or 'mcpu=...' to COPT and FOPT in this section
  CC    = gcc
  CF    = g77
  COPT  = -m$(BITS) -std=c99 -O3 -ffast-math -funroll-loops
  FOPT  = -O
  CWARN = -Wall -W -Wpointer-arith -Wcast-qual \
          -Wwrite-strings -Wstrict-prototypes -Wno-uninitialized \
          -Wno-unknown-pragmas -Wno-comment -Wno-unused-parameter
endif
ifeq ($(COMPILER),intel)
  CC    = icc
  CF    = ifort
  COPT  = -fp-model fast -O3 -ipo -no-prec-div -xP
  FOPT  = -O3 -ipo -no-prec-div -xP
  CWARN = -Wall -Wno-uninitialized -Wno-comment \
            -wd279,810,869,981,1418,1419,1572
endif
ifeq ($(COMPILER),compaq)
  # This compiler was not tested since 2007. In particular, warning options may
  # not fit exactly the C99 standard, to which the code was transferred.
  # Its support for 64 bit compilations is also undefined.
  # If you happen to use this compiler, please report results to the authors.
  #
  # You can add option '-arch host' to COPT and FOPT in this section
  CC    = cc
  CF    = f77
  COPT  = -fast
  CWARN = -w0 -msg_disable nestedcomment,unknownpragma,unreachcode
endif
ifeq ($(COMPILER),ibm)
  # -O5 implies "-arch=auto", which tunes compilation exclusively for the host
  # machine. However, it will not work in some configurations.
  # Then use '-O3 -qarch=... -qtune=...' instead
  CC      = xlc
  CF      = xlf
  COPT    = -O3 -qcache=auto -q$(BITS) -qlanglvl=extc99
  # DFREDIRECT = 2>nul -qipa=level=2 -qhot
  DEPFLAG = -qmakedep=gcc -qsyntaxonly
  CWARN   = -qsuppress=1506-224:1506-342:1500-036
endif
ifeq ($(COMPILER),other)
endif
# if 'release' turn off warningns
ifdef RELEASE
  CWARN = -w
  LWARN = -w
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

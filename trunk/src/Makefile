# Main makefile for ADDA package
# Requires GNU make to execute. Actual compiling goes in folders 'seq' and 'mpi' for sequential and
# parallel (MPI) version respectively
# $Author$
# $Date::                            $
#
# Copyright (C) 2006-2008 University of Amsterdam
# Copyright (C) 2009,2010 Institute of Chemical Kinetics and Combustion & University of Amsterdam
# This file is part of ADDA.
#
# ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
# the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along with ADDA. If not, see
# <http://www.gnu.org/licenses/>.

#===================================================================================================
# !!! This makefile accepts environmental variables CFLAGS, FFLAGS, LDFLAGS, and LDLIBS. First two
# are passed to C and Fortran compiler respectively, and last two are passed to linker.
#===================================================================================================

#===================================================================================================
# Fixed definitions for non-compilation targets 
#===================================================================================================

SHELL    = /bin/bash
# Separate directories to store Makefiles and .o and .d files for sequential and MPI executables
SEQ      = seq
MPI      = mpi
# Names of executables
PROGSEQ  = adda
PROGMPI  = adda_mpi
# Targets implying compilation (empty one is always included)
NONTRIVIAL = all seq mpi

#===================================================================================================
# !!! Start of control section. Flags and options here can be modified by user. However, the default
# set of options may work out-of-box on some systems.
# First section contains general variables used by non-compilation targets (clean, install, etc.) 
#===================================================================================================

# nothing here yet

#===================================================================================================
# Below are flags and options specific for compilation. If compiling with MPI support please also 
# look at control section in file 'mpi/Makefile'.

# Everything from here down to the main action part of this makefile is skipped if targets specified 
# in the command line to this makefile do not imply compilation.
ifneq ($(strip $(if $(MAKECMDGOALS), \
        $(if $(strip $(foreach var,$(NONTRIVIAL),$(findstring $(var),$(MAKECMDGOALS)))),1,),1)),)
#===================================================================================================

# --Optional flags-- Uncomment those that you find appropriate; information can be found in
# corresponding headers

# debug.h
#CFLAGS += -DDEBUG

# fft.h
#CFLAGS += -DFFT_TEMPERTON

# prec_timing.h
#CFLAGS += -DPRECISE_TIMING

# io.h
#CFLAGS += -DNOT_USE_LOCK
#CFLAGS += -DONLY_LOCKFILE

# This should be uncommented, if compiling of Fortran sources cause problems. However, currently a
# number of ADDA (optional) features rely on Fortran sources, in particular:
# IGT, Temperton FFT.
# This features will not work if line below is uncommented.
#CFLAGS += -DNO_FORTRAN

# --Compilers-- Choose one of the following
# gnu - tested for gcc 3.2.3 - 4.3.2
# intel - tested on icc 9.0 - 11.0
# compaq - tested on Compaq C V6.5-303 (dtk) - last tested in 2007
# ibm - tested on xlc 8.0 - last tested in 2008
# hpux - tested on ia64
# other -
COMPILER = gnu

# Additional options for compiler. For instance, one may use -m32 to force 32 bit compilation in 64
# bit environment
CFLAGS +=

# --Warnings-- If uncommented, all warning are suppressed
RELEASE = on

# --FFTW3 paths--
# Specify path to headers and libraries of FFTW3. Some systems do not need them at all, some specify
# special global variables (first or second 2 lines), on some - FFTW3 is installed under user
# account (next 2 lines). Under Windows it may be required to specify paths manually (last 2 lines).
# Relative (to location of this Makefile) paths should be immediately transformed into absolute
# ones using "$(abspath ...)".
#FFTW3_INC_PATH := $(FFTW_INC)
#FFTW3_LIB_PATH := $(FFTW_LIB)
#FFTW3_INC_PATH := $(FFTWINCLUDE)
#FFTW3_LIB_PATH := $(FFTWLIB)
#FFTW3_INC_PATH := $(HOME)/include
#FFTW3_LIB_PATH := $(HOME)/lib
#FFTW3_INC_PATH := "$(abspath ./../lib)"
#FFTW3_LIB_PATH := "$(abspath ./../lib)"

# --C standard--
# ADDA code relies on certain parts of C99 standard, therefore support of this standard is checked
# at compile time (and error is produced if the test fails). However, if you believe that your
# compiler supports all the required features (listed in const.h), but do not define itself
# conforming to C99, you may uncomment the following option to override the test. Do it at your own
# risk!
#CFLAGS += -DOVERRIDE_STDC_TEST

#===================================================================================================
# !!! End of control section. Everything below is not designed to be modified by user. However,
# advanced users may wish to modify some compilers  flags below, especially when using 'other'
# compiler.

# Unconditional variables
#===================================================================================================

# CSOURCE contains all C files in source folder, but we specify them explicitly to avoid possible
# problems if other files somehow appear in this folder.
CSOURCE = ADDAmain.c  calculator.c make_particle.c CalculateE.c GenerateB.c debug.c timing.c \
          linalg.c iterative.c memory.c comm.c matvec.c param.c fft.c crosssec.c Romberg.c io.c \
          prec_time.c vars.c mt19937ar.c sinint.c
# Some other fortran files may be added below
FSOURCE = d07hre.f d09hre.f d113re.f d132re.f dadhre.f dchhre.f dcuhre.f dfshre.f dinhre.f \
          drlhre.f dtrhre.f propaesplibreintadda.f

COBJECTS := $(CSOURCE:.c=.o)
CDEPEND  := $(CSOURCE:.c=.d)
# FOBJECTS is defined after all additions to FSOURCE

# Path to search for source files from the child Makefiles (in child folders)
# Those are used for vpath directives in child Makefiles
PARENT   = ..
FFOLDER  = fort
CPATH    := $(PARENT)
HPATH    := $(PARENT)
FPATH    := $(PARENT)/$(FFOLDER)

LDLIBS   = -lm
DEPFLAG  = -M
DFFLAG   = -MF # this is required
MFILE    = Makefile # this file; child makefiles are named the same
MFILES   := $(MFILE) ../$(MFILE)
# Fortran sources generate a lot of warnings
FWARN    = -w

#===================================================================================================
# Conditional variables that depend on the values of optional flags.
#===================================================================================================

# adds link to fftw3 library and modifies paths when above is commented
ifeq ($(findstring -DFFT_TEMPERTON,$(CFLAGS)),)
  $(info Compiling with FFTW3)
  LDLIBS += -lfftw3
  ifneq ($(FFTW3_INC_PATH),)
    CFLAGS += -I$(FFTW3_INC_PATH)
  endif
  ifneq ($(FFTW3_LIB_PATH),)
    LDFLAGS += -L$(FFTW3_LIB_PATH)
  endif
else
  $(info Compiling with Temperton FFT)
  ifeq ($(findstring -DNO_FORTRAN,$(CFLAGS)),)
    FSOURCE += cfft99D.f
  else
    $(error Temperton FFT (FFT_TEMPERTON) relies on Fortran sources, hence it is incompatible with \
    NO_FORTRAN)
  endif
endif
# essentials of NO_FORTRAN option
ifneq ($(findstring -DNO_FORTRAN,$(CFLAGS)),)
  $(info Compiling without Fortran sources)
  FSOURCE = 
endif
# compiler warnings and optimization flags;
# when release, warning are turned off below
FOPT = $(COPT)
# This is for additional libraries that may be needed when using C linker on a Fortran sources
# (using Fortran compiler may also cause some problems, i.e. for MPI mode). Particular values are
# assigned for each compiler below
FLIBS +=
ifeq ($(COMPILER),gnu)
  # You can add '-march=...' or 'mcpu=...' to COPT and FOPT in this section
  CC    = gcc
  CF    = g77
  COPT  = -std=c99 -O3 -ffast-math -funroll-loops
  FOPT  = -O
  FLIBS += -lg2c
  CWARN = -Wall -W -pedantic -Wcast-qual -Wpointer-arith \
          -Wwrite-strings -Wstrict-prototypes -Wno-uninitialized \
          -Wno-unknown-pragmas -Wno-unused-parameter
else ifeq ($(COMPILER),intel)
  CC    = icc
  CF    = ifort
  COPT  = -std=c99 -O3
  FOPT  = -O3
  FLIBS += -lifcore
  CWARN = -Wall -Wno-uninitialized -wd810,869,981,1418,1419,1572,2259
  # if IPO is used, corresponding flags should be added to linker options: LDFLAGS += ...
else ifeq ($(COMPILER),compaq)
  # This compiler was not tested since 2007. In particular, warning options may not fit exactly the
  # C99 standard, to which the code was transferred. Its support for 64 bit compilations is also
  # undefined. If you happen to use this compiler, please report results to the authors.
  #
  # You can add option '-arch host' to COPT and FOPT in this section
  CC    = cc
  CF    = f77
  COPT  = -fast
  CWARN = -w0 -msg_disable nestedcomment,unknownpragma,unreachcode
else ifeq ($(COMPILER),ibm)
  # This compiler was not tested since 2008. In particular, it is not clear, whether and what FLIBS
  # should be used. If you happen to use this compiler, please report results to the authors.
  #
  # -O5 implies "-arch=auto", which tunes compilation exclusively for the host
  # machine. However, it will not work in some configurations.
  # Then use '-O3 -qarch=... -qtune=...' instead
  CC      = xlc
  CF      = xlf
  COPT    = -O3 -qcache=auto -qlanglvl=extc99
  # DFREDIRECT = 2>nul -qipa=level=2 -qhot
  DEPFLAG = -qmakedep=gcc -qsyntaxonly
  CWARN   = -qsuppress=1506-224:1506-342:1500-036
else ifeq ($(COMPILER),hpux)
  CC    = cc
  CF    = f90
  FOPT  = +O3  +DD64
  COPT  =  -AC99 +O3  +DD64
  CWARN =
  CFLAGS += -DNOT_USE_LOCK
else ifeq ($(COMPILER),other)

else
  $(error Unknown compiler set '$(COMPILER)')
endif
$(info Using compiler set '$(COMPILER)')

# if Fortran sources are used, Fortran libraries are added
ifneq ($(strip $(FSOURCE)),)
  LDLIBS += $(FLIBS)
endif
# if 'release' turn off warningns
ifdef RELEASE
  $(info All warnings are suppressed)
  CWARN = -w
  LDFLAGS += -w
endif
# Finalize option flags (almost)
CFLAGS += $(COPT) $(CWARN)
FFLAGS += $(FOPT) $(FWARN)

FOBJECTS = $(FSOURCE:.f=.o)

#===================================================================================================
#end of check for targets implying compilation
endif
# Main action part
#===================================================================================================

.EXPORT_ALL_VARIABLES:
.PHONY: seq mpi all clean cleanseq cleanmpi

all: seq mpi

seq:
	@echo "Compiling sequential version of ADDA"
	$(MAKE) -C $(SEQ)

mpi:
	@echo "Compiling MPI version of ADDA"
	$(MAKE) -C $(MPI)

clean: cleanseq cleanmpi

# Moving 'clean' commands to child makefiles is problematic due to included dependency files
cleanseq:
	@echo "Removing sequential compiled files"
	cd $(SEQ) && rm -f *.o *.d $(PROGSEQ) $(PROGSEQ).exe

cleanmpi:
	@echo "Removing MPI compiled files"
	cd $(MPI) && rm -f *.o *.d $(PROGMPI) $(PROGMPI).exe

# Main makefile for ADDA package
# Requires GNU make to execute. Actual compiling goes in folders 'seq' and 'mpi' for sequential and
# parallel (MPI) version respectively
# $Date::                            $
#
# Copyright (C) 2006-2011 ADDA contributors
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
# Values specified below are added to the values given by the environment.
# These variables can also be given in the command line to make, but it will override all values
# specified in this Makefile. This is especially critical for CFLAGS which should be overridden only
# with much care, including all relevant flags in the new value.
#===================================================================================================

#===================================================================================================
# Fixed definitions for non-compilation targets
#===================================================================================================

SHELL      := /bin/sh
# Separate directories to store Makefiles and .o and .d files for sequential and MPI executables
SEQ        := seq
MPI        := mpi
OCL        := ocl
# Names of executables
PROGSEQ    := adda
PROGMPI    := adda_mpi
PROGOCL     := adda_ocl
# Targets implying compilation (empty one is always included)
NONTRIVIAL := all seq mpi ocl
# Files to store last-used command lines for linker, C and Fortran compilers
LDOPTSFILE := .ldopts
COPTSFILE  := .copts
FOPTSFILE  := .fopts
# Common makefile to be included from the specific makefiles for seq and mpi version
COMMONMK   := ../common.mk

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
ifneq ($(if $(MAKECMDGOALS),$(if $(filter $(NONTRIVIAL),$(MAKECMDGOALS)),1,),1),)
#===================================================================================================

# ---Optional flags---
# Uncomment those that you find appropriate; additional information is in specified headers.
# One can also export environmental variable OPTIONS before executing make, or define its value in
# the command line of make. For instance, 'make OPTIONS=DEBUG ...' or 'make OPTIONS+=DEBUG ...'. If
# several options need to be given, they should be given as one argument in quotes with its parts
# separated by spaces, e.g. 'make OPTIONS="DEBUG FFT_TEMPERTON" ...'. OPTIONS that are uncommented
# below are appended to the list specified elsewhere.
# Full list of possible options is the following:
VALID_OPTS := DEBUG DEBUGFULL FFT_TEMPERTON PRECISE_TIMING NOT_USE_LOCK ONLY_LOCKFILE NO_FORTRAN \
              OVERRIDE_STDC_TEST

# Debug mode. By default, release configuration is used (no debug, no warnings, maximum
# optimization). DEBUG turns on producing debugging symbols (-g) and warnings and brings
# optimization down to O2 (this is required to produce all possible warnings by the compiler).
# DEBUGFULL turns off optimization completely (for more accurate debugging symbols) and turns on
# additional diagnostic messages in the source code (debug.c/h).
#override OPTIONS += DEBUG
#override OPTIONS += DEBUGFULL

# Temperton FFT (fft.h).
#override OPTIONS += FFT_TEMPERTON

# Precise timing (prec_timing.h).
#override OPTIONS += PRECISE_TIMING

# Controls the mode of file locking, if any (io.h). Use at maximum one of the following options.
#override OPTIONS += NOT_USE_LOCK
#override OPTIONS += ONLY_LOCKFILE

# This should be uncommented, if compiling of Fortran sources cause problems. However, currently a
# number of ADDA (optional) features rely on Fortran sources, in particular: IGT, Temperton FFT.
# This features will not work if line below is uncommented.
#override OPTIONS += NO_FORTRAN

# ADDA code relies on certain parts of C99 standard, therefore support of this standard is checked
# at compile time (and error is produced if the test fails). However, if you believe that your
# compiler supports all the required features (listed in const.h), but do not define itself
# conforming to C99, you may uncomment the following option to override the test. Do it at your own
# risk!
#override OPTIONS += OVERRIDE_STDC_TEST

# ---Compilers---
# Choose one of the following. Can also be specified from command line to make (see explanation
# above for OPTIONS), overriding definition below.
# gnu - tested for gcc 3.2.3 - 4.3.2
# intel - tested on icc 9.0 - 11.0
# compaq - tested on Compaq C V6.5-303 (dtk) - last tested in 2007
# ibm - tested on xlc 8.0 - last tested in 2008
# hpux - tested on ia64
# other -
COMPILER := gnu

# Additional options for compiler. For instance, one may use -m32 to force 32 bit compilation in 64
# bit environment. Flags specified below are appended to the ones specified in the environment or
# command line of make (see explanation above for OPTIONS).
override EXTRA_CFLAGS +=

# --FFTW3 paths--
# Specify path to headers and libraries of FFTW3. Some systems do not need them at all, some specify
# special global variables (first or second 2 lines), on some - FFTW3 is installed under user
# account (next 2 lines). Under Windows it may be required to specify paths manually (last 2 lines).
# Relative (to location of this Makefile) paths should be immediately transformed into absolute
# ones using "$(abspath ...)". This two variables (FFTW3_INC_PATH and FFTW3_LIB_PATH) can also be
# defined in the enviroment or in the command line of make (see explanation above for OPTIONS).
#FFTW3_INC_PATH := $(FFTW_INC)
#FFTW3_LIB_PATH := $(FFTW_LIB)
#FFTW3_INC_PATH := $(FFTWINCLUDE)
#FFTW3_LIB_PATH := $(FFTWLIB)
#FFTW3_INC_PATH := $(HOME)/include
#FFTW3_LIB_PATH := $(HOME)/lib
#FFTW3_INC_PATH := "$(abspath ./../lib)"
#FFTW3_LIB_PATH := "$(abspath ./../lib)"

#===================================================================================================
# !!! End of control section. Everything below is not designed to be modified by user. However,
# advanced users may wish to modify some compilers  flags below, especially when using 'other'
# compiler.

# Unconditional variables
#===================================================================================================

# CSOURCE contains all C files in source folder, but we specify them explicitly to avoid possible
# problems if other files somehow appear in this folder.
CSOURCE := ADDAmain.c  calculator.c make_particle.c CalculateE.c GenerateB.c debug.c timing.c \
           linalg.c iterative.c memory.c comm.c matvec.c param.c fft.c crosssec.c Romberg.c io.c \
           prec_time.c vars.c mt19937ar.c sinint.c
# Fortran files are all that are present in src/fort folder, but given explicitly
FSOURCE := d07hre.f d09hre.f d113re.f d132re.f dadhre.f dchhre.f dcuhre.f dfshre.f dinhre.f \
           drlhre.f dtrhre.f propaesplibreintadda.f

COBJECTS := $(CSOURCE:.c=.o)
CDEPEND  := $(CSOURCE:.c=.d)
# FOBJECTS is defined after all additions to FSOURCE

# Path to search for source files from the child Makefiles (in child folders)
# Those are used for vpath directives in child Makefiles
PARENT   := ..
FFOLDER  := fort
CPPFOLDER:= cpp
CPATH    := $(PARENT)
HPATH    := $(PARENT)
FPATH    := $(PARENT)/$(FFOLDER)
CPPPATH  := $(PARENT)/$(CPPFOLDER)

LDLIBS   := -lm
DEPFLAG  := -M
DFFLAG   := -MF # this is required
# Fortran sources generate a lot of warnings, we do not plan to investigate them
FWARN    := -w

#===================================================================================================
# Conditional variables that depend on the values of optional flags.
#===================================================================================================
#
# Process OPTIONS
# Produce warnings if any of OPTIONS are not present in VALID_OPTS list
$(foreach var,$(filter-out $(VALID_OPTS),$(OPTIONS)),\
              $(warning The option '$(var)' is not recognized. Please check the spelling.))

$(info --- Compilation options: ---)
ifneq ($(filter DEBUGFULL,$(OPTIONS)),)
  $(info Full debug mode)
  DBGLVL := 2
  CDEFS  += -DDEBUGFULL
  ifneq ($(filter DEBUG,$(OPTIONS)),)
    $(warning DEBUG has no effect when DEBUGFULL is enabled)
  endif
else ifneq ($(filter DEBUG,$(OPTIONS)),)
  $(info Debug mode)
  DBGLVL := 1
  CDEFS  += -DDEBUG
else
  $(info Release mode)
  DBGLVL := 0
endif
ifneq ($(filter FFT_TEMPERTON,$(OPTIONS)),)
  $(info Temperton FFT)
  CDEFS += -DFFT_TEMPERTON
  ifeq ($(filter NO_FORTRAN,$(OPTIONS)),)
    FSOURCE += cfft99D.f
  else
    $(error Temperton FFT (FFT_TEMPERTON) relies on Fortran sources, hence it is incompatible with \
    NO_FORTRAN)
  endif
else
  $(info FFTW3)
  LDLIBS += -lfftw3
  ifdef FFTW3_INC_PATH
    CFLAGS += -I$(FFTW3_INC_PATH)
  endif
  ifdef FFTW3_LIB_PATH
    LDFLAGS += -L$(FFTW3_LIB_PATH)
  endif
endif
ifneq ($(filter PRECISE_TIMING,$(OPTIONS)),)
  $(info Precise timing)
  CDEFS += -DPRECISE_TIMING
endif
ifneq ($(filter NOT_USE_LOCK,$(OPTIONS)),)
  $(info No locks at all)
  CDEFS += -DNOT_USE_LOCK
  ifneq ($(filter ONLY_LOCKFILE,$(OPTIONS)),)
    $(warning ONLY_LOCKFILE has no effect when NOT_USE_LOCK is enabled)
  endif
else ifneq ($(filter ONLY_LOCKFILE,$(OPTIONS)),)
  $(info Only lock files (without system locks))
  CDEFS += -DDONLY_LOCKFILE
endif
ifneq ($(filter NO_FORTRAN,$(OPTIONS)),)
  $(info Without Fortran sources)
  CDEFS += -DNO_FORTRAN
  FSOURCE =
endif
ifneq ($(filter OVERRIDE_STDC_TEST,$(OPTIONS)),)
  $(info Overriding test for C99 conformance)
  CDEFS += -DOVERRIDE_STDC_TEST
endif
# Process EXTRA_CFLAGS
ifneq ($(strip $(EXTRA_CFLAGS)),)
  $(info Extra CFLAGS: '$(EXTRA_CFLAGS)')
  CFLAGS += $(EXTRA_CFLAGS)
endif

# By default Fortran optimization options are the same as C, but specific options can be used for
# some compilers below.
FOPT = $(COPT)
# This is for additional libraries that may be needed when using C linker on a Fortran sources
# (using Fortran compiler may also cause some problems, i.e. for MPI mode). Particular values are
# assigned for each compiler below.
FLIBS +=
ifeq ($(COMPILER),gnu)
  # You can add '-march=...' or 'mcpu=...' to COPT and FOPT in this section
  CC    := gcc
  CCPP  := g++
  CF    := g77
  CSTD  := -std=c99
  CDBG  := -g
  COPT1 := -O2
  COPT2 := -O3 -ffast-math -funroll-loops
  CWARN := -Wall -Wextra -pedantic -Wcast-qual -Wpointer-arith -Wwrite-strings -Wstrict-prototypes \
           -Wstrict-aliasing=1 -Wshadow -Wcast-align -Wnested-externs -Wcomment -Wno-unknown-pragmas
  FOPT  := -O
  FLIBS += -lg2c

else ifeq ($(COMPILER),intel)
  CC    := icc
  CF    := ifort
  CSTD  := -std=c99 -vec-report0 # the last flag is used to always remove vectorization remarks
  CDBG  := -g
  COPT1 := -O2
  COPT2 := -O3
  CWARN := -Wall -Wcheck -diag-disable 981,1418,1419,1572,2259
  FWARN += -vec-report0
  FLIBS += -lifcore
  # if IPO is used, corresponding flags should be added to linker options: LDFLAGS += ...
else ifeq ($(COMPILER),compaq)
  # This compiler was not tested since 2007. In particular, warning options may not fit exactly the
  # C99 standard, to which the code was transferred. Its support for 64 bit compilations is also
  # undefined. If you happen to use this compiler, please report results to the authors.
  #
  # You can add option '-arch host' to COPT and FOPT in this section
  CC    := cc
  CF    := f77
  #CSTD  = -std=c99
  CDBG  := -g
  COPT1 := -O2
  COPT2 := -fast
  CWARN := -w0 -msg_disable nestedcomment,unknownpragma,unreachcode
else ifeq ($(COMPILER),ibm)
  # This compiler was not tested since 2008. In particular, it is not clear, whether and what FLIBS
  # should be used. If you happen to use this compiler, please report results to the authors.
  #
  # -O5 implies "-arch=auto", which tunes compilation exclusively for the host
  # machine. However, it will not work in some configurations.
  # Then use '-O3 -qarch=... -qtune=...' instead
  CC    := xlc
  CF    := xlf
  CSTD  := -qlanglvl=extc99
  CDBG  := -g
  COPT1 := -O2
  COPT2 := -O3 -qcache=auto
  # DFREDIRECT = 2>nul -qipa=level=2 -qhot
  DEPFLAG := -qmakedep=gcc -qsyntaxonly
  CWARN   := -qsuppress=1506-224:1506-342:1500-036
else ifeq ($(COMPILER),hpux)
  CC    := cc
  CF    := f90
  CSTD  := -AC99
  CDBG  := -g
  COPT1 := +O2 +DD64
  COPT2 := +O3 +DD64
  CWARN :=
  CFLAGS += -DNOT_USE_LOCK
else ifeq ($(COMPILER),other)
# add here definitions corresponding to 'other' compiler, if you want to use it.
else
  $(error Unknown compiler set '$(COMPILER)')
endif
$(info Compiler set '$(COMPILER)')

# if Fortran sources are used, Fortran libraries are added
ifneq ($(strip $(FSOURCE)),)
  LDLIBS += $(FLIBS)
endif
# if 'release' turn off warningns
ifeq ($(DBGLVL),0)
  CWARN   := -w
  LDFLAGS += -w
  COPT    := $(COPT2)
  CDBG    :=
else ifeq ($(DBGLVL),1)
  COPT := $(COPT1)
else ifeq ($(DBGLVL),2)
  COPT :=
endif
# Finalize option flags; these aggregates are used as a whole furher on
LDFLAGS += $(LDLIBS)
CFLAGS  += $(CSTD) $(CDBG) $(COPT) $(CWARN) $(CDEFS)
FFLAGS  += $(FOPT) $(FWARN)

FOBJECTS := $(FSOURCE:.f=.o)

$(info ----------------------------)

#===================================================================================================
#end of check for targets implying compilation
endif
# Main action part
#===================================================================================================

.EXPORT_ALL_VARIABLES:
.PHONY: seq mpi all clean cleanseq cleanmpi ocl

all: seq mpi ocl

seq:
	@echo "Compiling sequential version of ADDA"
	$(MAKE) -C $(SEQ)

mpi:
	@echo "Compiling MPI version of ADDA"
	$(MAKE) -C $(MPI)

ocl:
	@echo "Compiling OpenCL version of ADDA"
	$(MAKE) -C $(OCL)


clean: cleanseq cleanmpi cleanocl

# 'clean' commands are executed here to keep them simple. Child makefiles are used only for
# compilation and thus contain quite heavy processing.
cleanseq:
	@echo "Removing sequential compiled files"
	cd $(SEQ) && rm -f *.o *.d $(LDOPTSFILE) $(COPTSFILE) $(FOPTSFILE) $(PROGSEQ) $(PROGSEQ).exe

cleanmpi:
	@echo "Removing MPI compiled files"
	cd $(MPI) && rm -f *.o *.d $(LDOPTSFILE) $(COPTSFILE) $(FOPTSFILE) $(PROGMPI) $(PROGMPI).exe

cleanocl:
	@echo "Removing OpenCL compiled files"
	cd $(OCL) && rm -f *.o *.d $(LDOPTSFILE) $(COPTSFILE) $(FOPTSFILE) $(PROGOCL) $(PROGOCL).exe
